import Foundation
#if canImport(Accelerate)
import Accelerate
#endif

public enum VivoQMError: Error, Sendable, CustomStringConvertible {
    case invalid(String), unsupported(String), capacity(String), convergence(String)
    public var description: String {
        switch self {
        case .invalid(let s): return "Invalid electronic-structure input: \(s)"
        case .unsupported(let s): return "Unsupported electronic-structure method: \(s)"
        case .capacity(let s): return "Electronic-structure resource limit: \(s)"
        case .convergence(let s): return "Electronic-structure convergence failure: \(s)"
        }
    }
}

/// FP64, row-major storage. Decoded data is revalidated at every numerical entry.
public struct VivoQMMatrix: Codable, Sendable, Equatable {
    public let rows: Int
    public let columns: Int
    public private(set) var values: [Double]

    public init(rows: Int, columns: Int, values: [Double]) throws {
        self.rows = rows; self.columns = columns; self.values = values
        try validate()
    }
    private enum CodingKeys: String, CodingKey { case rows, columns, values }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy:CodingKeys.self)
        rows = try c.decode(Int.self,forKey:.rows)
        columns = try c.decode(Int.self,forKey:.columns)
        values = try c.decode([Double].self,forKey:.values)
        try validate()
    }
    internal init(_ rows: Int, _ columns: Int, repeating value: Double = 0) {
        precondition(rows > 0 && columns > 0 && rows <= 4096 && columns <= 4096)
        self.rows = rows; self.columns = columns
        values = Array(repeating: value, count: rows * columns)
    }
    public func validate() throws {
        guard (1...4096).contains(rows), (1...4096).contains(columns),
              values.count == rows * columns, values.allSatisfy(\.isFinite) else {
            throw VivoQMError.invalid("matrix dimensions, storage, or finite values")
        }
    }
    public subscript(_ row: Int, _ column: Int) -> Double {
        get { values[row * columns + column] }
        set { values[row * columns + column] = newValue }
    }
    public static func identity(_ n: Int) throws -> Self {
        guard (1...4096).contains(n) else { throw VivoQMError.capacity("identity dimension") }
        var result = Self(n, n)
        for i in 0..<n { result[i, i] = 1 }
        return result
    }
    public var transposed: Self {
        var result = Self(columns, rows)
        for i in 0..<rows { for j in 0..<columns { result[j, i] = self[i, j] } }
        return result
    }
    public var frobeniusNorm: Double { sqrt(values.reduce(0) { $0 + $1 * $1 }) }
    public func scaled(_ scale: Double) throws -> Self {
        try Self(rows: rows, columns: columns, values: values.map { $0 * scale })
    }
    public func adding(_ other: Self, scale: Double = 1) throws -> Self {
        try validate(); try other.validate()
        guard rows == other.rows, columns == other.columns, scale.isFinite else {
            throw VivoQMError.invalid("matrix addition shapes")
        }
        return try Self(rows: rows, columns: columns,
                        values: zip(values, other.values).map { $0 + scale * $1 })
    }
    public func multiplied(by other: Self) throws -> Self {
        try validate(); try other.validate()
        guard columns == other.rows else { throw VivoQMError.invalid("matrix product shapes") }
        var result = Self(rows, other.columns)
        #if canImport(Accelerate)
        values.withUnsafeBufferPointer { a in
            other.values.withUnsafeBufferPointer { b in
                result.values.withUnsafeMutableBufferPointer { c in
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                                Int32(rows), Int32(other.columns), Int32(columns),
                                1, a.baseAddress!, Int32(columns), b.baseAddress!,
                                Int32(other.columns), 0, c.baseAddress!, Int32(other.columns))
                }
            }
        }
        #else
        for i in 0..<rows { for k in 0..<columns {
            let a = self[i, k]
            for j in 0..<other.columns { result[i, j] += a * other[k, j] }
        } }
        #endif
        try result.validate()
        return result
    }
    /// Returns C^T A C; C may be rectangular.
    public func transformed(by c: Self) throws -> Self {
        try c.transposed.multiplied(by: multiplied(by: c))
    }
    public func requireSymmetric(tolerance: Double = 1e-10) throws {
        try validate()
        guard rows == columns, tolerance.isFinite, tolerance > 0 else { throw VivoQMError.invalid("non-square symmetric matrix or invalid symmetry tolerance") }
        let scale = max(1, values.map(abs).max() ?? 0)
        for i in 0..<rows { for j in 0..<i {
            guard abs(self[i, j] - self[j, i]) <= tolerance * scale else {
                throw VivoQMError.invalid("matrix is not symmetric")
            }
        } }
    }

    /// Bounded cyclic Jacobi reference eigensolver. Columns are eigenvectors.
    /// Dense contractions use Accelerate on Apple; this FP64 reference never
    /// silently downcasts a secular equation to Metal FP32.
    public func symmetricEigen(tolerance: Double = 1e-12,
                               maximumSweeps: Int = 100) throws -> VivoQMEigensystem {
        try requireSymmetric()
        guard rows <= 512, tolerance.isFinite, tolerance > 0,
              (1...1000).contains(maximumSweeps) else {
            throw VivoQMError.capacity("Jacobi supports at most 512 rows and 1000 sweeps")
        }
        let n = rows
        var a = self
        var v = try Self.identity(n)
        let threshold = tolerance * max(1, frobeniusNorm)
        for _ in 0..<maximumSweeps {
            var maximum = 0.0
            for p in 0..<n { for q in (p + 1)..<n {
                let apq = a[p, q]
                maximum = max(maximum, abs(apq))
                if abs(apq) <= threshold { continue }
                let tau = (a[q, q] - a[p, p]) / (2 * apq)
                let t = (tau >= 0 ? 1.0 : -1.0) / (abs(tau) + hypot(1, tau))
                let c = 1 / sqrt(1 + t * t), s = t * c
                let app = a[p, p], aqq = a[q, q]
                a[p, p] = app - t * apq; a[q, q] = aqq + t * apq
                a[p, q] = 0; a[q, p] = 0
                for k in 0..<n where k != p && k != q {
                    let kp = a[k, p], kq = a[k, q]
                    a[k, p] = c * kp - s * kq; a[p, k] = a[k, p]
                    a[k, q] = s * kp + c * kq; a[q, k] = a[k, q]
                }
                for k in 0..<n {
                    let kp = v[k, p], kq = v[k, q]
                    v[k, p] = c * kp - s * kq; v[k, q] = s * kp + c * kq
                }
            } }
            if maximum <= threshold {
                let order = (0..<n).sorted {
                    a[$0, $0] == a[$1, $1] ? $0 < $1 : a[$0, $0] < a[$1, $1]
                }
                var sorted = Self(n, n)
                for (j, original) in order.enumerated() {
                    // Fix the arbitrary sign deterministically.
                    let pivot = (0..<n).max { abs(v[$0, original]) < abs(v[$1, original]) }!
                    let phase = v[pivot, original] < 0 ? -1.0 : 1.0
                    for i in 0..<n { sorted[i, j] = phase * v[i, original] }
                }
                return .init(values: order.map { a[$0, $0] }, vectors: sorted)
            }
        }
        throw VivoQMError.convergence("symmetric eigensolver exhausted its sweep budget")
    }

    internal func solve(_ rhs: [Double], pivotTolerance: Double = 1e-14) throws -> [Double] {
        try validate()
        guard rows == columns, rhs.count == rows, rhs.allSatisfy(\.isFinite) else {
            throw VivoQMError.invalid("linear solve shapes")
        }
        let n = rows
        var a = self, b = rhs
        let scale = max(1, values.map(abs).max() ?? 0)
        for k in 0..<n {
            let pivot = (k..<n).max { abs(a[$0, k]) < abs(a[$1, k]) }!
            guard abs(a[pivot, k]) > pivotTolerance * scale else {
                throw VivoQMError.convergence("singular linear system")
            }
            if pivot != k {
                for j in 0..<n { let t = a[k, j]; a[k, j] = a[pivot, j]; a[pivot, j] = t }
                b.swapAt(k, pivot)
            }
            for i in (k + 1)..<n {
                let factor = a[i, k] / a[k, k]
                for j in k..<n { a[i, j] -= factor * a[k, j] }
                b[i] -= factor * b[k]
            }
        }
        var x = Array(repeating: 0.0, count: n)
        for i in stride(from: n - 1, through: 0, by: -1) {
            var sum = b[i]
            for j in (i + 1)..<n { sum -= a[i, j] * x[j] }
            x[i] = sum / a[i, i]
        }
        guard x.allSatisfy(\.isFinite) else { throw VivoQMError.convergence("non-finite linear solve") }
        return x
    }
}

public struct VivoQMEigensystem: Sendable {
    public let values: [Double]
    public let vectors: VivoQMMatrix
}

public enum VivoQMUnits {
    public static let bohrInNM = 0.0529177210544
    public static let hartreeInKJPerMol = 2625.4996394799
}
