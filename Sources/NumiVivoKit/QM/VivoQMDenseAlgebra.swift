import Foundation
#if canImport(Accelerate)
import Accelerate
#endif

/// Real row-major FP64 matrices. Decoding validates shape and finiteness.
public struct VivoQMMatrix: Codable, Sendable, Equatable {
    public let rows: Int
    public let columns: Int
    public internal(set) var values: [Double]
    public init(rows: Int, columns: Int, values: [Double]) throws {
        let (count, overflow) = rows.multipliedReportingOverflow(by: columns)
        guard rows >= 0, columns >= 0, !overflow, values.count == count, values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("matrix shape or nonfinite values")
        }
        self.rows = rows; self.columns = columns; self.values = values
    }
    init(_ rows: Int, _ columns: Int, repeating: Double = 0) {
        self.rows = rows; self.columns = columns; values = .init(repeating: repeating, count: rows * columns)
    }
    enum CodingKeys: String, CodingKey { case rows, columns, values }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(rows: c.decode(Int.self, forKey: .rows), columns: c.decode(Int.self, forKey: .columns),
                      values: c.decode([Double].self, forKey: .values))
    }
    public subscript(_ r: Int, _ c: Int) -> Double {
        get { values[r * columns + c] }
        set { values[r * columns + c] = newValue }
    }
    public static func identity(_ n: Int) throws -> Self {
        guard n > 0, n <= 8192 else { throw VivoChemistryError.resourceLimit("identity dimension") }
        var result = Self(n,n); for i in 0..<n { result[i,i] = 1 }; return result
    }
    public var transposed: Self {
        var t = Self(columns,rows)
        for i in 0..<rows { for j in 0..<columns { t[j,i] = self[i,j] } }; return t
    }
    public var frobeniusNorm: Double { sqrt(values.reduce(0) { $0 + $1 * $1 }) }
    public func multiplied(by b: Self) throws -> Self {
        guard columns == b.rows, values.allSatisfy(\.isFinite), b.values.allSatisfy(\.isFinite),
              rows <= Int(Int32.max), columns <= Int(Int32.max), b.columns <= Int(Int32.max),
              rows <= 67_108_864 / max(1,b.columns) else {
            throw VivoChemistryError.invalid("matrix product shape or values")
        }
        var out = Self(rows,b.columns)
        if rows == 0 || columns == 0 || b.columns == 0 { return out }
        #if canImport(Accelerate)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    Int32(rows), Int32(b.columns), Int32(columns), 1,
                    values, Int32(columns), b.values, Int32(b.columns), 0, &out.values, Int32(b.columns))
        #else
        for i in 0..<rows { for k in 0..<columns {
            let a = self[i,k]
            for j in 0..<b.columns { out[i,j] += a * b[k,j] }
        } }
        #endif
        guard out.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("matrix product overflow") }
        return out
    }
    func scaled(_ factor: Double) -> Self { var o = self; for i in o.values.indices { o.values[i] *= factor }; return o }
    func adding(_ b: Self, scale: Double = 1) throws -> Self {
        guard rows == b.rows, columns == b.columns else { throw VivoChemistryError.invalid("matrix addition shape") }
        var out = self; for i in values.indices { out.values[i] += scale * b.values[i] }; return out
    }
    public func congruence(_ c: Self) throws -> Self { try c.transposed.multiplied(by: self).multiplied(by: c) }
}

public enum VivoQMDenseAlgebra {
    /// Cyclic Jacobi reference eigensolver, explicitly CPU FP64 on all platforms.
    public static func symmetricEigen(_ input: VivoQMMatrix, tolerance: Double = 1e-12,
                                      maximumSweeps: Int = 100) throws -> (values: [Double], vectors: VivoQMMatrix) {
        let n = input.rows
        guard n > 0, n == input.columns, n <= 8192, tolerance.isFinite, tolerance > 0,
              maximumSweeps > 0, input.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("eigensolver input") }
        let scale = max(1, input.values.map(abs).max() ?? 1)
        for i in 0..<n { for j in 0..<i where abs(input[i,j] - input[j,i]) > tolerance * scale * 10 {
            throw VivoChemistryError.invalid("matrix is not symmetric")
        } }
        var a = input, v = try VivoQMMatrix.identity(n)
        var converged = n == 1
        for _ in 0..<maximumSweeps {
            var largest = 0.0
            for p in 0..<n { for q in (p+1)..<n {
                let apq = a[p,q]; largest = max(largest, abs(apq))
                if abs(apq) <= tolerance * scale { continue }
                let tau = (a[q,q] - a[p,p]) / (2 * apq)
                let t = (tau >= 0 ? 1.0 : -1.0) / (abs(tau) + hypot(1,tau))
                let c = 1 / sqrt(1 + t*t), s = t*c
                a[p,p] -= t * apq; a[q,q] += t * apq; a[p,q] = 0; a[q,p] = 0
                for k in 0..<n where k != p && k != q {
                    let x = a[k,p], y = a[k,q]
                    a[k,p] = c*x - s*y; a[p,k] = a[k,p]
                    a[k,q] = s*x + c*y; a[q,k] = a[k,q]
                }
                for k in 0..<n {
                    let x = v[k,p], y = v[k,q]; v[k,p] = c*x-s*y; v[k,q] = s*x+c*y
                }
            } }
            if largest <= tolerance * scale { converged = true; break }
        }
        guard converged, a.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("symmetric eigensolver") }
        let order = (0..<n).sorted { a[$0,$0] < a[$1,$1] }
        var sorted = VivoQMMatrix(n,n)
        for j in 0..<n { for i in 0..<n { sorted[i,j] = v[i,order[j]] } }
        return (order.map { a[$0,$0] }, sorted)
    }
    /// Partial-pivot Gaussian elimination, used for small DIIS systems only.
    static func solve(_ input: VivoQMMatrix, rhs: [Double], pivotTolerance: Double = 1e-14) throws -> [Double] {
        let n = input.rows
        guard n > 0, n == input.columns, rhs.count == n, rhs.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("linear solve") }
        var a = input, b = rhs
        for k in 0..<n {
            var p = k; for i in (k+1)..<n where abs(a[i,k]) > abs(a[p,k]) { p = i }
            guard abs(a[p,k]) > pivotTolerance else { throw VivoChemistryError.convergence("singular linear system") }
            if p != k {
                for j in k..<n { let x=a[k,j]; a[k,j]=a[p,j]; a[p,j]=x }
                b.swapAt(k,p)
            }
            for i in (k+1)..<n {
                let f = a[i,k]/a[k,k]; a[i,k] = 0
                for j in (k+1)..<n { a[i,j] -= f*a[k,j] }; b[i] -= f*b[k]
            }
        }
        var x = [Double](repeating:0,count:n)
        for i in stride(from:n-1,through:0,by:-1) {
            var r = b[i]; for j in (i+1)..<n { r -= a[i,j]*x[j] }; x[i] = r/a[i,i]
        }
        guard x.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("nonfinite linear solution") }; return x
    }
    /// exp(-K), where K is real antisymmetric; no unconstrained orbital update.
    public static func orbitalRotation(generator k: VivoQMMatrix) throws -> VivoQMMatrix {
        let n = k.rows
        guard n > 0, n == k.columns, k.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("orbital generator") }
        for i in 0..<n { for j in 0..<n where abs(k[i,j]+k[j,i]) > 1e-12 {
            throw VivoChemistryError.invalid("orbital generator is not antisymmetric")
        } }
        let norm = k.frobeniusNorm
        guard norm < 1e6 else { throw VivoChemistryError.invalid("orbital generator norm") }
        let squarings = max(0, Int(ceil(log2(max(1,norm / 0.25)))))
        let a = k.scaled(-1 / pow(2,Double(squarings)))
        var sum = try VivoQMMatrix.identity(n), term = sum
        for order in 1...40 {
            term = try term.multiplied(by:a).scaled(1/Double(order)); sum = try sum.adding(term)
            if term.frobeniusNorm < 1e-16 { break }
        }
        for _ in 0..<squarings { sum = try sum.multiplied(by:sum) }
        let defect = try sum.transposed.multiplied(by:sum).adding(.identity(n),scale:-1).frobeniusNorm
        guard defect < 1e-9 else { throw VivoChemistryError.convergence("orbital rotation lost orthogonality") }; return sum
    }
}
