import Foundation

/// Raw Cartesian moments, not traceless quadrupoles: integral rho*x_a*x_b.
/// Units are e, e*Bohr, e*Bohr^2. The six second moments are xx,xy,xz,yy,yz,zz.
public struct VivoCartesianMultipole: Codable, Sendable, Equatable {
    public var positionBohr: VivoVector3D
    public var chargeE: Double
    public var dipoleEBohr: VivoVector3D
    public var secondMomentsEBohr2: [Double]
    public init(positionBohr: VivoVector3D, chargeE: Double, dipoleEBohr: VivoVector3D = .zero,
                secondMomentsEBohr2: [Double] = Array(repeating: 0, count: 6)) {
        self.positionBohr = positionBohr; self.chargeE = chargeE; self.dipoleEBohr = dipoleEBohr
        self.secondMomentsEBohr2 = secondMomentsEBohr2
    }
    public func validate() throws {
        guard positionBohr.isFinite, chargeE.isFinite, dipoleEBohr.isFinite,
              secondMomentsEBohr2.count == 6, secondMomentsEBohr2.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("Cartesian multipole shape or values")
        }
    }
    var moments: [Double] { [chargeE,dipoleEBohr.x,dipoleEBohr.y,dipoleEBohr.z] + secondMomentsEBohr2 }
}
public enum VivoEwaldChargeConvention: String, Codable, Sendable {
    case requireNeutral
    /// Explicit tin-foil, uniform compensating background. Required for charged
    /// intermediate source channels even when their assembled physical cell is neutral.
    case uniformNeutralizingBackground
}
public struct VivoPeriodicElectrostaticConfiguration: Codable, Sendable, Equatable {
    public var alphaPerBohr: Double
    public var realCutoffBohr: Double
    /// Fixed integer reciprocal modes. Cell deformation never silently replans them.
    public var reciprocalHalfWidths: [Int]
    public var chargeConvention: VivoEwaldChargeConvention
    public var neutralToleranceE: Double
    public var maximumImagePairs: Int
    public init(alphaPerBohr: Double = 0.25, realCutoffBohr: Double = 24,
                reciprocalHalfWidths: [Int] = [8,8,8],
                chargeConvention: VivoEwaldChargeConvention = .requireNeutral,
                neutralToleranceE: Double = 1e-8, maximumImagePairs: Int = 1_000_000) {
        self.alphaPerBohr = alphaPerBohr; self.realCutoffBohr = realCutoffBohr
        self.reciprocalHalfWidths = reciprocalHalfWidths; self.chargeConvention = chargeConvention
        self.neutralToleranceE = neutralToleranceE; self.maximumImagePairs = maximumImagePairs
    }
    public func validate() throws {
        guard alphaPerBohr.isFinite, alphaPerBohr > 0, realCutoffBohr.isFinite, realCutoffBohr > 0,
              reciprocalHalfWidths.count == 3, reciprocalHalfWidths.allSatisfy({ (1...128).contains($0) }),
              neutralToleranceE.isFinite, neutralToleranceE > 0, neutralToleranceE <= 1e-3,
              maximumImagePairs > 0 else { throw VivoChemistryError.invalid("Ewald configuration") }
    }
}
public struct VivoMultipolePairScale: Codable, Sendable, Equatable {
    public let first: Int
    public let second: Int
    public let scale: Double
    public init(first: Int, second: Int, scale: Double) { self.first = first; self.second = second; self.scale = scale }
}
public struct VivoPeriodicElectrostaticResult: Codable, Sendable, Equatable {
    public let energyHartree: Double
    public let realEnergyHartree: Double
    public let reciprocalEnergyHartree: Double
    public let selfEnergyHartree: Double
    public let backgroundEnergyHartree: Double
    public let exceptionEnergyHartree: Double
    public let forcesHartreePerBohr: [VivoVector3D]
    /// dE/d(q,mu_x,mu_y,mu_z,Qxx,Qxy,Qxz,Qyy,Qyz,Qzz); these are NOT moments.
    public let momentDerivatives: [[Double]]
    /// dE/d epsilon_ab for r'=(I+epsilon)r and h'=(I+epsilon)h.
    /// Multipoles are held fixed in the laboratory frame. A density/local-frame
    /// caller must add its moment-response chain rule. Not a particle-only torque.
    public let affineStrainDerivativeHartree: VivoQMMatrix
    public let imagePairEvaluations: Int
    public let reciprocalModeCount: Int
}

/// Immutable reciprocal profile. A fixed mesh and active Fourier mode set are
/// retained during cell moves. Mesh errors are separate from electronic residuals.
public struct VivoMultipoleMeshConfiguration: Codable, Sendable, Equatable {
    public let gridDimensions: [Int]
    public let interpolationOrder: Int
    public let numericalPrecision: String
    public init(gridDimensions: [Int] = [64,64,64]) {
        self.gridDimensions = gridDimensions; interpolationOrder = 6
        numericalPrecision = "metal-fp32-variational-mesh-v1"
    }
    public func validate() throws {
        guard gridDimensions.count == 3,gridDimensions.allSatisfy({ $0 >= 16 && $0 <= 512 && ($0 & ($0-1)) == 0 }),
              interpolationOrder == 6,numericalPrecision == "metal-fp32-variational-mesh-v1" else {
            throw VivoChemistryError.invalid("multipolar mesh grid, assignment order or precision")
        }
    }
}
public struct VivoReciprocalElectrostaticResult: Sendable {
    public let energyHartree: Double
    public let forcesHartreePerBohr: [VivoVector3D]
    public let momentDerivatives: [[Double]]
    public let affineStrainDerivativeHartree: VivoQMMatrix
    public let modeCount: Int
}
/// Geometry-local shared mesh resource. The implementation owns and synchronizes
/// its buffers. It may not silently substitute a different reciprocal algorithm.
public struct VivoReciprocalElectrostaticOperator: Sendable {
    public let configuration: VivoMultipoleMeshConfiguration
    let evaluate: @Sendable ([VivoCartesianMultipole],VivoPeriodicCell,VivoPeriodicElectrostaticConfiguration,VivoChemistryBudget) throws -> VivoReciprocalElectrostaticResult
    init(configuration: VivoMultipoleMeshConfiguration,
         evaluate: @escaping @Sendable ([VivoCartesianMultipole],VivoPeriodicCell,VivoPeriodicElectrostaticConfiguration,VivoChemistryBudget) throws -> VivoReciprocalElectrostaticResult) {
        self.configuration = configuration;self.evaluate = evaluate
    }
}

/// Complete bounded FP64 Ewald reference through raw Cartesian quadrupoles.
/// Analytic spatial derivatives through rank five come from the same radial
/// kernel as the energy. This is a numerical reference, not a PME performance claim.
public enum VivoPeriodicElectrostatics {
    static let powers = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[2,0,0],[1,1,0],[1,0,1],[0,2,0],[0,1,1],[0,0,2]]
    static let factors = [1.0,1,1,1,0.5,1,1,0.5,1,0.5]
    static let signs = [1.0,-1,-1,-1,1,1,1,1,1,1]
    struct Pair {
        var energy: Double
        var gradient: VivoVector3D
        var left: [Double]
        var right: [Double]
    }
    static func pair(_ a: [Double], _ b: [Double], tensor: VivoEwaldRadialTensor) -> Pair {
        var result = Pair(energy: 0, gradient: .zero, left: Array(repeating: 0,count: 10), right: Array(repeating: 0,count: 10))
        for i in 0..<10 { for j in 0..<10 {
            let x = powers[i][0]+powers[j][0], y = powers[i][1]+powers[j][1], z = powers[i][2]+powers[j][2]
            let weight = factors[i]*factors[j]*signs[j]
            let v = weight*tensor[x,y,z]
            result.energy += a[i]*b[j]*v
            result.left[i] += b[j]*v; result.right[j] += a[i]*v
            let product = weight*a[i]*b[j]
            if product != 0 {
                result.gradient = result.gradient + VivoVector3D(tensor[x+1,y,z],tensor[x,y+1,z],tensor[x,y,z+1])*product
            }
        } }
        return result
    }
    public static func evaluate(sources: [VivoCartesianMultipole], cell: VivoPeriodicCell,
                                configuration cfg: VivoPeriodicElectrostaticConfiguration = .init(),
                                primaryPairScales: [VivoMultipolePairScale] = [],
                                reciprocalOperator: VivoReciprocalElectrostaticOperator? = nil,
                                budget: VivoChemistryBudget = .init()) throws -> VivoPeriodicElectrostaticResult {
        try cfg.validate(); try budget.validate()
        guard cell.isValid else { throw VivoChemistryError.invalid("Ewald periodic cell") }
        for source in sources { try source.validate() }
        let n = sources.count
        _ = try budget.elements([n,16], simultaneousArrays: 3)
        let charge = sources.reduce(0) { $0+$1.chargeE }
        guard charge.isFinite, cfg.chargeConvention != .requireNeutral || abs(charge) <= cfg.neutralToleranceE else {
            throw VivoChemistryError.invalid("nonneutral Ewald source set requires an explicit background convention")
        }
        let lattice = [cell.a,cell.b,cell.c].map { $0/VivoAtomicUnits.bohrInNM }
        let determinant = lattice[0].dot(lattice[1].cross(lattice[2])), volume = abs(determinant)
        guard volume.isFinite, volume > 0 else { throw VivoChemistryError.invalid("Ewald cell volume in atomic units") }
        let reciprocal = [lattice[1].cross(lattice[2])/determinant,lattice[2].cross(lattice[0])/determinant,lattice[0].cross(lattice[1])/determinant]
        let coordinates = sources.map(\.positionBohr), moments = sources.map(\.moments)
        var forces = [VivoVector3D](repeating: .zero,count: n)
        var derivatives = [[Double]](repeating: Array(repeating: 0,count: 10),count: n)
        var strain = VivoQMMatrix(3,3)
        var real = 0.0, reciprocalEnergy = 0.0, selfEnergy = 0.0, exceptionEnergy = 0.0, work = 0, modes = 0
        func xyz(_ v: VivoVector3D) -> [Double] { [v.x,v.y,v.z] }
        func accumulate(_ interaction: Pair, i: Int, j: Int, delta: VivoVector3D, scale: Double) {
            forces[i] = forces[i] - interaction.gradient*scale
            forces[j] = forces[j] + interaction.gradient*scale
            for a in 0..<10 {
                derivatives[i][a] += interaction.left[a]*scale; derivatives[j][a] += interaction.right[a]*scale
            }
            let g = xyz(interaction.gradient), d = xyz(delta)
            for a in 0..<3 { for b in 0..<3 { strain[a,b] += scale*g[a]*d[b] } }
        }
        // Center the finite image enumeration on a fractional reduction, then use
        // reciprocal-vector bounds to include EVERY image within the real cutoff.
        // This is not an orthogonal-cell assumption or a fixed 27-image search.
        for i in 0..<n { for j in 0...i {
            var reduced = coordinates[i]-coordinates[j]
            let fractions = reciprocal.map { $0.dot(reduced).rounded() }
            for a in 0..<3 { reduced = reduced-lattice[a]*fractions[a] }
            var ranges: [ClosedRange<Int>] = [], candidateCount = 1
            for axis in 0..<3 {
                let center = reciprocal[axis].dot(reduced), radius = reciprocal[axis].norm*cfg.realCutoffBohr
                let low = ceil(center-radius-1e-12), high = floor(center+radius+1e-12)
                guard low.isFinite, high.isFinite, low >= -4096, high <= 4096 else {
                    throw VivoChemistryError.resourceLimit("Ewald real-image integer range")
                }
                if low > high { candidateCount = 0; break }
                let count = Int(high-low)+1, product = candidateCount.multipliedReportingOverflow(by: count)
                guard !product.overflow, product.partialValue <= min(cfg.maximumImagePairs,budget.maximumOperatorApplications)-work else {
                    throw VivoChemistryError.resourceLimit("Ewald real-image work capacity")
                }
                candidateCount = product.partialValue; ranges.append(Int(low)...Int(high))
            }
            if candidateCount == 0 { continue }
            work += candidateCount
            for x in ranges[0] { for y in ranges[1] { for z in ranges[2] {
                if i == j && x == 0 && y == 0 && z == 0 { continue }
                let delta = reduced-lattice[0]*Double(x)-lattice[1]*Double(y)-lattice[2]*Double(z)
                let r = delta.norm
                if r > cfg.realCutoffBohr { continue }
                guard r > 1e-10 else { throw VivoChemistryError.invalid("coincident periodic multipole centers") }
                let interaction = pair(moments[i],moments[j],tensor: try .screened(delta,alpha: cfg.alphaPerBohr))
                let scale = i == j ? 0.5 : 1.0
                real += scale*interaction.energy; accumulate(interaction,i: i,j: j,delta: delta,scale: scale)
            } } }
        } }
        let widths = cfg.reciprocalHalfWidths
        let modeProduct = (2*widths[0]+1)*(2*widths[1]+1)*(2*widths[2]+1)-1
        let reciprocalWork = modeProduct.multipliedReportingOverflow(by: max(n,1))
        guard !reciprocalWork.overflow, reciprocalWork.partialValue <= budget.maximumOperatorApplications-work else {
            throw VivoChemistryError.resourceLimit("Ewald reciprocal mode/source capacity")
        }
        let alpha2 = cfg.alphaPerBohr*cfg.alphaPerBohr
        if let reciprocalOperator {
            let result = try reciprocalOperator.evaluate(sources,cell,cfg,budget)
            guard result.forcesHartreePerBohr.count == n,result.momentDerivatives.count == n,
                  result.momentDerivatives.allSatisfy({ $0.count == 10 }),
                  result.affineStrainDerivativeHartree.rows == 3,result.affineStrainDerivativeHartree.columns == 3 else {
                throw VivoChemistryError.invalid("reciprocal operator result shape")
            }
            reciprocalEnergy = result.energyHartree;modes = result.modeCount
            for i in 0..<n {
                forces[i] = forces[i]+result.forcesHartreePerBohr[i]
                for c in 0..<10 { derivatives[i][c] += result.momentDerivatives[i][c] }
            }
            strain = try strain.adding(result.affineStrainDerivativeHartree)
        } else {
        for x in -widths[0]...widths[0] { for y in -widths[1]...widths[1] { for z in -widths[2]...widths[2] {
            if x == 0 && y == 0 && z == 0 { continue }; modes += 1
            let k = (reciprocal[0]*Double(x)+reciprocal[1]*Double(y)+reciprocal[2]*Double(z))*(2*Double.pi)
            let kv = xyz(k), k2 = k.squaredNorm
            let factor = 2*Double.pi/volume*exp(-k2/(4*alpha2))/k2
            let coefficient: [SIMD2<Double>] = [.init(1,0),.init(0,k.x),.init(0,k.y),.init(0,k.z),
                .init(-0.5*k.x*k.x,0),.init(-k.x*k.y,0),.init(-k.x*k.z,0),
                .init(-0.5*k.y*k.y,0),.init(-k.y*k.z,0),.init(-0.5*k.z*k.z,0)]
            var phase = [SIMD2<Double>](), amplitudes = [SIMD2<Double>](), total = SIMD2<Double>.zero
            for i in 0..<n {
                let theta = k.dot(coordinates[i]), e = SIMD2<Double>(cos(theta),sin(theta))
                var local = SIMD2<Double>.zero
                for component in 0..<10 { local += coefficient[component]*moments[i][component] }
                let amplitude = complexProduct(e,local)
                phase.append(e); amplitudes.append(amplitude); total += amplitude
            }
            let norm2 = total.x*total.x+total.y*total.y, energy = factor*norm2
            reciprocalEnergy += energy
            for a in 0..<3 { for b in 0..<3 {
                strain[a,b] += energy*((a == b ? -1 : 0)+(2/k2+1/(2*alpha2))*kv[a]*kv[b])
            } }
            for i in 0..<n {
                let projection = total.x*amplitudes[i].y-total.y*amplitudes[i].x
                forces[i] = forces[i]+k*(2*factor*projection)
                for component in 0..<10 {
                    let d = complexProduct(phase[i],coefficient[component])
                    derivatives[i][component] += 2*factor*(total.x*d.x+total.y*d.y)
                }
                let q = sources[i].secondMomentsEBohr2, mu = xyz(sources[i].dipoleEBohr)
                let qk = [q[0]*k.x+q[1]*k.y+q[2]*k.z,q[1]*k.x+q[3]*k.y+q[4]*k.z,q[2]*k.x+q[4]*k.y+q[5]*k.z]
                for a in 0..<3 { for b in 0..<3 {
                    let change = complexProduct(phase[i],.init(kv[a]*qk[b],-kv[a]*mu[b]))
                    strain[a,b] += 2*factor*(total.x*change.x+total.y*change.y)
                } }
            }
        } } }
        }
        let smooth = VivoEwaldRadialTensor.selfKernel(alpha: cfg.alphaPerBohr)
        for i in 0..<n {
            let correction = pair(moments[i],moments[i],tensor: smooth)
            selfEnergy -= 0.5*correction.energy
            for component in 0..<10 { derivatives[i][component] -= 0.5*(correction.left[component]+correction.right[component]) }
        }
        let background = -Double.pi*charge*charge/(2*alpha2*volume)
        for i in 0..<n { derivatives[i][0] -= Double.pi*charge/(alpha2*volume) }
        for a in 0..<3 { strain[a,a] -= background }
        var pairs = Set<String>()
        for exception in primaryPairScales {
            let i = exception.first, j = exception.second
            guard i >= 0, j >= 0, i < n, j < n, i != j, exception.scale.isFinite, exception.scale >= 0,
                  pairs.insert("\(min(i,j)):\(max(i,j))").inserted else { throw VivoChemistryError.invalid("multipole primary exception") }
            // A molecular primary-pair exception is evaluated in its nearest image.
            let delta = try cell.minimumImage((coordinates[i]-coordinates[j])*VivoAtomicUnits.bohrInNM)/VivoAtomicUnits.bohrInNM
            let correction = pair(moments[i],moments[j],tensor: try .screened(delta,alpha: 0))
            let scale = exception.scale-1
            exceptionEnergy += scale*correction.energy; accumulate(correction,i: i,j: j,delta: delta,scale: scale)
        }
        let energy = real+reciprocalEnergy+selfEnergy+background+exceptionEnergy
        guard [energy,real,reciprocalEnergy,selfEnergy,background,exceptionEnergy].allSatisfy(\.isFinite),
              forces.allSatisfy(\.isFinite), derivatives.joined().allSatisfy(\.isFinite),strain.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite multipolar Ewald result")
        }
        return .init(energyHartree: energy, realEnergyHartree: real, reciprocalEnergyHartree: reciprocalEnergy,
            selfEnergyHartree: selfEnergy, backgroundEnergyHartree: background, exceptionEnergyHartree: exceptionEnergy,
            forcesHartreePerBohr: forces, momentDerivatives: derivatives, affineStrainDerivativeHartree: strain,
            imagePairEvaluations: work, reciprocalModeCount: modes)
    }
    static func complexProduct(_ a: SIMD2<Double>, _ b: SIMD2<Double>) -> SIMD2<Double> {
        .init(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x)
    }
}

/// Taylor coefficients in three Cartesian displacements, total degree <=5.
/// For g(t)=erfc(alpha*sqrt(t))/sqrt(t),
/// g^(n)(t)=(-1)^n Gamma(n+1/2,alpha^2*t)/(sqrt(pi)*t^(n+1/2)).
/// The upper-gamma recurrence avoids subtraction of nearly equal erf/r and 1/r.
struct VivoEwaldRadialTensor {
    private let derivatives: [Double]
    private static let exponents: [(Int,Int,Int)] = (0...5).flatMap { degree in
        (0...degree).flatMap { x in (0...(degree-x)).map { y in (x,y,degree-x-y) } }
    }
    private static func index(_ x: Int,_ y: Int,_ z: Int) -> Int { 36*x+6*y+z }
    subscript(_ x: Int,_ y: Int,_ z: Int) -> Double { derivatives[Self.index(x,y,z)] }
    private init(delta: VivoVector3D, radialCoefficients: [Double]) {
        var value = [Double](repeating: 0,count: 216), power = value
        power[0] = 1
        let terms: [(Int,Int,Int,Double)] = [(1,0,0,2*delta.x),(0,1,0,2*delta.y),(0,0,1,2*delta.z),(2,0,0,1),(0,2,0,1),(0,0,2,1)]
        for degree in 0...5 {
            for (x,y,z) in Self.exponents { let i = Self.index(x,y,z); value[i] += radialCoefficients[degree]*power[i] }
            if degree == 5 { break }
            var next = [Double](repeating: 0,count: 216)
            for (x,y,z) in Self.exponents {
                let p = power[Self.index(x,y,z)]
                if p == 0 { continue }
                for (a,b,c,w) in terms where x+y+z+a+b+c <= 5 {
                    next[Self.index(x+a,y+b,z+c)] += p*w
                }
            }
            power = next
        }
        let factorial = [1.0,1,2,6,24,120]
        for (x,y,z) in Self.exponents { value[Self.index(x,y,z)] *= factorial[x]*factorial[y]*factorial[z] }
        derivatives = value
    }
    static func screened(_ delta: VivoVector3D, alpha: Double) throws -> Self {
        let t = delta.squaredNorm, r = sqrt(t), x = alpha*alpha*t
        guard delta.isFinite,t.isFinite,t > 1e-20,alpha.isFinite,alpha >= 0,x.isFinite else {
            throw VivoChemistryError.invalid("multipolar radial kernel geometry")
        }
        var gamma = sqrt(Double.pi)*erfc(alpha*r), factorial = 1.0, tPower = r
        var coefficients: [Double] = []
        for n in 0...5 {
            coefficients.append((n%2 == 0 ? 1 : -1)*gamma/(sqrt(Double.pi)*tPower*factorial))
            gamma = (Double(n)+0.5)*gamma+pow(x,Double(n)+0.5)*exp(-x)
            tPower *= t; factorial *= Double(n+1)
        }
        return .init(delta: delta,radialCoefficients: coefficients)
    }
    static func selfKernel(alpha: Double) -> Self {
        var coefficients = [Double](repeating: 0,count: 6), power = 1.0, factorial = 1.0
        for n in 0...5 {
            coefficients[n] = 2*alpha/sqrt(Double.pi)*power/(factorial*Double(2*n+1))
            power *= -alpha*alpha; factorial *= Double(n+1)
        }
        return .init(delta: .zero,radialCoefficients: coefficients)
    }
}
