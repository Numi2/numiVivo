import Foundation

public struct VivoElectronicCenterGradient: Codable, Sendable, Equatable {
    public let nuclearGradientsHartreePerBohr: [VivoVector3D]
    public let pointChargeGradientsHartreePerBohr: [VivoVector3D]
}

/// Analytic Cartesian Gaussian derivatives use the existing integral engine:
/// d phi_l / d A = 2*a*phi_(l+1) - l*phi_(l-1). Contracted normalization and
/// exponents are held fixed. No numerical displacement or second SCF is used.
enum VivoGaussianDerivative {
    typealias PrimitivePair = (Double,[Int],SIMD3<Double>,Double,[Int],SIMD3<Double>) -> Double
    static func pair(_ a: VivoCartesianOrbital, _ b: VivoCartesianOrbital,
                     ra: SIMD3<Double>, rb: SIMD3<Double>, axis: Int, onFirst: Bool,
                     operator integral: PrimitivePair) -> Double {
        var value = 0.0
        for i in a.weights.indices { for j in b.weights.indices {
            let ea = a.primitiveExponents[i], eb = b.primitiveExponents[j], weight = a.weights[i]*b.weights[j]
            var la = a.angular, lb = b.angular
            if onFirst {
                la[axis] += 1
                value += weight*2*ea*integral(ea,la,ra,eb,lb,rb)
                if a.angular[axis] > 0 {
                    la[axis] -= 2
                    value -= weight*Double(a.angular[axis])*integral(ea,la,ra,eb,lb,rb)
                }
            } else {
                lb[axis] += 1
                value += weight*2*eb*integral(ea,la,ra,eb,lb,rb)
                if b.angular[axis] > 0 {
                    lb[axis] -= 2
                    value -= weight*Double(b.angular[axis])*integral(ea,la,ra,eb,lb,rb)
                }
            }
        } }
        return value
    }
    static func contracted(_ a: VivoCartesianOrbital,_ b: VivoCartesianOrbital,
                           ra: SIMD3<Double>,rb: SIMD3<Double>,operator integral: PrimitivePair) -> Double {
        var value = 0.0
        for i in a.weights.indices { for j in b.weights.indices {
            value += a.weights[i]*b.weights[j]*integral(a.primitiveExponents[i],a.angular,ra,b.primitiveExponents[j],b.angular,rb)
        } }
        return value
    }
    /// Potential of a point source with charge and Cartesian dipole. Source
    /// derivatives follow translation invariance: d/dC = -(d/dA+d/dB).
    /// The raised primitives retain the ORIGINAL contracted normalization.
    static func multipolarPotential(_ ea: Double,_ la: [Int],_ ra: SIMD3<Double>,
                                    _ eb: Double,_ lb: [Int],_ rb: SIMD3<Double>,
                                    center: SIMD3<Double>,moments: [Double]) -> Double {
        var result = moments[0]*VivoGaussianIntegralEngine.primitivePotential(ea,la,ra,eb,lb,rb,center)
        for axis in 0..<3 where moments[axis+1] != 0 {
            var aa=la,bb=lb;aa[axis]+=1;bb[axis]+=1
            var derivative = -2*ea*VivoGaussianIntegralEngine.primitivePotential(ea,aa,ra,eb,lb,rb,center)
                - 2*eb*VivoGaussianIntegralEngine.primitivePotential(ea,la,ra,eb,bb,rb,center)
            if la[axis]>0 { aa[axis]-=2;derivative+=Double(la[axis])*VivoGaussianIntegralEngine.primitivePotential(ea,aa,ra,eb,lb,rb,center) }
            if lb[axis]>0 { bb[axis]-=2;derivative+=Double(lb[axis])*VivoGaussianIntegralEngine.primitivePotential(ea,la,ra,eb,bb,rb,center) }
            result += moments[axis+1]*derivative
        }
        return result
    }
    static func add(_ values: inout [VivoVector3D], center: Int, axis: Int, value: Double) {
        switch axis { case 0: values[center].x += value; case 1: values[center].y += value; default: values[center].z += value }
    }
}

public enum VivoHartreeFockNuclearDerivatives {
    /// Includes one- and two-electron derivatives, nuclear/charge reactions, and
    /// the energy-weighted overlap (Pulay) term. For an embedded SCF, reference
    /// orbital energies MUST come from the complete embedded Fock operator.
    public static func evaluate(integrals ao: VivoAOIntegrals, reference: VivoHartreeFockResult,
                                budget: VivoChemistryBudget = .init()) throws -> VivoElectronicCenterGradient {
        let n = ao.count
        guard reference.alphaElectrons == ao.sourceSystem.alphaElectrons,
              reference.betaElectrons == ao.sourceSystem.betaElectrons,
              reference.alphaCoefficients.rows == n, reference.alphaCoefficients.columns == n,
              reference.betaCoefficients.rows == n, reference.betaCoefficients.columns == n,
              reference.alphaOrbitalEnergies.count == n, reference.betaOrbitalEnergies.count == n,
              (reference.alphaOrbitalEnergies+reference.betaOrbitalEnergies).allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("HF derivative orbital/electron identity")
        }
        var weighted = VivoQMMatrix(n,n)
        for p in 0..<n { for q in 0..<n {
            for i in 0..<reference.alphaElectrons {
                weighted[p,q] += reference.alphaCoefficients[p,i]*reference.alphaCoefficients[q,i]*reference.alphaOrbitalEnergies[i]
            }
            for i in 0..<reference.betaElectrons {
                weighted[p,q] += reference.betaCoefficients[p,i]*reference.betaCoefficients[q,i]*reference.betaOrbitalEnergies[i]
            }
        } }
        return try evaluate(integrals: ao,alphaDensity: reference.alphaDensity,betaDensity: reference.betaDensity,
                            energyWeightedDensity: weighted,budget: budget)
    }

    public static func evaluate(integrals ao: VivoAOIntegrals, alphaDensity da: VivoQMMatrix,
                                betaDensity db: VivoQMMatrix, energyWeightedDensity weighted: VivoQMMatrix,
                                budget: VivoChemistryBudget = .init()) throws -> VivoElectronicCenterGradient {
        try ao.validate(budget: budget)
        let system = ao.sourceSystem, n = ao.count, orbitals = ao.orbitals
        for matrix in [da,db,weighted] {
            guard matrix.rows == n,matrix.columns == n,matrix.values.allSatisfy(\.isFinite),
                  try matrix.adding(matrix.transposed,scale: -1).frobeniusNorm < 1e-8 else {
                throw VivoChemistryError.invalid("HF derivative density shape, symmetry or values")
            }
        }
        let density = try da.adding(db)
        var nuclei = [VivoVector3D](repeating: .zero,count: system.nuclei.count)
        var charges = [VivoVector3D](repeating: .zero,count: system.pointCharges.count)
        var work = 0
        func reserve(_ count: Int) throws {
            guard count >= 0,count <= budget.maximumOperatorApplications-work else {
                throw VivoChemistryError.resourceLimit("analytic Gaussian derivative primitive-work budget")
            }
            work += count
        }
        for p in 0..<n { for q in 0..<n {
            let a = orbitals[p],b = orbitals[q],ia = a.nucleusIndex,ib = b.nucleusIndex
            let ra = system.nuclei[ia].positionBohr,rb = system.nuclei[ib].positionBohr
            let occupation = density[p,q],pulay = weighted[p,q]
            if occupation == 0 && pulay == 0 { continue }
            for axis in 0..<3 {
                for first in [true,false] {
                    try reserve(4*a.weights.count*b.weights.count)
                    let ds = VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: first,
                        operator: VivoGaussianIntegralEngine.primitiveOverlap)
                    let dt = VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: first,
                        operator: VivoGaussianIntegralEngine.primitiveKinetic)
                    VivoGaussianDerivative.add(&nuclei,center: first ? ia : ib,axis: axis,value: occupation*dt-pulay*ds)
                }
                for center in 0..<(system.nuclei.count+system.pointCharges.count) {
                    let nuclear = center < system.nuclei.count
                    let charge = nuclear ? Double(system.nuclei[center].atomicNumber) : system.pointCharges[center-system.nuclei.count].chargeE
                    if charge == 0 || occupation == 0 { continue }
                    let rc = nuclear ? system.nuclei[center].positionBohr : system.pointCharges[center-system.nuclei.count].positionBohr
                    try reserve(4*a.weights.count*b.weights.count)
                    let integral: VivoGaussianDerivative.PrimitivePair = { ea,la,aa,eb,lb,bb in
                        VivoGaussianIntegralEngine.primitivePotential(ea,la,aa,eb,lb,bb,rc)
                    }
                    let ga = -charge*occupation*VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: true,operator: integral)
                    let gb = -charge*occupation*VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: false,operator: integral)
                    VivoGaussianDerivative.add(&nuclei,center: ia,axis: axis,value: ga)
                    VivoGaussianDerivative.add(&nuclei,center: ib,axis: axis,value: gb)
                    if nuclear { VivoGaussianDerivative.add(&nuclei,center: center,axis: axis,value: -ga-gb) }
                    else { VivoGaussianDerivative.add(&charges,center: center-system.nuclei.count,axis: axis,value: -ga-gb) }
                }
            }
        } }
        // Full chemist-order contraction is explicit. This bounded reference can
        // later share shell screening without changing its derivative convention.
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            let coefficient = 0.5*(density[p,q]*density[r,s]-da[p,r]*da[q,s]-db[p,r]*db[q,s])
            if coefficient == 0 { continue }
            let row = [orbitals[p],orbitals[q],orbitals[r],orbitals[s]]
            let positions = row.map { system.nuclei[$0.nucleusIndex].positionBohr }
            // A common displacement of all four centers leaves this integral unchanged.
            if Set(row.map(\.nucleusIndex)).count == 1 { continue }
            for i in row[0].weights.indices { for j in row[1].weights.indices {
                for k in row[2].weights.indices { for l in row[3].weights.indices {
                    try reserve(24)
                    let indices = [i,j,k,l],exponents = (0..<4).map { row[$0].primitiveExponents[indices[$0]] }
                    let weight = coefficient*(0..<4).reduce(1.0) { $0*row[$1].weights[indices[$1]] }
                    if weight == 0 { continue }
                    func value(_ angular: [[Int]]) -> Double {
                        VivoGaussianIntegralEngine.primitiveERI(exponents[0],angular[0],positions[0],exponents[1],angular[1],positions[1],
                            exponents[2],angular[2],positions[2],exponents[3],angular[3],positions[3])
                    }
                    for center in 0..<4 { for axis in 0..<3 {
                        var angular = row.map(\.angular),g = 0.0
                        angular[center][axis] += 1
                        g = 2*exponents[center]*value(angular)
                        if row[center].angular[axis] > 0 {
                            angular[center][axis] -= 2
                            g -= Double(row[center].angular[axis])*value(angular)
                        }
                        VivoGaussianDerivative.add(&nuclei,center: row[center].nucleusIndex,axis: axis,value: weight*g)
                    } }
                } }
            } }
        } } } }
        for i in system.nuclei.indices {
            let position = system.nuclei[i].positionBohr,charge = Double(system.nuclei[i].atomicNumber)
            for j in 0..<i {
                let d = position-system.nuclei[j].positionBohr,r = vivoQMNorm(d)
                let g = VivoVector3D(d.x,d.y,d.z)*(-charge*Double(system.nuclei[j].atomicNumber)/(r*r*r))
                nuclei[i] = nuclei[i]+g; nuclei[j] = nuclei[j]-g
            }
            for j in system.pointCharges.indices where system.pointCharges[j].chargeE != 0 {
                let d = position-system.pointCharges[j].positionBohr,r = vivoQMNorm(d)
                let g = VivoVector3D(d.x,d.y,d.z)*(-charge*system.pointCharges[j].chargeE/(r*r*r))
                nuclei[i] = nuclei[i]+g; charges[j] = charges[j]-g
            }
        }
        guard nuclei.allSatisfy(\.isFinite),charges.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite analytic HF center derivative")
        }
        return .init(nuclearGradientsHartreePerBohr: nuclei,pointChargeGradientsHartreePerBohr: charges)
    }
}
