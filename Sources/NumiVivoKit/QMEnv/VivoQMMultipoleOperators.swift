import Foundation

/// Density-linear distributed RAW moments about each nucleus. An AO product is
/// assigned equally to its two basis centers; if both are on one center its full
/// product belongs there. This partition is explicit and sums to the exact AO
/// density, including off-diagonal products and their charge/overlap response.
struct VivoQMMultipoleOperators {
    let system: VivoElectronicSystem
    let orbitals: [VivoCartesianOrbital]
    let matrices: [[VivoQMMatrix]]
    var count: Int { orbitals.count }
    static func ownership(_ center: Int,_ a: VivoCartesianOrbital,_ b: VivoCartesianOrbital) -> Double {
        (a.nucleusIndex == center ? 0.5 : 0)+(b.nucleusIndex == center ? 0.5 : 0)
    }
    static func primitiveMoment(_ a: Double,_ la: [Int],_ ra: SIMD3<Double>,
                                _ b: Double,_ lb: [Int],_ rb: SIMD3<Double>,
                                origin: SIMD3<Double>, powers: [Int]) -> Double {
        // (r-C)^m = ((r-A)+(A-C))^m. m<=2 in each component here.
        var value = 0.0
        for x in 0...powers[0] { for y in 0...powers[1] { for z in 0...powers[2] {
            let raised = [x,y,z]
            var angular = la,weight = 1.0
            for axis in 0..<3 {
                angular[axis] += raised[axis]
                let binomial = powers[axis] == 2 && raised[axis] == 1 ? 2.0 : 1.0
                weight *= binomial*pow(ra[axis]-origin[axis],Double(powers[axis]-raised[axis]))
            }
            value += weight*VivoGaussianIntegralEngine.primitiveOverlap(a,angular,ra,b,lb,rb)
        } } }
        return value
    }
    init(integrals: VivoAOIntegrals,budget: VivoChemistryBudget) throws {
        let n = integrals.count,system = integrals.sourceSystem,orbitals = integrals.orbitals
        guard system.pointCharges.isEmpty else { throw VivoChemistryError.invalid("periodic moment operators require isolated QM integrals") }
        _ = try budget.elements([system.nuclei.count,10,n,n])
        var operators = [[VivoQMMatrix]](repeating: Array(repeating: VivoQMMatrix(n,n),count: 10),count: system.nuclei.count)
        var work = 0
        for p in 0..<n { for q in 0...p {
            let a = orbitals[p],b = orbitals[q],ra = system.nuclei[a.nucleusIndex].positionBohr,rb = system.nuclei[b.nucleusIndex].positionBohr
            for center in Set([a.nucleusIndex,b.nucleusIndex]).sorted() {
                let weight = Self.ownership(center,a,b),origin = system.nuclei[center].positionBohr
                for component in 0..<10 {
                    let cost = 4*a.weights.count*b.weights.count
                    guard cost <= budget.maximumOperatorApplications-work else { throw VivoChemistryError.resourceLimit("distributed multipole integral budget") }
                    work += cost
                    let value = weight*VivoGaussianDerivative.contracted(a,b,ra: ra,rb: rb) { ea,la,aa,eb,lb,bb in
                        Self.primitiveMoment(ea,la,aa,eb,lb,bb,origin: origin,powers: VivoPeriodicElectrostatics.powers[component])
                    }
                    operators[center][component][p,q] = value; operators[center][component][q,p] = value
                }
            }
        } }
        self.system = system; self.orbitals = orbitals; matrices = operators
    }
    func sources(density: VivoQMMatrix) throws -> [VivoCartesianMultipole] {
        guard density.rows == count,density.columns == count,density.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("multipole AO density")
        }
        return system.nuclei.indices.map { center in
            var moments = (0..<10).map { component in
                -zip(density.values,matrices[center][component].values).reduce(0.0) { $0+$1.0*$1.1 }
            }
            moments[0] += Double(system.nuclei[center].atomicNumber)
            let p = system.nuclei[center].positionBohr
            return .init(positionBohr: .init(p.x,p.y,p.z),chargeE: moments[0],dipoleEBohr: .init(moments[1],moments[2],moments[3]),
                secondMomentsEBohr2: Array(moments.dropFirst(4)))
        }
    }
    func fock(momentDerivatives lambda: [[Double]]) throws -> VivoQMMatrix {
        guard lambda.count == system.nuclei.count,lambda.allSatisfy({ $0.count == 10 && $0.allSatisfy(\.isFinite) }) else {
            throw VivoChemistryError.invalid("multipolar embedding potential shape")
        }
        var potential = VivoQMMatrix(count,count)
        for center in matrices.indices { for component in 0..<10 {
            for i in potential.values.indices { potential.values[i] -= lambda[center][component]*matrices[center][component].values[i] }
        } }
        return potential
    }
    /// dE/dR from AO-product and moving moment-origin response, holding the AO
    /// density matrix fixed. Source-center forces are accounted for separately.
    func geometryGradient(density: VivoQMMatrix,momentDerivatives lambda: [[Double]],
                          budget: VivoChemistryBudget) throws -> [VivoVector3D] {
        _ = try fock(momentDerivatives: lambda)
        var gradient = [VivoVector3D](repeating: .zero,count: system.nuclei.count),work = 0
        for p in 0..<count { for q in 0..<count {
            let a = orbitals[p],b = orbitals[q],ia = a.nucleusIndex,ib = b.nucleusIndex
            let ra = system.nuclei[ia].positionBohr,rb = system.nuclei[ib].positionBohr
            for center in Set([ia,ib]).sorted() {
                let weight = -Self.ownership(center,a,b)*density[p,q],origin = system.nuclei[center].positionBohr
                for component in 0..<10 {
                    let scalar = weight*lambda[center][component]
                    if scalar == 0 { continue }
                    let powers = VivoPeriodicElectrostatics.powers[component]
                    let integral: VivoGaussianDerivative.PrimitivePair = { ea,la,aa,eb,lb,bb in
                        Self.primitiveMoment(ea,la,aa,eb,lb,bb,origin: origin,powers: powers)
                    }
                    for axis in 0..<3 {
                        let cost = 20*a.weights.count*b.weights.count
                        guard cost <= budget.maximumOperatorApplications-work else { throw VivoChemistryError.resourceLimit("multipole geometry derivative budget") }
                        work += cost
                        let ga = scalar*VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: true,operator: integral)
                        let gb = scalar*VivoGaussianDerivative.pair(a,b,ra: ra,rb: rb,axis: axis,onFirst: false,operator: integral)
                        VivoGaussianDerivative.add(&gradient,center: ia,axis: axis,value: ga)
                        VivoGaussianDerivative.add(&gradient,center: ib,axis: axis,value: gb)
                        if powers[axis] > 0 {
                            var lower = powers; lower[axis] -= 1
                            let gc = -scalar*Double(powers[axis])*VivoGaussianDerivative.contracted(a,b,ra: ra,rb: rb) { ea,la,aa,eb,lb,bb in
                                Self.primitiveMoment(ea,la,aa,eb,lb,bb,origin: origin,powers: lower)
                            }
                            VivoGaussianDerivative.add(&gradient,center: center,axis: axis,value: gc)
                        }
                    }
                }
            }
        } }
        return gradient
    }
}
