import Foundation
#if NUMIVIVO_SCOPED_NUMERICS
@testable import NumiVivoNumerics
#endif

/// Executes native production numerical sources against independently generated
/// PySCF fixtures. A mismatch is a failing exit status, not a skipped comparison.
@main struct ExternalConformance {
    struct Fixture: Codable {
        let identifier: String
        let system: VivoElectronicSystem
        let basis: VivoGaussianBasis
        let auxiliaryBasis: VivoGaussianBasis
        let overlap: VivoQMMatrix
        let coreHamiltonian: VivoQMMatrix
        let fittedERI: [Double]
        let embedded: VivoEmbeddedHamiltonian
        let reference: Reference
    }
    struct Reference: Codable {
        let hf: Double; let mp2: Double; let ccsd: Double
        let fciRoots: [Double]; let riHF: Double; let riMP2: Double; let riRank: Int
        let smoothCPCM: Double?; let smoothPolarization: Double?
        let saConverged: Bool?; let saEnergy: Double?; let saStates: [Double]?
    }
    struct Check: Codable {
        let system: String; let observable: String; let actual: Double
        let expected: Double; let tolerance: Double
        var error: Double { abs(actual - expected) }
        var passed: Bool { actual.isFinite && expected.isFinite && error <= tolerance }
    }
    struct Report: Codable {
        let schema: String; let passed: Bool; let checks: [Check]
        let failures: [String]; let scope: String
    }
    static func main() throws {
        guard CommandLine.arguments.count == 3 else { throw VivoChemistryError.invalid("expected fixture directory and output JSON") }
        let root = URL(fileURLWithPath: CommandLine.arguments[1])
        var checks: [Check] = [], failures: [String] = []
        for name in ["h2", "lih", "water"] {
            do {
                let f = try JSONDecoder().decode(Fixture.self, from: Data(contentsOf: root.appendingPathComponent(name+".json")))
                guard f.identifier==name,f.reference.fciRoots.count==3 else {throw VivoChemistryError.invalid("fixture identity/root count")}
                if name=="lih" {
                    guard f.reference.saConverged==true,f.reference.saEnergy != nil,f.reference.saStates?.count==2 else{throw VivoChemistryError.invalid("required independent state-average references missing")}
                } else {
                    guard f.reference.smoothCPCM != nil,f.reference.smoothPolarization != nil else{throw VivoChemistryError.invalid("required independent solvent references missing")}
                }
                let budget = VivoChemistryBudget(maximumBytes:512*1024*1024,maximumBasisFunctions:256,
                    maximumDeterminants:100_000,maximumOperatorApplications:1_000_000_000)
                func check(_ observable:String,_ value:Double,_ expected:Double,_ tolerance:Double) {
                    let c = Check(system:name,observable:observable,actual:value,expected:expected,tolerance:tolerance)
                    checks.append(c)
                    print("\(c.passed ? "PASS" : "FAIL") \(name) \(observable) error=\(c.error) tolerance=\(tolerance)")
                }
                func difference(_ a:[Double],_ b:[Double]) throws -> Double {
                    guard a.count == b.count else { throw VivoChemistryError.invalid("conformance array dimensions") }
                    return zip(a,b).map { abs($0-$1) }.max() ?? 0
                }
                let ao = try VivoGaussianIntegralEngine.compute(system:f.system,basis:f.basis,budget:budget)
                check("overlap",try difference(ao.overlap.values,f.overlap.values),0,1e-9)
                check("coreHamiltonian",try difference(ao.coreHamiltonian.values,f.coreHamiltonian.values),0,1e-9)
                let hf = try VivoHartreeFock.solve(system:f.system,integrals:ao,budget:budget)
                check("RHF",hf.energyHartree,f.reference.hf,1e-8)
                let nativeH = try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:hf.alphaCoefficients,
                    alphaElectrons:f.system.alphaElectrons,betaElectrons:f.system.betaElectrons,
                    orbitalIdentifiers:f.embedded.orbitalIdentifiers,energyReference:f.embedded.energyReference,budget:budget)
                let mp2 = try VivoRestrictedMP2.solve(nativeH,budget:budget)
                check("MP2",mp2.totalEnergyHartree,f.reference.mp2,1e-8)
                let cc = try VivoTensorCCSD.solve(nativeH,budget:budget)
                guard cc.converged else { throw VivoChemistryError.convergence("native tensor CCSD") }
                check("tensorCCSD",cc.energyHartree,f.reference.ccsd,1e-8)
                let roots = try VivoDirectCI.solve(nativeH,configuration:.init(roots:3,maximumSubspace:60),budget:budget)
                for (i,r) in roots.roots.enumerated() { check("FCI-root-\(i)",r.energyHartree,f.reference.fciRoots[i],1e-8) }
                let ri = try VivoDensityFitting.compute(system:f.system,basis:f.basis,auxiliaryBasis:f.auxiliaryBasis,budget:budget)
                var eri: [Double] = []
                for p in 0..<ao.count { for q in 0..<ao.count { for r in 0..<ao.count { for s in 0..<ao.count {
                    eri.append(ri.factors.eri(p,q,r,s))
                } } } }
                check("RI-ERI",try difference(eri,f.fittedERI),0,1e-8)
                check("RI-rank",Double(ri.factors.rank),Double(f.reference.riRank),0)
                let rif = try VivoFactorizedHartreeFock.solve(ri,budget:budget)
                check("RI-HF",rif.scf.energyHartree,f.reference.riHF,1e-8)
                let rimp = try VivoFactorizedMP2.solve(rif,budget:budget)
                check("RI-MP2",rimp.totalEnergyHartree,f.reference.riMP2,1e-8)
                if let expected = f.reference.smoothCPCM {
                    let solvated = try VivoSmoothSolventSCF.solve(system:f.system,integrals:ao,budget:budget)
                    check("smooth-CPCM-RHF",solvated.scf.energyHartree,expected,1e-8)
                    guard let polarization = f.reference.smoothPolarization else { throw VivoChemistryError.invalid("missing independent solvent reference") }
                    check("smooth-CPCM-polarization",solvated.solvent.polarizationEnergyHartree,polarization,1e-8)
                    check("smooth-CPCM-residual",solvated.solvent.trueRelativeResidual,0,5e-10)
                }
                if let expected = f.reference.saEnergy {
                    guard f.reference.saConverged == true else { throw VivoChemistryError.convergence("oracle state-average") }
                    let cfg = VivoMultiStateCASSCFConfiguration(weights:[0.5,0.5],rootLabels:["ground","excited"],
                        optimization:.init(maximumMacroIterations:200,maximumEnergyEvaluations:4000,gradientTolerance:1e-6),followRoots:true)
                    let cas = try VivoMultiStateCASSCF.solve(nativeH,partition:.init(doublyOccupiedCore:[0],active:[1,2]),configuration:cfg,budget:budget)
                    guard cas.converged else { throw VivoChemistryError.convergence("native state-average \(cas.termination)") }
                    check("SA-CASSCF",cas.weightedEnergyHartree,expected,5e-7)
                    if let energies = f.reference.saStates {
                        let actual = cas.states.map(\.energyHartree).sorted()
                        for i in actual.indices { check("SA-CASSCF-state-\(i)",actual[i],energies.sorted()[i],5e-6) }
                    }
                }
            } catch { failures.append("\(name): \(error)"); print("ERROR \(name): \(error)") }
        }
        let report = Report(schema:"numivivo.org/native-external-conformance/v1",passed:failures.isEmpty && checks.count==45 && checks.allSatisfy(\.passed),
            checks:checks,failures:failures,scope:"Production numerical sources versus pinned independent PySCF; not a complete Apple package or paper-reaction qualification")
        let encoder=JSONEncoder();encoder.outputFormatting=[.prettyPrinted,.sortedKeys]
        try encoder.encode(report).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.atomic)
        if !report.passed { exit(1) }
    }
}
