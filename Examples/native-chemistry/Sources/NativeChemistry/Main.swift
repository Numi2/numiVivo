import Foundation
#if !NUMIVIVO_PORTABLE_CHEMISTRY
import NumiVivoKit
#endif

private struct Request: Codable {
    let system: VivoElectronicSystem
    let basis: VivoGaussianBasis
    let scf: VivoSCFConfiguration
    let budget: VivoChemistryBudget
}
private struct Report: Codable {
    let schema: String
    let status: String
    let request: Request
    let hartreeFock: VivoHartreeFockResult
    let mp2: VivoMP2Result?
    let ci: VivoCIResult
    let embeddedHamiltonian: VivoEmbeddedHamiltonian
    let orbitalEntropiesNats: [Double]
    let mutualInformation: VivoQMMatrix
}
@main struct NativeChemistry {
    static func main() throws {
        guard CommandLine.arguments.count==3 else {
            throw VivoChemistryError.invalid("usage: native-chemistry <request.json|--h2> <report.json>")
        }
        let request:Request
        if CommandLine.arguments[1]=="--h2" {
            request = .init(system:.init(nuclei:[.init(atomicNumber:1,positionBohr:.zero),.init(atomicNumber:1,positionBohr:.init(0,0,1.4))],
                                           alphaElectrons:1,betaElectrons:1),basis:.hydrogenSTO3G(nucleusIndices:[0,1]),scf:.init(),budget:.init())
        } else {
            let input=URL(fileURLWithPath:CommandLine.arguments[1]),attributes=try FileManager.default.attributesOfItem(atPath:input.path)
            guard let size=attributes[.size] as? NSNumber,size.intValue<=64*1024*1024 else { throw VivoChemistryError.resourceLimit("example request file exceeds 64 MiB") }
            request=try JSONDecoder().decode(Request.self,from:Data(contentsOf:input))
        }
        guard request.scf.reference == .restricted else {
            throw VivoChemistryError.unsupported("this example's shared-spatial-orbital post-HF pipeline requires RHF; call native UHF separately")
        }
        let ao=try VivoGaussianIntegralEngine.compute(system:request.system,basis:request.basis,budget:request.budget)
        let hf=try VivoHartreeFock.solve(system:request.system,integrals:ao,configuration:request.scf,budget:request.budget)
        let h=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:hf.alphaCoefficients,
            alphaElectrons:request.system.alphaElectrons,betaElectrons:request.system.betaElectrons,
            orbitalIdentifiers:(0..<ao.count).map{"canonical-mo-\($0)"},energyReference:"E_electronic+E_nucleus-nucleus+E_nucleus-MM; no frozen-MM-self",budget:request.budget)
        let mp2=try VivoRestrictedMP2.solve(h,budget:request.budget),ci=try VivoConfigurationInteraction.solve(h,budget:request.budget)
        let information=try VivoCIDensityMatrices.orbitalInformation(ci.state,budget:request.budget)
        let report=Report(schema:"numivivo.org/native-chemistry-example/v1",status:"native-computation; not a reaction-barrier or chemical-accuracy qualification",
            request:request,hartreeFock:hf,mp2:mp2,ci:ci,embeddedHamiltonian:h,
            orbitalEntropiesNats:information.entropy,mutualInformation:information.mutualInformation)
        let encoder=JSONEncoder();encoder.outputFormatting=[.prettyPrinted,.sortedKeys,.withoutEscapingSlashes]
        try encoder.encode(report).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.atomic)
        print("RHF: \(hf.energyHartree) Hartree")
        print("MP2: \(mp2.totalEnergyHartree) Hartree")
        print("FCI: \(ci.energyHartree) Hartree")
        print("SCF residual: \(hf.finalCommutatorNorm)")
    }
}
