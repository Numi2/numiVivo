import Foundation
@testable import NumiVivoKit
@main struct LifecycleMain {
    static func main() throws {
        guard CommandLine.arguments.count == 4 else { throw VivoOmicsError.invalid("gaussian-lifecycle ORIGINAL_COHORT OUTPUT OWNER_RESULT") }
        let root = URL(fileURLWithPath: CommandLine.arguments[1]), out = URL(fileURLWithPath: CommandLine.arguments[2]), qualified = URL(fileURLWithPath: CommandLine.arguments[3])
        let input = root.appendingPathComponent("pca")
        let prior = try VivoPCAIntegration.read(VivoH5ADPCAReceipt.self,input,"receipt.json",maximum:65_536)
        let originalPlan = try VivoPCAIntegration.read(VivoPCAIntegrationPlan.self,root.appendingPathComponent("mnn"),"plan.json",maximum:65_536)
        guard var options = originalPlan.mnn else { throw VivoOmicsError.invalid("MNN plan required") }
        options.kernel = .tiledGaussian
        let start = Date()
        let receipt = try VivoPCAIntegration.publish(input: input,plan: .init(mnn:options),implementation:prior.implementation,to:out)
        let publishedSeconds = Date().timeIntervalSince(start)
        let verifyStart = Date()
        guard try VivoPCAIntegration.verify(out,implementation:prior.implementation) == receipt else { throw VivoOmicsError.invalid("lifecycle replay") }
        for name in ["scores.bin","anchors.bin","report.json"] {
            guard try VivoOmicsFileSnapshot.fingerprint(out.appendingPathComponent(name),maximumBytes:1_000_000_000) == VivoOmicsFileSnapshot.fingerprint(qualified.appendingPathComponent(name),maximumBytes:1_000_000_000) else { throw VivoOmicsError.invalid("owner and public route differ") }
        }
        let record: [String:Any] = ["status":"passed","publicationSeconds":publishedSeconds,"verificationSeconds":Date().timeIntervalSince(verifyStart),"allOwnerScoresAnchorsAndReportExact":true,"completeParentCountAndPCAReconstruction":true,"originalImplementationTagPreserved":true,"scope":"Public API lifecycle on original complete Hagai, with actual tested executable bound separately; other cohorts qualified through native numerical owner."]
        try JSONSerialization.data(withJSONObject:record,options:[.prettyPrinted,.sortedKeys]).write(to:out.deletingLastPathComponent().appendingPathComponent("lifecycle-checks.json"))
        print(String(decoding:try JSONSerialization.data(withJSONObject:record,options:[.sortedKeys]),as:UTF8.self))
    }
}
