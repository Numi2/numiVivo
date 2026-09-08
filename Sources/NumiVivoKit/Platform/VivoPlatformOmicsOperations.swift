import Foundation

public enum VivoPlatformOmicsOperations {
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [VivoPlatformOperations.pure(identifier: "vivo.platform.singlecell", id: id,
            inputs: ["input": VivoSingleCellArtifacts.inputKind],
            outputs: [.init(name: "report", kind: "vivo.singlecell-report-v1")],
            summary: "Exact count MEX import, sample-aware QC, optional log normalization and replicate-aware raw pseudobulk; no differential expression or inferred cell type.",
            configure: VivoPlatformOperations.empty, calculate: { _, inputs, _ in
                let input = try VivoPlatformOperations.input(VivoSingleCellInputBundle.self, "input", inputs)
                // The manifest can lower but cannot silently raise native count,
                // feature, cell and source-byte admission limits.
                return ["report": try VivoCanonicalJSON.encode(VivoSingleCellCampaign.evaluate(input))]
            })]
    }
}
