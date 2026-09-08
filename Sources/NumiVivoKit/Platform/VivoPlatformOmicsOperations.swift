import Foundation

public enum VivoPlatformOmicsOperations {
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [VivoPlatformOperations.pure(identifier: "vivo.platform.singlecell", id: id,
            inputs: ["input": VivoSingleCellArtifacts.inputKind],
            outputs: [.init(name: "report", kind: "vivo.singlecell-report-v1")],
            summary: "Exact count MEX/gzip import, sample-aware QC metrics, optional log normalization and replicate-aware raw pseudobulk.",
            configure: VivoPlatformOperations.empty, calculate: { _, inputs, _ in
                let input = try VivoPlatformOperations.input(VivoSingleCellInputBundle.self, "input", inputs)
                return ["report": try VivoCanonicalJSON.encode(VivoSingleCellCampaign.evaluate(input))]
            }),
         VivoPlatformOperations.pure(identifier: "vivo.platform.singlecell-analyze", id: id,
            inputs: ["input": VivoSingleCellArtifacts.inputKind, "plan": VivoSingleCellAnalysisArtifacts.planKind],
            outputs: [.init(name: "report", kind: "vivo.singlecell-cohort-report-v1")],
            summary: "Explicit cell selection, source-row lineage, feature statistics and optional donor-aware pseudobulk linear contrasts.",
            configure: VivoPlatformOperations.empty, calculate: { _, inputs, _ in
                let input = try VivoPlatformOperations.input(VivoSingleCellInputBundle.self, "input", inputs)
                guard let bytes = inputs["plan"] else { throw VivoOmicsError.invalid("missing analysis plan") }
                let plan = try VivoSingleCellAnalysisArtifacts.decodePlan(bytes)
                let source = try VivoSingleCellCampaign.evaluate(input)
                return ["report": try VivoCanonicalJSON.encode(VivoSingleCellCohortAnalysis.run(source.dataset, plan: plan))]
            })]
    }
}
