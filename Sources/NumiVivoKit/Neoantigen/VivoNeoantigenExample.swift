import Foundation

public enum VivoNeoantigenExample {
    /// Invented records for software checks only. Not HCC1395, patient data,
    /// experimentally validated epitopes, or output from an executed pVACseq run.
    public static func inputs() throws -> (manifest: VivoNeoantigenCase, tsv: Data) {
        let header = ["Chromosome", "Start", "Stop", "Reference", "Variant", "Transcript", "Gene Name",
                      "HLA Allele", "Peptide Length", "MT Epitope Seq", "WT Epitope Seq", "Median MT IC50 Score",
                      "Tumor DNA Depth", "Tumor DNA VAF", "Tumor RNA Depth", "Tumor RNA VAF", "Normal Depth", "Normal VAF", "Gene Expression", "Evaluation"]
        let rows = [
            ["1", "100", "101", "A", "T", "DEMO_TX_1", "DEMO_GENE_1", "HLA-A*02:01", "9", "ACDEFGHIK", "ACDEYGHIK", "125", "80", "0.4", "50", "0.3", "60", "0", "12", "Accept"],
            ["2", "200", "201", "C", "G", "DEMO_TX_2", "DEMO_GENE_2", "HLA-A*02:01", "9", "LMNPQRSTV", "LMNAQRSTV", "30", "70", "0.2", "NA", "NA", "65", "0", "NA", "Accept"],
            ["3", "300", "301", "G", "A", "DEMO_TX_3", "DEMO_GENE_3", "HLA-A*02:01", "9", "WYACDEFGH", "WYACDEFGH", "NA", "80", "0.3", "30", "0", "60", "0.1", "0", "Review"]
        ]
        let tsv = Data(([header] + rows).map { $0.joined(separator: "\t") }.joined(separator: "\n").appending("\n").utf8)
        func hash(_ text: String) throws -> String { try VivoNeoantigenWorkbench.digest(Data(text.utf8)) }
        let reference = try hash("synthetic reference descriptor; no reference genome supplied")
        let samples = try VivoNeoantigenCase.Sample.Role.allCases.map { role in
            VivoNeoantigenCase.Sample(role: role, sampleID: "demo-" + role.rawValue, subjectID: "synthetic-subject",
                assembly: "GRCh38", referenceSHA256: reference, sourceSHA256: try hash("synthetic sample descriptor " + role.rawValue))
        }
        let run = VivoNeoantigenCase.ExternalRun(tool: "pvacseq", version: "0.0.0+synthetic",
            format: VivoNeoantigenWorkbench.adapterFormat, executionImageSHA256: try hash("no container executed"),
            arguments: ["synthetic-fixture", "not-an-executed-command"], modelVersions: ["demonstration": "0.0.0-synthetic"],
            sampleInputs: Dictionary(uniqueKeysWithValues: samples.map { ($0.sampleID, $0.sourceSHA256) }),
            annotatedVCFSHA256: try hash("no VCF processed"), reportSHA256: try VivoNeoantigenWorkbench.digest(tsv))
        let manifest = VivoNeoantigenCase(schema: "numivivo.org/neoantigen-case/v1", caseID: "synthetic-demo",
            subjectID: "synthetic-subject", dataClass: .synthetic,
            sourceCitation: "Invented software fixture. No biological measurements, real reference case or pVACseq execution.",
            assembly: "GRCh38", referenceSHA256: reference, annotationRelease: "synthetic-annotation-v1",
            annotationSHA256: try hash("synthetic annotation descriptor"), hlaNomenclatureRelease: "synthetic-hla-v1",
            hlaSampleID: "demo-matchedNormalDNA", hlaAlleles: ["HLA-A*02:01"], samples: samples, externalRun: run)
        return (manifest, tsv)
    }
}
