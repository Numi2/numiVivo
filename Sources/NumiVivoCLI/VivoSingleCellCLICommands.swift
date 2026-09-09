import Foundation
import NumiVivoKit

struct VivoSingleCellCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["singlecell-pca-integrate", "singlecell-pca-integrate-verify", "singlecell-graph-embed", "singlecell-graph-embed-verify", "singlecell-graph-cluster", "singlecell-graph-cluster-verify", "singlecell-pca-neighbors", "singlecell-pca-neighbors-verify", "singlecell-h5ad-pca-query", "singlecell-h5ad-pca-query-verify", "singlecell-h5ad-pca", "singlecell-h5ad-pca-verify", "singlecell-h5ad-store", "singlecell-count-store-verify", "singlecell-count-store-normalize", "singlecell-count-store-normalize-verify", "multiassay-h5mu-import", "multiassay-h5mu-write", "multiassay-10x-import", "multiassay-verify", "singlecell-composition-prepare", "singlecell-composition-fit", "singlecell-composition-predict", "singlecell-composition-verify", "singlecell-composition-prediction-verify", "singlecell-perturbation-fit", "singlecell-perturbation-predict", "singlecell-perturbation-verify", "singlecell-perturbation-prediction-verify", "singlecell-reference-fit", "singlecell-reference-map", "singlecell-reference-verify", "singlecell-reference-map-verify", "singlecell-h5ad-project", "singlecell-h5ad-project-verify", "singlecell-h5ad-pseudobulk", "singlecell-h5ad-pseudobulk-verify", "singlecell-h5ad-annotate", "singlecell-h5ad-import", "singlecell-h5ad-write", "singlecell-run", "singlecell-verify", "singlecell-export", "singlecell-mex", "singlecell-help", "singlecell-example",
         "singlecell-analyze", "singlecell-analysis-verify", "singlecell-analysis-export", "singlecell-analysis-mex", "singlecell-analysis-tables"].contains(name ?? "")
    }
    private func canonicalURL(_ url: URL) throws -> URL {
        let manager = FileManager.default
        var ancestor = url.absoluteURL, suffix: [String] = []
        guard ancestor.isFileURL, ancestor.path.utf8.count <= 8192 else { throw VivoOmicsError.invalid("bounded local output path required") }
        while !manager.fileExists(atPath: ancestor.path) {
            if let attributes = try? manager.attributesOfItem(atPath: ancestor.path),
               attributes[.type] as? FileAttributeType == .typeSymbolicLink { throw VivoOmicsError.invalid("dangling symbolic link in output path") }
            let parent = ancestor.deletingLastPathComponent()
            guard parent.path != ancestor.path, suffix.count < 4096 else { throw VivoOmicsError.invalid("output has no resolvable ancestor") }
            suffix.append(ancestor.lastPathComponent); ancestor = parent
        }
        var resolved = ancestor.resolvingSymlinksInPath().standardizedFileURL
        for part in suffix.reversed() { resolved.appendPathComponent(part) }
        return resolved.standardizedFileURL
    }
    private func load<T: Decodable>(_ type: T.Type, _ url: URL) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoSingleCellCampaignIO.readDocument(url, maximumBytes: 128 * 1_024))
    }
    private func printJSON<T: Encodable>(_ value: T) throws {
        FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(value)); FileHandle.standardOutput.write(Data("\n".utf8))
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoOmicsError.invalid("missing command") }
            if command == "singlecell-help" {
                guard arguments.count == 1 else { throw VivoOmicsError.invalid("singlecell-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            if command == "singlecell-example" {
                guard arguments.count == 3, arguments[1] == "--output", arguments[2] != "-" else {
                    throw VivoOmicsError.invalid("singlecell-example --output <new-directory>")
                }
                let destination = try canonicalURL(URL(fileURLWithPath: arguments[2]))
                var files = try VivoSingleCellMEXExchange.files(VivoSingleCellExamples.pairedCounts())
                files["analysis.json"] = try VivoCanonicalJSON.encode(VivoSingleCellExamples.pairedPlan())
                try VivoOmicsDirectoryExport.write(files, to: destination)
                try printJSON(["status": "written-synthetic-example", "directory": destination.path]); return 0
            }
            if command == "singlecell-h5ad-store" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("singlecell-h5ad-store <source.h5ad> --plan <mapping.json> --output <new-store>") }
                let plan = try VivoCanonicalJSON.decode(VivoH5ADImportPlan.self, from: VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]), maximumBytes: 2_097_152))
                try printJSON(VivoH5ADCountStore.publish(source: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-count-store-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-count-store-verify <store>") }
                try printJSON(VivoH5ADCountStore.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-count-store-normalize" {
                guard arguments.count == 6, arguments[2] == "--target", arguments[4] == "--output", let target = Double(arguments[3]) else { throw VivoOmicsError.invalid("singlecell-count-store-normalize <store> --target <positive-total> --output <new-directory>") }
                try printJSON(VivoH5ADCountStore.normalize(URL(fileURLWithPath: arguments[1]), target: target, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-count-store-normalize-verify" {
                guard arguments.count == 4, arguments[2] == "--store" else { throw VivoOmicsError.invalid("singlecell-count-store-normalize-verify <normalized> --store <raw-store>") }
                try printJSON(VivoH5ADCountStore.verifyNormalization(URL(fileURLWithPath: arguments[1]), store: URL(fileURLWithPath: arguments[3]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "multiassay-h5mu-import" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("multiassay-h5mu-import <source.h5mu> --plan <mapping.json> --output <new-bundle>") }
                let bytes = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]), maximumBytes: 2_097_152)
                let plan = try VivoCanonicalJSON.decode(VivoH5MUMultiAssayPlan.self, from: bytes)
                try printJSON(VivoMultiAssayIO.importH5MU(source: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "multiassay-h5mu-write" {
                guard arguments.count == 4, arguments[2] == "--output" else { throw VivoOmicsError.invalid("multiassay-h5mu-write <dataset.json> --output <new.h5mu>") }
                let bytes = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[1]), maximumBytes: 536_870_912)
                let data = try VivoCanonicalJSON.decode(VivoMultiAssayDataset.self, from: bytes)
                try VivoMultiAssayH5MU.write(data, to: URL(fileURLWithPath: arguments[3])); try printJSON(["status": "exported-h5mu"]); return 0
            }
            if command == "multiassay-10x-import" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("multiassay-10x-import <matrix.h5> --plan <plan.json> --output <new-bundle>") }
                let plan = try load(VivoTenXMultiAssayPlan.self, URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoMultiAssayIO.importTenX(source: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "multiassay-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("multiassay-verify <bundle>") }
                try printJSON(VivoMultiAssayIO.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-composition-prepare" {
                guard arguments.count==6,arguments[2]=="--plan",arguments[4]=="--output" else { throw VivoOmicsError.invalid("singlecell-composition-prepare <pseudobulk-bundle> --plan <selection.json> --output <training.json>") }
                let plan=try load(VivoCompositionPreparation.self,URL(fileURLWithPath: arguments[3]))
                try VivoComposition.prepare(from: URL(fileURLWithPath: arguments[1]),plan: plan,implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[5])))
                try printJSON(["status":"prepared-composition-training"]);return 0
            }
            if command == "singlecell-composition-fit" {
                guard arguments.count==4,arguments[2]=="--output" else { throw VivoOmicsError.invalid("singlecell-composition-fit <training.json> --output <model-bundle>") }
                try printJSON(VivoComposition.fit(source: URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[3]))));return 0
            }
            if command == "singlecell-composition-predict" {
                guard arguments.count==6,arguments[2]=="--plan",arguments[4]=="--output" else { throw VivoOmicsError.invalid("singlecell-composition-predict <model-bundle> --plan <queries.json> --output <prediction-bundle>") }
                let plan=try load(VivoCompositionQueryPlan.self,URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoComposition.predict(reference: URL(fileURLWithPath: arguments[1]),plan: plan,implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[5]))));return 0
            }
            if command == "singlecell-composition-verify" || command == "singlecell-composition-prediction-verify" {
                guard arguments.count==2 else { throw VivoOmicsError.invalid("composition verifier requires one bundle") }
                if command == "singlecell-composition-verify" { _=try VivoComposition.verifyModel(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint()) }
                else { _=try VivoComposition.verifyPrediction(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint()) }
                try printJSON(["status":"verified-composition-bundle"]);return 0
            }
            if command == "singlecell-perturbation-fit" {
                guard arguments.count==6,arguments[2]=="--plan",arguments[4]=="--output" else { throw VivoOmicsError.invalid("singlecell-perturbation-fit <training.h5ad> --plan <fit.json> --output <new-bundle>") }
                let plan=try load(VivoPerturbationPlan.self,URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoPerturbation.fit(source: URL(fileURLWithPath: arguments[1]),plan: plan,implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[5]))));return 0
            }
            if command == "singlecell-perturbation-predict" {
                guard arguments.count==8,arguments[2]=="--plan",arguments[4]=="--reference",arguments[6]=="--output" else { throw VivoOmicsError.invalid("singlecell-perturbation-predict <query.h5ad> --plan <query.json> --reference <reference-bundle> --output <new-bundle>") }
                let plan=try load(VivoPerturbationQueryPlan.self,URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoPerturbation.map(source: URL(fileURLWithPath: arguments[1]),plan: plan,reference: URL(fileURLWithPath: arguments[5]),implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[7]))));return 0
            }
            if command == "singlecell-perturbation-verify" || command == "singlecell-perturbation-prediction-verify" {
                guard arguments.count==2 else { throw VivoOmicsError.invalid("reference verifier requires one bundle directory") }
                if command == "singlecell-perturbation-verify" { _=try VivoPerturbation.verifyModel(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint()) }
                else { _=try VivoPerturbation.verifyPrediction(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint()) }
                try printJSON(["status":"verified-reference-bundle"]);return 0
            }
            if command == "singlecell-reference-fit" {
                guard arguments.count==6,arguments[2]=="--plan",arguments[4]=="--output" else { throw VivoOmicsError.invalid("singlecell-reference-fit <training.h5ad> --plan <fit.json> --output <new-bundle>") }
                let plan=try load(VivoSingleCellReferencePlan.self,URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoSingleCellReference.fit(source: URL(fileURLWithPath: arguments[1]),plan: plan,implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[5]))));return 0
            }
            if command == "singlecell-reference-map" {
                guard arguments.count==8,arguments[2]=="--plan",arguments[4]=="--reference",arguments[6]=="--output" else { throw VivoOmicsError.invalid("singlecell-reference-map <query.h5ad> --plan <query.json> --reference <reference-bundle> --output <new-bundle>") }
                let plan=try load(VivoSingleCellReferenceQueryPlan.self,URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoSingleCellReference.map(source: URL(fileURLWithPath: arguments[1]),plan: plan,reference: URL(fileURLWithPath: arguments[5]),implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[7]))));return 0
            }
            if command == "singlecell-reference-verify" || command == "singlecell-reference-map-verify" {
                guard arguments.count==2 else { throw VivoOmicsError.invalid("reference verifier requires one bundle directory") }
                if command == "singlecell-reference-verify" { _=try VivoSingleCellReference.verifyReference(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint()) }
                else { _=try VivoSingleCellReference.verifyMapping(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint()) }
                try printJSON(["status":"verified-reference-bundle"]);return 0
            }
            if command == "singlecell-h5ad-project" {
                guard arguments.count==6,arguments[2]=="--plan",arguments[4]=="--output" else {
                    throw VivoOmicsError.invalid("singlecell-h5ad-project <source.h5ad> --plan <projection.json> --output <new-directory>")
                }
                let bytes=try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]),maximumBytes: 64*1_024*1_024)
                let plan=try VivoCanonicalJSON.decode(VivoH5ADProjectionPlan.self,from: bytes)
                try printJSON(VivoH5ADProjection.publish(source: URL(fileURLWithPath: arguments[1]),plan: plan,
                    implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[5]))));return 0
            }
            if command == "singlecell-h5ad-project-verify" {
                guard arguments.count==2 else { throw VivoOmicsError.invalid("singlecell-h5ad-project-verify <bundle-directory>") }
                let report=try VivoH5ADProjection.verify(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint())
                try printJSON(report);return 0
            }
            if command == "singlecell-pca-integrate" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("singlecell-pca-integrate <PCA-bundle> --plan <plan.json> --output <new-bundle>") }
                let plan = try load(VivoPCAIntegrationPlan.self, URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoPCAIntegration.publish(input: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-pca-integrate-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-pca-integrate-verify <bundle>") }
                try printJSON(VivoPCAIntegration.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-graph-embed" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("singlecell-graph-embed <binary-graph> --plan <plan.json> --output <new-bundle>") }
                let plan = try load(VivoPCAGraphEmbeddingPlan.self, URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoPCAGraphEmbedding.publish(input: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-graph-embed-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-graph-embed-verify <bundle>") }
                try printJSON(VivoPCAGraphEmbedding.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-graph-cluster" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("singlecell-graph-cluster <binary-graph> --plan <plan.json> --output <new-bundle>") }
                let plan = try load(VivoPCAGraphClusteringPlan.self, URL(fileURLWithPath: arguments[3]))
                try printJSON(VivoPCAGraphClustering.publish(input: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-graph-cluster-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-graph-cluster-verify <bundle>") }
                try printJSON(VivoPCAGraphClustering.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-pca-neighbors" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("singlecell-pca-neighbors <pca-bundle> --plan <plan.json> --output <new-bundle>") }
                let bytes = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]), maximumBytes: 65_536)
                let plan = try VivoCanonicalJSON.decode(VivoPCANeighborPlan.self, from: bytes)
                try printJSON(VivoPCANeighborBundle.publish(input: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-pca-neighbors-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-pca-neighbors-verify <bundle>") }
                try printJSON(VivoPCANeighborBundle.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-h5ad-pca-query" {
                guard arguments.count == 8, arguments[2] == "--plan", arguments[4] == "--reference", arguments[6] == "--output" else { throw VivoOmicsError.invalid("singlecell-h5ad-pca-query <query.h5ad> --plan <plan.json> --reference <pca-bundle> --output <new-bundle>") }
                let bytes = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]), maximumBytes: 2_097_152)
                let plan = try VivoCanonicalJSON.decode(VivoH5ADPCAQueryPlan.self, from: bytes)
                try printJSON(VivoH5ADPCAQuery.publish(source: URL(fileURLWithPath: arguments[1]), plan: plan, reference: URL(fileURLWithPath: arguments[5]), implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[7])))); return 0
            }
            if command == "singlecell-h5ad-pca-query-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-h5ad-pca-query-verify <bundle>") }
                try printJSON(VivoH5ADPCAQuery.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-h5ad-pca" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else { throw VivoOmicsError.invalid("singlecell-h5ad-pca <source.h5ad> --plan <plan.json> --output <new-bundle>") }
                let bytes = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]), maximumBytes: 2_097_152)
                let plan = try VivoCanonicalJSON.decode(VivoH5ADPCAPlan.self, from: bytes)
                try printJSON(VivoH5ADPCA.publish(source: URL(fileURLWithPath: arguments[1]), plan: plan, implementation: VivoWorkflowCLIImplementation.fingerprint(), to: canonicalURL(URL(fileURLWithPath: arguments[5])))); return 0
            }
            if command == "singlecell-h5ad-pca-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-h5ad-pca-verify <bundle>") }
                try printJSON(VivoH5ADPCA.verify(URL(fileURLWithPath: arguments[1]), implementation: VivoWorkflowCLIImplementation.fingerprint())); return 0
            }
            if command == "singlecell-h5ad-pseudobulk-verify" {
                guard arguments.count == 2 else { throw VivoOmicsError.invalid("singlecell-h5ad-pseudobulk-verify <bundle-directory>") }
                let report=try VivoH5ADPseudobulk.verify(URL(fileURLWithPath: arguments[1]),implementation: VivoWorkflowCLIImplementation.fingerprint())
                try printJSON(["status":"verified-streamed-pseudobulk","sourceCells":String(report.metadata.cells.count),"nonzeros":String(report.canonicalNonzeros)])
                return 0
            }
            if command == "singlecell-h5ad-pseudobulk" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else {
                    throw VivoOmicsError.invalid("singlecell-h5ad-pseudobulk <source.h5ad> --plan <stream-plan.json> --output <new-bundle-directory>")
                }
                let plan=try load(VivoH5ADPseudobulkPlan.self,URL(fileURLWithPath: arguments[3]))
                let receipt=try VivoH5ADPseudobulk.publish(source: URL(fileURLWithPath: arguments[1]),plan: plan,
                    implementation: VivoWorkflowCLIImplementation.fingerprint(),to: canonicalURL(URL(fileURLWithPath: arguments[5])))
                try printJSON(receipt); return 0
            }
            if command == "singlecell-h5ad-annotate" {
                guard arguments.count == 6, arguments[2] == "--plan", arguments[4] == "--output" else {
                    throw VivoOmicsError.invalid("singlecell-h5ad-annotate <source.h5ad> --plan <annotations.json> --output <new.h5ad>")
                }
                let bytes = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: arguments[3]), maximumBytes: 64 * 1_024 * 1_024)
                let plan = try VivoCanonicalJSON.decode(VivoH5ADAnnotationPlan.self, from: bytes)
                let destination = try canonicalURL(URL(fileURLWithPath: arguments[5]))
                try printJSON(VivoSingleCellH5AD.annotate(URL(fileURLWithPath: arguments[1]), plan: plan,
                    implementation: VivoWorkflowCLIImplementation.fingerprint(), to: destination))
                return 0
            }
            if command == "singlecell-h5ad-import" || command == "singlecell-h5ad-write" {
                let importing = command == "singlecell-h5ad-import"
                guard arguments.count == (importing ? 6 : 4), arguments[2] == (importing ? "--plan" : "--output"),
                      !importing || arguments[4] == "--output" else {
                    throw VivoOmicsError.invalid("singlecell-h5ad-import <source.h5ad> --plan <mapping.json> --output <new-directory>; singlecell-h5ad-write <dataset.json> --output <new.h5ad>")
                }
                let source = URL(fileURLWithPath: arguments[1])
                let destination = try canonicalURL(URL(fileURLWithPath: arguments[importing ? 5 : 3]))
                if importing {
                    let plan = try load(VivoH5ADImportPlan.self, URL(fileURLWithPath: arguments[3]))
                    let document = try VivoSingleCellH5AD.read(source, plan: plan)
                    var files = try VivoSingleCellMEXExchange.files(document.dataset)
                    files["original.h5ad"] = document.source
                    files["dataset.json"] = try VivoCanonicalJSON.encode(document.dataset)
                    files["h5ad-mapping.json"] = try VivoCanonicalJSON.encode(plan)
                    struct ImportReceipt: Encodable {
                        let schema = "numivivo.org/h5ad-import/v1"
                        let source: VivoFingerprint; let dataset: VivoFingerprint; let mapping: VivoFingerprint
                        let implementation: VivoFingerprint
                        let hdf5Version: String
                    }
                    files["h5ad-import.json"] = try VivoCanonicalJSON.encode(ImportReceipt(
                        source: VivoCanonicalJSON.fingerprint(document.source),
                        dataset: VivoCanonicalJSON.fingerprint(files["dataset.json"]!),
                        mapping: VivoCanonicalJSON.fingerprint(files["h5ad-mapping.json"]!),
                        implementation: VivoWorkflowCLIImplementation.fingerprint(), hdf5Version: document.hdf5Version))
                    try VivoOmicsDirectoryExport.write(files, to: destination)
                    try printJSON(["status": "imported-count-projection-with-original-anndata", "directory": destination.path])
                } else {
                    let bytes = try VivoSingleCellCampaignIO.readDocument(source, maximumBytes: 128 * 1_024 * 1_024)
                    let dataset = try VivoCanonicalJSON.decode(VivoSingleCellDataset.self, from: bytes)
                    try VivoSingleCellH5AD.write(dataset, to: destination)
                    try printJSON(["status": "written-native-count-anndata", "file": destination.path])
                }
                return 0
            }
            guard arguments.count >= 4, !arguments[1].hasPrefix("--"), arguments[1] != "-" else {
                throw VivoOmicsError.invalid("one manifest/receipt file and --store <directory> are required")
            }
            let inputURL = URL(fileURLWithPath: arguments[1]).standardizedFileURL
            var options: [String: String] = [:], index = 2
            while index < arguments.count {
                let key = arguments[index]
                guard ["--store", "--output", "--plan"].contains(key), options[key] == nil,
                      index + 1 < arguments.count, !arguments[index + 1].isEmpty,
                      !arguments[index + 1].hasPrefix("--") else { throw VivoOmicsError.invalid("unknown, duplicate or incomplete option") }
                options[key] = arguments[index + 1]; index += 2
            }
            guard (command == "singlecell-analyze") == (options["--plan"] != nil) else {
                throw VivoOmicsError.invalid("--plan is required only for singlecell-analyze")
            }
            guard let storePath = options["--store"], storePath != "-" else { throw VivoOmicsError.invalid("--store <directory> is required") }
            let storeURL = URL(fileURLWithPath: storePath).standardizedFileURL
            let directoryCommand = ["singlecell-mex", "singlecell-analysis-mex", "singlecell-analysis-tables"].contains(command)
            if directoryCommand, options["--output"] == nil || options["--output"] == "-" {
                throw VivoOmicsError.invalid("directory exports require --output <new-directory>")
            }
            var outputURL: URL?
            if let output = options["--output"], output != "-" {
                let url = try canonicalURL(URL(fileURLWithPath: output)), root = try canonicalURL(storeURL).path
                var inputs = [inputURL]
                if let plan = options["--plan"] { inputs.append(URL(fileURLWithPath: plan)) }
                guard !inputs.contains(where: { $0.resolvingSymlinksInPath().standardizedFileURL == url }),
                      url.path != root, !url.path.hasPrefix(root == "/" ? "/" : root + "/"),
                      !FileManager.default.fileExists(atPath: url.path) else {
                    throw VivoOmicsError.invalid("output must be new, outside the artifact store, and not an input")
                }
                outputURL = url
            }
            let implementation = try VivoWorkflowCLIImplementation.fingerprint()
            let store = try VivoArtifactStore(rootURL: storeURL, createIfNeeded: command == "singlecell-run")
            let output: Data
            switch command {
            case "singlecell-run":
                let input = try VivoSingleCellCampaignIO.snapshot(manifestURL: inputURL)
                output = try await VivoCanonicalJSON.encode(VivoSingleCellArtifacts.publish(input: input, implementation: implementation, store: store))
            case "singlecell-analyze":
                let receipt = try load(VivoSingleCellRunReceipt.self, inputURL)
                let plan = try VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: options["--plan"]!), maximumBytes: 2 * 1_024 * 1_024)
                output = try await VivoCanonicalJSON.encode(VivoSingleCellAnalysisArtifacts.publish(counts: receipt, planBytes: plan, implementation: implementation, store: store))
            case "singlecell-verify", "singlecell-export", "singlecell-mex":
                let receipt = try load(VivoSingleCellRunReceipt.self, inputURL)
                let report = try await VivoSingleCellArtifacts.verify(receipt: receipt, implementation: implementation, store: store)
                if command == "singlecell-mex" {
                    var files = try VivoSingleCellMEXExchange.files(report.dataset)
                    files["source-receipt.json"] = try VivoCanonicalJSON.encode(receipt)
                    try VivoOmicsDirectoryExport.write(files, to: outputURL!)
                    try printJSON(["status": "written-verified-counts", "directory": outputURL!.path]); return 0
                }
                if command == "singlecell-export" { output = try VivoCanonicalJSON.encode(report) }
                else {
                    struct Verification: Encodable {
                        let schemaVersion: Int; let status: String; let evidence: VivoOmicsEvidence
                        let cells: Int; let features: Int; let pseudobulkGroups: Int; let numericalProfile: String
                    }
                    output = try VivoCanonicalJSON.encode(Verification(schemaVersion: 1, status: "verified-native-reconstruction-not-biological-validation",
                        evidence: report.dataset.evidence, cells: report.dataset.cells.count, features: report.dataset.features.count,
                        pseudobulkGroups: report.pseudobulk.groups.count, numericalProfile: report.numericalProfile))
                }
            case "singlecell-analysis-verify", "singlecell-analysis-export", "singlecell-analysis-mex", "singlecell-analysis-tables":
                let receipt = try load(VivoSingleCellAnalysisReceipt.self, inputURL)
                let report = try await VivoSingleCellAnalysisArtifacts.verify(receipt, implementation: implementation, store: store)
                if directoryCommand {
                    var files: [String: Data]
                    if command == "singlecell-analysis-tables" { files = try VivoSingleCellAnalysisTables.files(report, receipt: receipt) }
                    else {
                        files = try VivoSingleCellMEXExchange.files(report.processed.dataset)
                        files["analysis-receipt.json"] = try VivoCanonicalJSON.encode(receipt)
                        files["processed-to-original-cell-indices.json"] = try VivoCanonicalJSON.encode(report.processed.sourceCellIndices)
                    }
                    try VivoOmicsDirectoryExport.write(files, to: outputURL!)
                    try printJSON(["status": "written-verified-analysis", "directory": outputURL!.path]); return 0
                }
                if command == "singlecell-analysis-export" { output = try VivoCanonicalJSON.encode(report) }
                else {
                    struct Verification: Encodable { let status: String; let acceptedCells: Int; let rejectedCells: Int; let contrasts: Int; let testedFeatures: [Int] }
                    output = try VivoCanonicalJSON.encode(Verification(status: "verified-native-analysis-not-biological-calibration",
                        acceptedCells: report.processed.dataset.cells.count,
                        rejectedCells: report.processed.decisions.count - report.processed.dataset.cells.count,
                        contrasts: report.contrasts.count, testedFeatures: report.contrasts.map(\.testedFeatures)))
                }
            default: throw VivoOmicsError.invalid("unsupported single-cell command")
            }
            try Task.checkCancellation()
            if let url = outputURL { try VivoKineticsDocumentIO.write(output, to: url, overwrite: false) }
            else { FileHandle.standardOutput.write(output); FileHandle.standardOutput.write(Data("\n".utf8)) }
            return 0
        } catch is CancellationError {
            FileHandle.standardError.write(Data("numivivo singlecell: cancelled; no success receipt published\n".utf8)); return 130
        } catch {
            FileHandle.standardError.write(Data("numivivo singlecell: \(error)\n".utf8)); return 65
        }
    }
    static let help = """
    NumiVivo native single-cell workflows
      singlecell-h5ad-store <source.h5ad> --plan <mapping.json> --output <new-store>
    singlecell-count-store-verify <store>
    singlecell-count-store-normalize <store> --target <positive-total> --output <new-directory>
    singlecell-count-store-normalize-verify <normalized> --store <raw-store>
    multiassay-h5mu-import <source.h5mu> --plan <mapping.json> --output <new-bundle>
      multiassay-h5mu-write <dataset.json> --output <new.h5mu>
      multiassay-10x-import <matrix.h5> --plan <plan.json> --output <new-bundle>
      multiassay-verify <bundle>
      singlecell-composition-prepare <pseudobulk-bundle> --plan <selection.json> --output <training.json>
      singlecell-composition-fit <training.json> --output <model-bundle>
      singlecell-composition-predict <model-bundle> --plan <queries.json> --output <prediction-bundle>
      singlecell-composition-verify <model-bundle>
      singlecell-composition-prediction-verify <prediction-bundle>
      singlecell-perturbation-fit <training.h5ad> --plan <fit.json> --output <new-model>
      singlecell-perturbation-predict <control.h5ad> --plan <query.json> --reference <model> --output <new-bundle>
      singlecell-perturbation-verify <model>
      singlecell-perturbation-prediction-verify <prediction-bundle>
      singlecell-reference-fit <training.h5ad> --plan <fit.json> --output <new-reference>
      singlecell-reference-map <query.h5ad> --plan <query.json> --reference <reference> --output <new-bundle>
      singlecell-reference-verify <reference>
      singlecell-reference-map-verify <mapping-bundle>
      singlecell-pca-integrate <PCA-bundle> --plan <plan.json> --output <new-bundle>
      singlecell-pca-integrate-verify <bundle>
      singlecell-graph-embed <binary-graph> --plan <plan.json> --output <new-bundle>
      singlecell-graph-embed-verify <bundle>
      singlecell-graph-cluster <binary-graph> --plan <plan.json> --output <new-bundle>
      singlecell-graph-cluster-verify <bundle>
      singlecell-pca-neighbors <pca-bundle> --plan <plan.json> --output <new-bundle>
      singlecell-pca-neighbors-verify <bundle>
      singlecell-h5ad-pca-query <query.h5ad> --plan <plan.json> --reference <pca-bundle> --output <new-bundle>
      singlecell-h5ad-pca-query-verify <bundle>
      singlecell-h5ad-pca <source.h5ad> --plan <plan.json> --output <new-bundle>
      singlecell-h5ad-pca-verify <bundle>
      singlecell-h5ad-pseudobulk <source.h5ad> --plan <stream-plan.json> --output <new-bundle-directory>
      singlecell-h5ad-project <source.h5ad> --plan <projection.json> --output <new-directory>
      singlecell-h5ad-project-verify <bundle-directory>
      singlecell-h5ad-pseudobulk-verify <bundle-directory>
      singlecell-h5ad-annotate <source.h5ad> --plan <annotations.json> --output <new.h5ad>
      singlecell-h5ad-import <source.h5ad> --plan <mapping.json> --output <new-directory>
      singlecell-h5ad-write <dataset.json> --output <new.h5ad>
      singlecell-example --output <new-example-directory>
      singlecell-run <manifest.json> --store <directory> [--output <new-receipt.json|->]
      singlecell-verify <receipt.json> --store <directory>
      singlecell-export <receipt.json> --store <directory> [--output <new-report.json|->]
      singlecell-mex <receipt.json> --store <directory> --output <new-MEX-directory>
      singlecell-analyze <receipt.json> --plan <analysis.json> --store <directory> [--output <new-analysis-receipt.json|->]
      singlecell-analysis-verify <analysis-receipt.json> --store <directory>
      singlecell-analysis-export <analysis-receipt.json> --store <directory> [--output <new-report.json|->]
      singlecell-analysis-mex <analysis-receipt.json> --store <directory> --output <new-MEX-directory>
      singlecell-analysis-tables <analysis-receipt.json> --store <directory> --output <new-table-directory>
    Source paths are relative to the manifest directory; symlinks and ../ are rejected.
    Input is integer/general MEX Gene Expression data, plain or gzip with bounded native decoding.
    Counts stay UInt64; QC records every cell decision and keeps feature identities.
    Expression uses explicit biological replicates/paired donors, not cells as independent samples.
    The moderated log-expression model is untrended: no voom weights, mixed model, or NB fit is claimed.
    Replay requires the recorded executable/OS. Regenerate count receipts after rebuilding the executable.
    H5AD uses native HDF5 (install hdf5 or set NUMIVIVO_HDF5_LIBRARY), with explicit X/raw/X/layer and design mapping; CSR/CSC and row-wise dense reads.
    Import retains original.h5ad unchanged alongside the count projection and MEX manifest; no ancillary AnnData fields are discarded.
    H5AD write creates a new count object from dataset.json; it does not copy source embeddings onto changed cells.
    No doublet correction, inferred cell annotation or biological calibration is performed.
    """ + "\n"
}
