import Foundation
import NumiVivoKit

@main struct H5ADCheck {
    static func main() {
        do { try run() }
        catch {
            FileHandle.standardError.write(Data("H5AD rejected: \(error)\n".utf8))
            exit(65)
        }
    }
    static func run() throws {
        let args = CommandLine.arguments
        let implementation=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
        if args.count == 5,args[1] == "count-store" {
            let plan=try JSONDecoder().decode(VivoH5ADImportPlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            let receipt=try VivoH5ADCountStore.publish(source: URL(fileURLWithPath: args[2]),plan: plan,implementation: implementation,to: URL(fileURLWithPath: args[4]))
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt));return
        }
        if args.count == 3,args[1] == "verify-count-store" {
            _=try VivoH5ADCountStore.verify(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count == 6, args[1] == "normalize-count-store", let target = Double(args[3]),
           let backend = VivoCountStoreNormalizationBackend(rawValue: args[4]) {
            let receipt = try VivoH5ADCountStore.normalize(URL(fileURLWithPath: args[2]), target: target, backend: backend,
                implementation: implementation, to: URL(fileURLWithPath: args[5]))
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt)); return
        }
        if args.count == 4, args[1] == "verify-normalized-count-store" {
            _ = try VivoH5ADCountStore.verifyNormalization(URL(fileURLWithPath: args[2]), store: URL(fileURLWithPath: args[3]), implementation: implementation)
            print("verified"); return
        }
        if args.count == 5, args[1] == "multiassay-h5mu-import" {
            let plan = try JSONDecoder().decode(VivoH5MUMultiAssayPlan.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _ = try VivoMultiAssayIO.importH5MU(source: URL(fileURLWithPath: args[2]), plan: plan, implementation: implementation, to: URL(fileURLWithPath: args[4])); print("imported"); return
        }
        if args.count == 4, args[1] == "multiassay-write" {
            let data = try JSONDecoder().decode(VivoMultiAssayDataset.self, from: Data(contentsOf: URL(fileURLWithPath: args[2])))
            try VivoMultiAssayH5MU.write(data, to: URL(fileURLWithPath: args[3])); print("exported"); return
        }
        if args.count == 5, args[1] == "multiassay-import" {
            let plan = try JSONDecoder().decode(VivoTenXMultiAssayPlan.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _ = try VivoMultiAssayIO.importTenX(source: URL(fileURLWithPath: args[2]), plan: plan, implementation: implementation, to: URL(fileURLWithPath: args[4])); print("imported"); return
        }
        if args.count == 3, args[1] == "multiassay-verify" {
            _ = try VivoMultiAssayIO.verify(URL(fileURLWithPath: args[2]), implementation: implementation); print("verified"); return
        }
        if args.count==5,args[1]=="composition-prepare" {
            let plan=try JSONDecoder().decode(VivoCompositionPreparation.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            try VivoComposition.prepare(from: URL(fileURLWithPath: args[2]),plan: plan,implementation: implementation,to: URL(fileURLWithPath: args[4]));print("prepared");return
        }
        if args.count==4,args[1]=="composition-fit" {
            _=try VivoComposition.fit(source: URL(fileURLWithPath: args[2]),implementation: implementation,to: URL(fileURLWithPath: args[3]));print("fit");return
        }
        if args.count==5,args[1]=="composition-predict" {
            let plan=try JSONDecoder().decode(VivoCompositionQueryPlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _=try VivoComposition.predict(reference: URL(fileURLWithPath: args[2]),plan: plan,implementation: implementation,to: URL(fileURLWithPath: args[4]));print("predicted");return
        }
        if args.count==3,args[1]=="composition-verify" {
            _=try VivoComposition.verifyModel(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count==3,args[1]=="composition-prediction-verify" {
            _=try VivoComposition.verifyPrediction(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count==5,args[1]=="perturbation-fit" {
            let plan=try JSONDecoder().decode(VivoPerturbationPlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _=try VivoPerturbation.fit(source: URL(fileURLWithPath: args[2]),plan: plan,implementation: implementation,to: URL(fileURLWithPath: args[4]));print("fit");return
        }
        if args.count==6,args[1]=="perturbation-predict" {
            let plan=try JSONDecoder().decode(VivoPerturbationQueryPlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _=try VivoPerturbation.map(source: URL(fileURLWithPath: args[2]),plan: plan,reference: URL(fileURLWithPath: args[4]),implementation: implementation,to: URL(fileURLWithPath: args[5]));print("mapped");return
        }
        if args.count==3,args[1]=="perturbation-verify" {
            _=try VivoPerturbation.verifyModel(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count==3,args[1]=="perturbation-prediction-verify" {
            _=try VivoPerturbation.verifyPrediction(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count==5,args[1]=="reference-fit" {
            let plan=try JSONDecoder().decode(VivoSingleCellReferencePlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _=try VivoSingleCellReference.fit(source: URL(fileURLWithPath: args[2]),plan: plan,implementation: implementation,to: URL(fileURLWithPath: args[4]));print("fit");return
        }
        if args.count==6,args[1]=="reference-map" {
            let plan=try JSONDecoder().decode(VivoSingleCellReferenceQueryPlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            _=try VivoSingleCellReference.map(source: URL(fileURLWithPath: args[2]),plan: plan,reference: URL(fileURLWithPath: args[4]),implementation: implementation,to: URL(fileURLWithPath: args[5]));print("mapped");return
        }
        if args.count==3,args[1]=="reference-verify" {
            _=try VivoSingleCellReference.verifyReference(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count==3,args[1]=="reference-map-verify" {
            _=try VivoSingleCellReference.verifyMapping(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count == 5,args[1] == "project" {
            let plan=try JSONDecoder().decode(VivoH5ADProjectionPlan.self,from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            let implementation=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
            let receipt=try VivoH5ADProjection.publish(source: URL(fileURLWithPath: args[2]),plan: plan,implementation: implementation,to: URL(fileURLWithPath: args[4]))
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt));return
        }
        if args.count == 3,args[1] == "verify-project" {
            let implementation=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
            _=try VivoH5ADProjection.verify(URL(fileURLWithPath: args[2]),implementation: implementation);print("verified");return
        }
        if args.count == 5, args[1] == "aggregate" {
            let plan = try JSONDecoder().decode(VivoH5ADPseudobulkPlan.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            let implementation = try VivoFingerprint(bytes: Array(repeating: 0, count: 32))
            let receipt = try VivoH5ADPseudobulk.publish(source: URL(fileURLWithPath: args[2]), plan: plan, implementation: implementation, to: URL(fileURLWithPath: args[4]))
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt))
            return
        }
        if args.count == 3, args[1] == "verify-aggregate" {
            let implementation = try VivoFingerprint(bytes: Array(repeating: 0, count: 32))
            _ = try VivoH5ADPseudobulk.verify(URL(fileURLWithPath: args[2]), implementation: implementation)
            print("verified")
            return
        }
        if args.count == 5, args[1] == "annotate" {
            let plan = try JSONDecoder().decode(VivoH5ADAnnotationPlan.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            let implementation = try VivoFingerprint(bytes: Array(repeating: 0, count: 32))
            let receipt = try VivoSingleCellH5AD.annotate(URL(fileURLWithPath: args[2]), plan: plan, implementation: implementation,
                                                        to: URL(fileURLWithPath: args[4]))
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt))
            return
        }
        guard args.count == 4 else { throw VivoOmicsError.invalid("expected source, plan and output arguments") }
        let source = URL(fileURLWithPath: args[1]), planURL = URL(fileURLWithPath: args[2])
        let plan = try JSONDecoder().decode(VivoH5ADImportPlan.self, from: Data(contentsOf: planURL))
        let document = try VivoSingleCellH5AD.read(source, plan: plan)
        let out = URL(fileURLWithPath: args[3])
        try FileManager.default.createDirectory(at: out, withIntermediateDirectories: false)
        try VivoCanonicalJSON.encode(["hdf5Version": document.hdf5Version]).write(to: out.appendingPathComponent("runtime.json"))
        try VivoCanonicalJSON.encode(document.dataset).write(to: out.appendingPathComponent("dataset.json"))
        try VivoCanonicalJSON.encode(VivoSingleCellAnalysis.quality(document.dataset)).write(to: out.appendingPathComponent("quality.json"))
        try VivoCanonicalJSON.encode(VivoSingleCellAnalysis.pseudobulk(document.dataset)).write(to: out.appendingPathComponent("pseudobulk.json"))
        try document.exportOriginal(to: out.appendingPathComponent("original.h5ad"))
        try VivoSingleCellH5AD.write(document.dataset, to: out.appendingPathComponent("native.h5ad"))
        do {
            try VivoSingleCellH5AD.write(document.dataset, to: out.appendingPathComponent("native.h5ad"))
            fatalError("overwrite accepted")
        } catch {}
    }
}
