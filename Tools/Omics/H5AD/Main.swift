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
        if args.count == 5, args[1] == "annotate" {
            let plan = try JSONDecoder().decode(VivoH5ADAnnotationPlan.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
            let implementation = try VivoFingerprint(bytes: Array(repeating: 0, count: 32))
            let receipt = try VivoSingleCellH5AD.annotate(URL(fileURLWithPath: args[2]), plan: plan, implementation: implementation,
                                                        to: URL(fileURLWithPath: args[4]))
            FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(receipt))
            return
        }
        let source = URL(fileURLWithPath: args[1]), planURL = URL(fileURLWithPath: args[2])
        let plan = try JSONDecoder().decode(VivoH5ADImportPlan.self, from: Data(contentsOf: planURL))
        let document = try VivoSingleCellH5AD.read(source, plan: plan)
        let out = URL(fileURLWithPath: args[3])
        try FileManager.default.createDirectory(at: out, withIntermediateDirectories: false)
        try VivoCanonicalJSON.encode(["hdf5Version": document.hdf5Version]).write(to: out.appendingPathComponent("runtime.json"))
        try VivoCanonicalJSON.encode(document.dataset).write(to: out.appendingPathComponent("dataset.json"))
        try document.exportOriginal(to: out.appendingPathComponent("original.h5ad"))
        try VivoSingleCellH5AD.write(document.dataset, to: out.appendingPathComponent("native.h5ad"))
        do {
            try VivoSingleCellH5AD.write(document.dataset, to: out.appendingPathComponent("native.h5ad"))
            fatalError("overwrite accepted")
        } catch {}
    }
}
