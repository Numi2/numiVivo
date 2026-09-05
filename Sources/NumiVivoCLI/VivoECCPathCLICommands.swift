import Foundation
import NumiVivoKit

struct VivoECCPathCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["chemistry-path-template","chemistry-path-run","chemistry-path-help"].contains(name ?? "")
    }
    private struct Identity: Codable {
        let executable: VivoFingerprint
        let operatingSystem: String
        let architecture: String
        let numericalBackend: String
    }
    private struct Node: Codable {
        let identifier: String
        let task: VivoFingerprint
        let receipt: VivoFingerprint
        let reused: Bool
    }
    private struct Receipt: Codable {
        let schema: String
        let input: VivoFingerprint
        let implementation: VivoFingerprint
        let nodes: [Node]
        let result: VivoFingerprint
        let resultKind: String
    }
    private func read(_ url: URL, limit: Int = 256*1024*1024) throws -> Data {
        let metadata=try url.resourceValues(forKeys:[.isRegularFileKey,.fileSizeKey])
        guard metadata.isRegularFile==true,let size=metadata.fileSize,size>=0,size<=limit else {
            throw VivoChemistryError.resourceLimit("path input must be a bounded regular file")
        }
        let data=try Data(contentsOf:url,options:.mappedIfSafe)
        guard data.count<=limit else { throw VivoChemistryError.resourceLimit("path input grew beyond its byte limit") }
        return data
    }
    private func write(_ data: Data, _ url: URL?) throws {
        if let url {
            try FileManager.default.createDirectory(at:url.deletingLastPathComponent(),withIntermediateDirectories:true)
            try data.write(to:url,options:.atomic)
        } else { FileHandle.standardOutput.write(data); FileHandle.standardOutput.write(Data("\n".utf8)) }
    }
    /// Same four-field canonical identity as the existing electronic commands.
    /// This preserves primitive task reuse between isolated and path execution.
    private func implementation() throws -> VivoFingerprint {
        guard let executable=Bundle.main.executableURL else { throw VivoChemistryError.invalid("cannot locate running binary for path cache identity") }
        #if arch(arm64)
        let architecture="arm64"
        #elseif arch(x86_64)
        let architecture="x86_64"
        #else
        let architecture="other"
        #endif
        #if canImport(Accelerate)
        let backend="cpu-fp64-accelerate"
        #else
        let backend="cpu-fp64-portable"
        #endif
        let identity=Identity(executable:try VivoCanonicalJSON.fingerprint(read(executable,limit:1024*1024*1024)),
            operatingSystem:ProcessInfo.processInfo.operatingSystemVersionString,architecture:architecture,numericalBackend:backend)
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity))
    }
    private func store<T:Encodable>(_ value: T, kind: String, in artifacts: VivoArtifactStore) async throws -> VivoFingerprint {
        try await artifacts.put(data:VivoCanonicalJSON.encode(value),kind:kind,mediaType:"application/json").fingerprint
    }
    func run(arguments: [String]) async -> Int32 {
        do {
            guard let command=arguments.first,Self.handles(command) else { throw VivoChemistryError.invalid("unknown path command") }
            if command=="chemistry-path-help" {
                guard arguments.count==1 else { throw VivoChemistryError.invalid("path help accepts no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            guard arguments.count>=2,!arguments[1].hasPrefix("--") else { throw VivoChemistryError.invalid("path command requires a template name or request file") }
            var options:[String:String]=[:],index=2
            while index<arguments.count {
                let key=arguments[index]
                guard ["--output","--store"].contains(key),options[key]==nil,index+1<arguments.count,
                      !arguments[index+1].hasPrefix("--") else { throw VivoChemistryError.invalid("unknown, duplicate or valueless path option") }
                options[key]=arguments[index+1]; index+=2
            }
            let output=options["--output"].map { URL(fileURLWithPath:$0) }
            if command=="chemistry-path-template" {
                guard arguments[1]=="h2-stretch",options["--store"]==nil else { throw VivoChemistryError.invalid("path template is h2-stretch; --store is not a template option") }
                try write(VivoCanonicalJSON.encode(VivoMolecularECCPath.hydrogenStretchTemplate()),output)
                return 0
            }
            let source=URL(fileURLWithPath:arguments[1]),storeURL=URL(fileURLWithPath:options["--store"] ?? ".numivivo/chemistry-artifacts")
            if let output {
                let sourceResolved=source.resolvingSymlinksInPath().standardizedFileURL
                let root=storeURL.resolvingSymlinksInPath().standardizedFileURL.path+"/"
                for url in [output,URL(fileURLWithPath:output.path+".receipt.json")] {
                    let resolved=url.resolvingSymlinksInPath().standardizedFileURL
                    guard resolved != sourceResolved,!resolved.path.hasPrefix(root) else {
                        throw VivoChemistryError.invalid("path output/receipt may not replace the request or a store object")
                    }
                    if let a=try? FileManager.default.attributesOfItem(atPath:sourceResolved.path),
                       let b=try? FileManager.default.attributesOfItem(atPath:resolved.path),
                       let inodeA=a[.systemFileNumber] as? NSNumber,let inodeB=b[.systemFileNumber] as? NSNumber,
                       let deviceA=a[.systemNumber] as? NSNumber,let deviceB=b[.systemNumber] as? NSNumber,
                       inodeA==inodeB,deviceA==deviceB {
                        throw VivoChemistryError.invalid("path output is a hard-link alias of the request")
                    }
                }
            }
            let request=try VivoCanonicalJSON.decode(VivoMolecularECCPathRequest.self,from:read(source))
            try request.validate()
            let identity=try implementation(),artifacts=try VivoArtifactStore(rootURL:storeURL)
            let requestID=try await store(request,kind:"vivo.molecular-ecc-path-request",in:artifacts)
            let basisID=try await store(request.basis,kind:"vivo.gaussian-basis",in:artifacts)
            var systemIDs:[VivoFingerprint]=[]
            for snapshot in request.snapshots { systemIDs.append(try await store(snapshot.system,kind:"vivo.electronic-system",in:artifacts)) }
            let plan=try VivoMolecularECCPathWorkflow.plan(request,requestArtifact:requestID,systemArtifacts:systemIDs,
                basisArtifact:basisID,implementationFingerprint:identity)
            let workflow=VivoChemistryWorkflow(store:artifacts),results=try await workflow.runDAG(plan.nodes)
            guard let final=results[plan.resultNode]?.outputs.first(where:{$0.name==plan.resultOutput}) else { throw VivoChemistryError.invalid("path workflow omitted its output") }
            let payload=try await workflow.payload(artifact:final.artifact,expectedKind:plan.resultKind)
            let nodes=results.keys.sorted().map { name in Node(identifier:name,task:results[name]!.taskFingerprint,
                receipt:results[name]!.receiptFingerprint,reused:results[name]!.reused) }
            let receipt=Receipt(schema:"numivivo.org/electronic-cli-receipt/v1",input:requestID,implementation:identity,
                nodes:nodes,result:final.artifact,resultKind:plan.resultKind)
            try write(payload,output)
            if let output { try write(VivoCanonicalJSON.encode(receipt),URL(fileURLWithPath:output.path+".receipt.json")) }
            FileHandle.standardError.write(Data("ECC electronic path: \(request.snapshots.count) geometries, \(nodes.filter(\.reused).count)/\(nodes.count) stages reused; \(final.artifact.hex). Not a characterized activation barrier.\n".utf8))
            return 0
        } catch {
            FileHandle.standardError.write(Data("ECC path failed: \(error)\n".utf8)); return 1
        }
    }
    private static let help="""
    Shared-orbital ECC electronic paths
      numivivo chemistry-path-template h2-stretch --output path.json
      numivivo chemistry-path-run path.json --output profile.json --store .numivivo/chemistry-artifacts

    One ordered atom mapping, basis and electron sector is required. Nuclear coordinates
    use Bohr. External embedding charges are frozen across the path. Physical Gaussian
    cross overlaps align the orbitals; one shared QIO rotation is used at all points.
    Bath dimension, subspace/state continuity, particle closure and shared stationarity
    are checked before a successful artifact is published, including on cache reuse.
    An electronic profile is not a transition-state, Gibbs barrier or rate certificate.
    H2 stretching is a conformance example, not the paper's Michael-addition reaction.

    """
}
