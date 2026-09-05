import Foundation
import NumiVivoKit

struct VivoReactionCLICommands {
    static func handles(_ name:String?) -> Bool {["reaction-template","reaction-run","reaction-help"].contains(name ?? "")}
    private struct Receipt:Codable {let request:VivoFingerprint;let task:VivoFingerprint;let receipt:VivoFingerprint;let result:VivoFingerprint;let reused:Bool}
    private struct Implementation:Codable {let executable:VivoFingerprint;let platform:String;let architecture:String}
    private func read(_ url:URL,limit:Int) throws -> Data {
        let m=try url.resourceValues(forKeys:[.isRegularFileKey,.fileSizeKey])
        guard m.isRegularFile==true,let size=m.fileSize,size>=0,size<=limit else {throw VivoChemistryError.resourceLimit("reaction input file bound")}
        let data=try Data(contentsOf:url,options:.mappedIfSafe)
        guard data.count<=limit else {throw VivoChemistryError.resourceLimit("reaction input grew past limit")};return data
    }
    private func write(_ data:Data,to url:URL?) throws {
        if let url {try FileManager.default.createDirectory(at:url.deletingLastPathComponent(),withIntermediateDirectories:true);try data.write(to:url,options:.atomic)}
        else {FileHandle.standardOutput.write(data);FileHandle.standardOutput.write(Data("\n".utf8))}
    }
    func run(arguments:[String]) async -> Int32 {
        do {
            guard let command=arguments.first,Self.handles(command) else {throw VivoChemistryError.invalid("reaction command")}
            if command=="reaction-help" {
                guard arguments.count==1 else {throw VivoChemistryError.invalid("reaction-help has no options")}
                FileHandle.standardOutput.write(Data(Self.help.utf8));return 0
            }
            guard arguments.count>=2,!arguments[1].hasPrefix("--") else {throw VivoChemistryError.invalid("reaction template name or request path required")}
            var options:[String:String]=[:],i=2
            while i<arguments.count {
                let key=arguments[i]
                guard ["--output","--store"].contains(key),options[key]==nil,i+1<arguments.count,!arguments[i+1].hasPrefix("--") else {throw VivoChemistryError.invalid("unknown or duplicate reaction option")}
                options[key]=arguments[i+1];i+=2
            }
            let output=options["--output"].map{URL(fileURLWithPath:$0)}
            if command=="reaction-template" {
                guard options["--store"]==nil else {throw VivoChemistryError.invalid("template does not use an artifact store")}
                try write(VivoCanonicalJSON.encode(VivoReactionQualificationWorkflow.template(arguments[1])),to:output);return 0
            }
            let source=URL(fileURLWithPath:arguments[1]).resolvingSymlinksInPath().standardizedFileURL
            let root=URL(fileURLWithPath:options["--store"] ?? ".numivivo/chemistry-artifacts").resolvingSymlinksInPath().standardizedFileURL
            if let output {for url in [output,URL(fileURLWithPath:output.path+".receipt.json")] {
                let target=url.resolvingSymlinksInPath().standardizedFileURL
                guard target != source,target != root,!target.path.hasPrefix(root.path+"/") else {throw VivoChemistryError.invalid("reaction output aliases source or artifact store")}
                if let a=try? FileManager.default.attributesOfItem(atPath:source.path),let b=try? FileManager.default.attributesOfItem(atPath:target.path),
                   let x=a[.systemFileNumber] as? NSNumber,let y=b[.systemFileNumber] as? NSNumber,
                   let dx=a[.systemNumber] as? NSNumber,let dy=b[.systemNumber] as? NSNumber,x==y,dx==dy {throw VivoChemistryError.invalid("hard-link reaction output alias")}
            }}
            let request=try VivoCanonicalJSON.decode(VivoReactionCalculationRequest.self,from:read(source,limit:256*1024*1024))
            try request.validate()
            guard let executable=Bundle.main.executableURL else {throw VivoChemistryError.invalid("executing binary identity unavailable")}
            #if arch(arm64)
            let arch="arm64"
            #elseif arch(x86_64)
            let arch="x86_64"
            #else
            let arch="other"
            #endif
            let identity=Implementation(executable:try VivoCanonicalJSON.fingerprint(read(executable,limit:1024*1024*1024)),platform:ProcessInfo.processInfo.operatingSystemVersionString,architecture:arch)
            let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity)),store=try VivoArtifactStore(rootURL:root)
            let input=try await store.put(data:VivoCanonicalJSON.encode(request),kind:"vivo.reaction-calculation-request",mediaType:"application/json")
            let operation=VivoReactionQualificationWorkflow.operation(implementationFingerprint:id)
            let task=VivoChemistryTask(operation:operation.identifier,version:operation.version,implementationFingerprint:id,
                inputs:[.init(name:"request",artifact:input.fingerprint,kind:input.kind)],configuration:.object([:]),outputs:operation.outputs,
                resources:.init(budget:request.budget,maximumInputBytes:request.budget.maximumBytes,maximumOutputBytes:request.budget.maximumBytes))
            let workflow=VivoChemistryWorkflow(store:store),result=try await workflow.run(task,using:operation)
            guard let artifact=result.outputs.first else {throw VivoChemistryError.invalid("reaction workflow omitted its result")}
            try write(await workflow.payload(artifact:artifact.artifact,expectedKind:artifact.kind),to:output)
            if let output {try write(VivoCanonicalJSON.encode(Receipt(request:input.fingerprint,task:result.taskFingerprint,receipt:result.receiptFingerprint,result:artifact.artifact,reused:result.reused)),to:URL(fileURLWithPath:output.path+".receipt.json"))}
            FileHandle.standardError.write(Data("Reaction calculation stored: \(artifact.artifact.hex), reused=\(result.reused). No kinetic parameter installed.\n".utf8))
            return 0
        } catch {FileHandle.standardError.write(Data("Reaction calculation failed: \(error)\n".utf8));return 1}
    }
    private static let help="""
    Native nuclear and harmonic reaction qualification
      numivivo reaction-template h3-saddle --output saddle.json
      numivivo reaction-run saddle.json --output saddle.result.json --store .numivivo/chemistry-artifacts
    Templates: h2-minimum, h3-saddle, h-atom, h2-solvated-path.
    Request calculations: qualify, solvatedPath, harmonicBarrier, descent.
    Coordinates: Bohr. Masses: explicit Da. Energies: Hartree. RRHO requires all
    free-molecule vibrational modes resolved and a stationary nuclear gradient.
    A local first-order saddle is not a demonstrated reaction connection. A
    harmonic Gibbs barrier estimate does not automatically become a kinetic rate.
    Correlated smooth-CPCM coupling uses an explicitly frozen RHF-reference field.
    H3/H2/STO-3G templates check numerical execution; they are not the paper reaction.

    """
}
