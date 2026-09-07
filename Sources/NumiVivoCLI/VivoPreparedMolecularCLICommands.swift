import Foundation
import NumiVivoKit

struct VivoPreparedMolecularCLICommands {
    static func handles(_ name:String?) -> Bool {
        ["molecule-prepare-template","molecule-prepare","molecule-sampling-template","molecule-sampling-analyze",
         "molecule-sampling-run","molecule-rate","molecule-help"].contains(name ?? "") || VivoMolecularSamplingExportCLICommands.handles(name)
    }
    private struct Implementation:Codable {let executable:VivoFingerprint;let platform:String;let architecture:String}
    private struct Receipt:Codable {let request:VivoFingerprint;let task:VivoFingerprint;let receipt:VivoFingerprint;let result:VivoFingerprint;let reused:Bool}
    private func load<T:Decodable & Sendable>(_ type:T.Type,_ path:String) throws -> T {
        try VivoValidatedArtifactLoader.decode(type,at:URL(fileURLWithPath:path),limits:.init(
            maximumBytes:256*1024*1024,maximumNodes:8_000_000,maximumArrayElements:2_000_000)).value
    }
    private func canonical(_ url:URL) -> URL {
        var root=url.standardizedFileURL,suffix:[String]=[]
        while !FileManager.default.fileExists(atPath:root.path),root.path != "/" {
            suffix.append(root.lastPathComponent);root.deleteLastPathComponent()
        }
        var result=root.resolvingSymlinksInPath().standardizedFileURL
        for part in suffix.reversed() {result.appendPathComponent(part)}
        return result.standardizedFileURL
    }
    private func protect(_ output:URL,inputs:[URL],store:URL) throws {
        let target=canonical(output),root=canonical(store),prefix=root.path == "/" ? "/":root.path+"/"
        guard target != root,!target.path.hasPrefix(prefix) else {throw VivoChemistryError.invalid("molecular output aliases artifact storage")}
        for input in inputs {
            guard target != canonical(input) else {throw VivoChemistryError.invalid("molecular output aliases an input")}
            if let a=try? FileManager.default.attributesOfItem(atPath:target.path),let b=try? FileManager.default.attributesOfItem(atPath:input.path),
               let ai=a[.systemFileNumber] as? NSNumber,let bi=b[.systemFileNumber] as? NSNumber,
               let ad=a[.systemNumber] as? NSNumber,let bd=b[.systemNumber] as? NSNumber,ai==bi,ad==bd {
                throw VivoChemistryError.invalid("molecular output is hard-linked to an input")
            }
        }
    }
    private func write<T:Encodable>(_ value:T,to path:String?) throws {
        let data=try VivoCanonicalJSON.encode(value)
        guard let path else {FileHandle.standardOutput.write(data);FileHandle.standardOutput.write(Data("\n".utf8));return}
        // Existing pinned-root, no-clobber writer: reruns do not destroy results.
        try VivoKineticsDocumentIO.write(data,to:URL(fileURLWithPath:path),overwrite:false)
    }
    private func fingerprint(_ text:String) throws -> VivoFingerprint {
        let chars=Array(text.utf8)
        guard chars.count==64,chars.allSatisfy({(48...57).contains($0) || (97...102).contains($0)}) else {
            throw VivoChemistryError.invalid("resume requires a lowercase SHA-256 checkpoint hash")
        }
        func digit(_ b:UInt8)->UInt8 {b<=57 ? b-48:b-97+10}
        return try .init(bytes:stride(from:0,to:64,by:2).map{digit(chars[$0])*16+digit(chars[$0+1])})
    }
    func run(arguments:[String]) async -> Int32 {
        if VivoMolecularSamplingExportCLICommands.handles(arguments.first) {
            return await VivoMolecularSamplingExportCLICommands().run(arguments: arguments)
        }
        do {
            guard let command=arguments.first,Self.handles(command) else {throw VivoChemistryError.invalid("molecular command")}
            if command=="molecule-help" {
                guard arguments.count==1 else {throw VivoChemistryError.invalid("molecule-help has no options")}
                FileHandle.standardOutput.write(Data(Self.help.utf8));return 0
            }
            guard arguments.count>=2,!arguments[1].hasPrefix("--") else {throw VivoChemistryError.invalid("molecular input JSON is required")}
            var options:[String:String]=[:],i=2
            while i<arguments.count {
                let key=arguments[i]
                guard key.hasPrefix("--"),options[key]==nil,i+1<arguments.count,!arguments[i+1].hasPrefix("--") else {
                    throw VivoChemistryError.invalid("duplicate or malformed molecular option")
                }
                options[key]=arguments[i+1];i+=2
            }
            var allowed:Set<String>=["--output","--store"]
            switch command {
            case "molecule-prepare-template":allowed.formUnion(["--microstate","--pH","--protonation-source","--forcefield"])
            case "molecule-sampling-template":allowed.formUnion(["--context","--distance","--temperature","--cutoff"])
            case "molecule-sampling-run":allowed.formUnion(["--resume", "--read-limits"])
            case "molecule-rate":allowed.insert("--kinetics")
            default:break
            }
            guard Set(options.keys).isSubset(of:allowed) else {throw VivoChemistryError.invalid("unknown molecular command option")}
            guard options["--read-limits"] == nil || options["--resume"] != nil else {
                throw VivoChemistryError.invalid("sampling --read-limits applies to --resume")
            }
            let resumeLimits = try VivoMolecularSamplingCLILimits.load(options["--read-limits"])
            guard options["--kinetics"] == nil || options["--output"] != nil else {
                throw VivoChemistryError.invalid("--kinetics requires --output for the separately retained kinetic pack")
            }
            let root=URL(fileURLWithPath:options["--store"] ?? ".numivivo/chemistry-artifacts")
            let inputs=([arguments[1]]+[options["--forcefield"],options["--kinetics"],options["--read-limits"]].compactMap{$0}).map{URL(fileURLWithPath:$0)}
            if let output=options["--output"] {
                for name in [output,output+".receipt.json",output+".kinetics.json"] {try protect(URL(fileURLWithPath:name),inputs:inputs,store:root)}
            }
            if command=="molecule-prepare-template" {
                guard let microstate=options["--microstate"],let source=options["--protonation-source"],
                      let text=options["--pH"],let pH=Double(text),pH.isFinite else {
                    throw VivoChemistryError.invalid("preparation template requires --microstate, --pH and --protonation-source")
                }
                let structure=try load(VivoMolecularStructure.self,arguments[1])
                _ = try VivoStructureValidator.validate(structure)
                let library=try options["--forcefield"].map{try load(VivoForceFieldLibrary.self,$0)}
                let request=VivoMolecularPreparationRequest(structure:structure,microstateIdentifier:microstate,
                    protonationSourceIdentifier:source,pH:pH,expectedFormalCharge:structure.atoms.reduce(0){$0+Int($1.formalCharge)},forceField:library)
                try write(request,to:options["--output"]);return 0
            }
            if command=="molecule-sampling-template" {
                let prepared=try load(VivoMolecularPreparationResult.self,arguments[1])
                guard let compiled=prepared.compiledForceField,!prepared.structure.conformers.isEmpty,
                      let context=options["--context"],let distance=options["--distance"] else {
                    throw VivoChemistryError.invalid("sampling template needs a compiled preparation, --context and --distance atomA,atomB")
                }
                let atoms=distance.split(separator:",").compactMap{UInt32($0)}
                let temperature=Double(options["--temperature"] ?? "298.15"),cutoff=Double(options["--cutoff"] ?? "1.0")
                guard atoms.count==2,let temperature,let cutoff else {throw VivoChemistryError.invalid("sampling distance, temperature or cutoff")}
                let md=VivoMDConfiguration(timeStepPS:0.001,cutoffNM:cutoff,electrostatics:prepared.structure.periodicCell == nil ? .cutoff:.pme,targetTemperatureK:temperature)
                let systemID=try compiled.system.fingerprint()
                let seeds:[UInt64]=[0x4e564d01,0x4e564d02,0x4e564d03,0x4e564d04]
                let states=try seeds.indices.map { i -> VivoClassicalInitialState in
                    let positions=prepared.structure.conformers[i%prepared.structure.conformers.count].positionsNM
                    let state=VivoClassicalInitialState(systemFingerprint:systemID,positionsNM:positions,periodicCell:prepared.structure.periodicCell)
                    try state.validate(particleCount:compiled.system.particles.count);return state
                }
                let request=VivoMolecularSamplingRunRequest(structure:prepared.structure,system:compiled.system,initialStates:states,
                    replicaSeeds:seeds,md:md,contextIdentifier:context,observables:[
                        .init(identifier:"mapped-distance-nm",kind:.distance(atomA:atoms[0],atomB:atoms[1]),maximumMeanStandardError:0.001),
                        .init(identifier:"potential-energy-kj-mol",kind:.potentialEnergyKJPerMol,maximumMeanStandardError:0.1)],
                    minimization:.init(),equilibrationSteps:10000,stepsPerBlock:5000,sampleEvery:100,maximumBlocks:40)
                try request.validate();try write(request,to:options["--output"])
                FileHandle.standardError.write(Data("Sampling template requires system-specific observables and sampling budgets; distinct velocity seeds do not establish conformational diversity.\n".utf8));return 0
            }
            let store=try VivoArtifactStore(rootURL:root)
            if command=="molecule-sampling-run" {
                let request=try load(VivoMolecularSamplingRunRequest.self,arguments[1])
                let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
                FileHandle.standardError.write(Data(("Durable checkpoint reference: "+VivoMolecularSamplingRunner.checkpointReferenceName(requestFingerprint:id)+"\n").utf8))
                let resume=try options["--resume"].map(fingerprint)
                let result=try await VivoMolecularSamplingRunner.run(request,store:store,resumeFrom:resume,resumeReadLimits:resumeLimits)
                _ = try await store.put(data:VivoCanonicalJSON.encode(result),kind:"molecular-sampling-receipt",mediaType:"application/json")
                try write(result,to:options["--output"]);return result.status == .converged ? 0:75
            }
            let calculation:VivoPreparedMolecularCalculation
            switch command {
            case "molecule-prepare":calculation = .prepare(request:try load(VivoMolecularPreparationRequest.self,arguments[1]))
            case "molecule-sampling-analyze":calculation = .sampling(request:try load(VivoMolecularSamplingRequest.self,arguments[1]))
            case "molecule-rate":calculation = .rate(request:try load(VivoPreparedReactionRateRequest.self,arguments[1]),
                kinetics:try options["--kinetics"].map{try load(VivoCovalentKineticPack.self,$0)})
            default:throw VivoChemistryError.invalid("unresolved molecular calculation")
            }
            guard let executable=Bundle.main.executableURL else {throw VivoChemistryError.invalid("executing binary identity unavailable")}
            let meta=try executable.resourceValues(forKeys:[.isRegularFileKey,.fileSizeKey])
            guard meta.isRegularFile==true,let size=meta.fileSize,size>0,size<=1024*1024*1024 else {throw VivoChemistryError.resourceLimit("executable identity byte bound")}
            #if arch(arm64)
            let architecture="arm64"
            #elseif arch(x86_64)
            let architecture="x86_64"
            #else
            let architecture="other"
            #endif
            let binary=try Data(contentsOf:executable,options:.mappedIfSafe)
            guard binary.count<=1024*1024*1024 else {throw VivoChemistryError.resourceLimit("executable identity grew past bound")}
            let implementation=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Implementation(
                executable:VivoCanonicalJSON.fingerprint(binary),platform:ProcessInfo.processInfo.operatingSystemVersionString,architecture:architecture)))
            let input=try await store.put(data:VivoCanonicalJSON.encode(calculation),kind:"vivo.prepared-molecular-request",mediaType:"application/json")
            let operation=VivoPreparedMolecularWorkflow.operation(implementationFingerprint:implementation),budget=calculation.budget
            let task=VivoChemistryTask(operation:operation.identifier,version:operation.version,implementationFingerprint:implementation,
                inputs:[.init(name:"request",artifact:input.fingerprint,kind:input.kind)],configuration:.object([:]),outputs:operation.outputs,
                resources:.init(budget:budget,maximumInputBytes:budget.maximumBytes,maximumOutputBytes:budget.maximumBytes))
            let workflow=VivoChemistryWorkflow(store:store),done=try await workflow.run(task,using:operation)
            guard let artifact=done.outputs.first else {throw VivoChemistryError.invalid("molecular task omitted result")}
            let value=try VivoCanonicalJSON.decode(VivoPreparedMolecularCalculationResult.self,from:await workflow.payload(artifact:artifact.artifact,expectedKind:artifact.kind))
            var exitCode:Int32=0
            switch value {
            case .prepared(let result):try write(result,to:options["--output"])
            case .sampling(let result):try write(result,to:options["--output"]);exitCode=result.converged ? 0:75
            case .rate(let result,let kinetics):
                try write(result,to:options["--output"])
                if let kinetics {
                    guard let output=options["--output"] else {throw VivoChemistryError.invalid("--kinetics requires --output for the updated kinetic pack")}
                    try write(kinetics,to:output+".kinetics.json")
                }
            }
            if let output=options["--output"] {
                try write(Receipt(request:input.fingerprint,task:done.taskFingerprint,receipt:done.receiptFingerprint,
                                  result:artifact.artifact,reused:done.reused),to:output+".receipt.json")
            }
            return exitCode
        } catch {FileHandle.standardError.write(Data("Molecular workflow failed: \(error)\n".utf8));return 65}
    }
    private static let help="""
    Prepared molecular workflows (canonical structures, force fields and artifact store)
      numivivo molecule-prepare-template structure.json --microstate STATE --pH 7.4 --protonation-source SOURCE --forcefield library.json --output preparation.json
      numivivo molecule-prepare preparation.json --output prepared.json
      numivivo molecule-sampling-template prepared.json --context CONTEXT --distance 0,1 --output sampling.json
      numivivo molecule-sampling-run sampling.json --store .numivivo/chemistry-artifacts --output sampling.receipt.json
      numivivo molecule-sampling-run sampling.json --resume CHECKPOINT_SHA256 [--read-limits limits.json] --output resumed.receipt.json
      numivivo molecule-sampling-analyze measurements.json --output diagnostics.json
      numivivo reaction-template h3-connected-rate --output connected.json
      numivivo reaction-run connected.json --output connected.result.json
      numivivo molecule-rate prepared-rate.json --output rate.json [--kinetics kinetic-pack.json]
    Outputs are no-clobber. Pure calculations include a sibling .receipt.json with
    input, task, result and implementation identities. Nonconverged sampling exits 75;
    inspect the receipt rather than treating process completion as convergence.
    Sampling templates are starting specifications, not established molecular sampling
    budgets. Edit observables, initial conformers, timestep, cutoff and run length.
    Checkpoint references survive interruption; --resume takes the referenced hash.
    The connected hydrogen template exercises execution, not realistic barrier accuracy.
    molecule-rate accepts only a mapped, explicitly parameterized, sampled unimolecular
    local model. It retains missing free-energy contributions and an assumed kinetic
    origin. Protein-environment promotion is rejected, not inferred from solvent closure.
    See Documentation/PreparedMolecularWorkflows.md for schemas and scientific limits.

    \(VivoMolecularSamplingExportCLICommands.help)

    """
}
