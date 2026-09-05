import Foundation
import NumiVivoKit

struct VivoElectronicCLICommands {
    static func handles(_ command:String?) -> Bool {
        ["chemistry-template","chemistry-run","chemistry-solve","chemistry-help"].contains(command ?? "")
    }
    private struct ImplementationIdentity: Codable {
        let executable: VivoFingerprint
        let operatingSystem: String
        let architecture: String
        let numericalBackend: String
    }
    private struct NodeReceipt: Codable {
        let identifier: String
        let task: VivoFingerprint
        let receipt: VivoFingerprint
        let reused: Bool
    }
    private struct RunReceipt: Codable {
        let schema: String
        let input: VivoFingerprint
        let implementation: VivoFingerprint
        let nodes: [NodeReceipt]
        let result: VivoFingerprint
        let resultKind: String
    }
    private struct Arguments {
        let command: String
        let positional: [String]
        let options: [String:String]
        init(_ arguments:[String]) throws {
            command=arguments.first ?? "chemistry-help"
            var positional:[String]=[], options:[String:String]=[:], i=1
            let allowed:Set<String> = ["--output","--store","--solver","--budget"]
            while i<arguments.count {
                let word=arguments[i]
                if word.hasPrefix("--") {
                    guard allowed.contains(word), options[word]==nil, i+1<arguments.count,
                          !arguments[i+1].hasPrefix("--") else { throw VivoChemistryError.invalid("unknown, duplicate or valueless option \(word)") }
                    options[word]=arguments[i+1]; i+=2
                } else { positional.append(word); i+=1 }
            }
            self.positional=positional; self.options=options
        }
        func allow(_ options:Set<String>) throws {
            guard Set(self.options.keys).isSubset(of:options) else { throw VivoChemistryError.invalid("option does not apply to \(command)") }
        }
    }
    private func read(_ path:String,maximumBytes:Int = 256*1024*1024) throws -> Data {
        let url=URL(fileURLWithPath:path)
        let values=try url.resourceValues(forKeys:[.isRegularFileKey,.fileSizeKey])
        guard values.isRegularFile==true, let size=values.fileSize, size>=0, size<=maximumBytes else {
            throw VivoChemistryError.resourceLimit("input must be a regular file within the byte budget: \(path)")
        }
        let data=try Data(contentsOf:url,options:.mappedIfSafe)
        guard data.count<=maximumBytes else { throw VivoChemistryError.resourceLimit("input grew beyond its byte budget") }
        return data
    }
    private func implementation() throws -> VivoFingerprint {
        let path:String
        if let executable=Bundle.main.executableURL { path=executable.path }
        else {
            let argument=CommandLine.arguments[0]
            if argument.contains("/") { path=URL(fileURLWithPath:argument).standardizedFileURL.path }
            else {
                let directories=(ProcessInfo.processInfo.environment["PATH"] ?? "").split(separator:":")
                guard let resolved=directories.map({URL(fileURLWithPath:String($0)).appendingPathComponent(argument).path})
                    .first(where:{FileManager.default.isExecutableFile(atPath:$0)}) else {
                    throw VivoChemistryError.invalid("cannot resolve the executing binary for cache identity")
                }
                path=resolved
            }
        }
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
        let binary=try VivoCanonicalJSON.fingerprint(read(path,maximumBytes:1024*1024*1024))
        let identity=ImplementationIdentity(executable:binary,operatingSystem:ProcessInfo.processInfo.operatingSystemVersionString,
                                          architecture:architecture,numericalBackend:backend)
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity))
    }
    private func put<T:Encodable>(_ value:T,kind:String,store:VivoArtifactStore) async throws -> VivoFingerprint {
        let artifact=try await store.put(data:VivoCanonicalJSON.encode(value),kind:kind,mediaType:"application/json")
        return artifact.fingerprint
    }
    private func output(_ data:Data,path:String?) throws {
        if let path {
            let url=URL(fileURLWithPath:path)
            try FileManager.default.createDirectory(at:url.deletingLastPathComponent(),withIntermediateDirectories:true)
            try data.write(to:url,options:.atomic)
        } else { FileHandle.standardOutput.write(data); FileHandle.standardOutput.write(Data("\n".utf8)) }
    }
    private func finish(payload:Data,receipt:RunReceipt,args:Arguments) throws {
        try output(payload,path:args.options["--output"])
        if let path=args.options["--output"] { try output(VivoCanonicalJSON.encode(receipt),path:path+".receipt.json") }
        let reused=receipt.nodes.filter(\.reused).count
        FileHandle.standardError.write(Data("Native chemistry complete: \(receipt.nodes.count) stages, \(reused) reused; result SHA-256 \(receipt.result.hex).\n".utf8))
    }
    func run(arguments:[String]) async -> Int32 {
        do {
            let args=try Arguments(arguments)
            if args.command == "chemistry-help" {
                guard args.positional.isEmpty, args.options.isEmpty else { throw VivoChemistryError.invalid("chemistry-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            if args.command == "chemistry-template" {
                try args.allow(["--output"])
                guard args.positional.count==1 else { throw VivoChemistryError.invalid("chemistry-template requires one template name") }
                try output(template(args.positional[0]),path:args.options["--output"]); return 0
            }
            guard args.positional.count==1 else { throw VivoChemistryError.invalid("one input JSON file is required") }
            let source=args.positional[0]
            switch args.command {
            case "chemistry-run": try args.allow(["--output","--store"])
            case "chemistry-solve":
                try args.allow(["--output","--store","--solver","--budget"])
                guard args.options["--solver"] != nil else { throw VivoChemistryError.invalid("chemistry-solve requires --solver solver.json") }
            default: throw VivoChemistryError.invalid("unknown electronic command")
            }
            if let destination=args.options["--output"] {
                let inputs=([source]+[args.options["--solver"],args.options["--budget"]].compactMap{$0}).map {
                    URL(fileURLWithPath:$0).resolvingSymlinksInPath().standardizedFileURL
                }
                for path in [destination,destination+".receipt.json"] {
                    guard !inputs.contains(URL(fileURLWithPath:path).resolvingSymlinksInPath().standardizedFileURL) else {
                        throw VivoChemistryError.invalid("output/receipt must not replace any input request, solver or budget file")
                    }
                }
            }
            let store=try VivoArtifactStore(rootURL:URL(fileURLWithPath:args.options["--store"] ?? ".numivivo/chemistry-artifacts"))
            let workflow=VivoChemistryWorkflow(store:store), identity=try implementation()
            switch args.command {
            case "chemistry-run":
                try args.allow(["--output","--store"])
                let request=try VivoElectronicRequestDocument.decode(read(source))
                let stored=try await store.put(data:request.canonicalData(),kind:request.artifactKind,mediaType:"application/json")
                let requestID=stored.fingerprint
                let system=try await put(request.system,kind:"vivo.electronic-system",store:store)
                let basis=try await put(request.basis,kind:"vivo.gaussian-basis",store:store)
                let plan=try request.plan(systemArtifact:system,basisArtifact:basis,implementationFingerprint:identity)
                let results=try await workflow.runDAG(plan.nodes)
                guard let result=results[plan.resultNode]?.outputs.first(where:{$0.name==plan.resultOutput}) else { throw VivoChemistryError.invalid("workflow omitted the requested result") }
                let payload=try await workflow.payload(artifact:result.artifact,expectedKind:plan.resultKind)
                let records=results.keys.sorted().map { key in
                    NodeReceipt(identifier:key,task:results[key]!.taskFingerprint,receipt:results[key]!.receiptFingerprint,reused:results[key]!.reused)
                }
                try finish(payload:payload,receipt:.init(schema:"numivivo.org/electronic-cli-receipt/v1",input:requestID,
                    implementation:identity,nodes:records,result:result.artifact,resultKind:result.kind),args:args)
            case "chemistry-solve":
                try args.allow(["--output","--store","--solver","--budget"])
                guard let solverPath=args.options["--solver"] else { throw VivoChemistryError.invalid("chemistry-solve requires --solver solver.json") }
                let budget:VivoChemistryBudget
                if let budgetPath=args.options["--budget"] { budget=try VivoCanonicalJSON.decode(VivoChemistryBudget.self,from:read(budgetPath)) }
                else { budget = .init() }
                try budget.validate()
                let h=try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from:read(source,maximumBytes:budget.maximumBytes))
                let solver=try VivoNativeSolverDispatch.decode(read(solverPath),implementationFingerprint:identity)
                try h.validate(budget:budget)
                let inputID=try await put(h,kind:"vivo.embedded-hamiltonian",store:store)
                let operation=solver.operation
                let config=solver.configuration
                let task=VivoChemistryTask(operation:operation.identifier,version:operation.version,implementationFingerprint:identity,
                    inputs:[.init(name:"hamiltonian",artifact:inputID,kind:"vivo.embedded-hamiltonian")],configuration:config,
                    outputs:operation.outputs,resources:.init(budget:budget,maximumInputBytes:budget.maximumBytes,maximumOutputBytes:budget.maximumBytes))
                let result=try await workflow.run(task,using:operation)
                guard let value=result.outputs.first else { throw VivoChemistryError.invalid("many-body result missing") }
                let payload=try await workflow.payload(artifact:value.artifact,expectedKind:value.kind)
                try finish(payload:payload,receipt:.init(schema:"numivivo.org/electronic-cli-receipt/v1",input:inputID,
                    implementation:identity,nodes:[.init(identifier:"many-body",task:result.taskFingerprint,receipt:result.receiptFingerprint,reused:result.reused)],
                    result:value.artifact,resultKind:value.kind),args:args)
            default: throw VivoChemistryError.invalid("unknown electronic command")
            }
            return 0
        } catch {
            FileHandle.standardError.write(Data("Native chemistry failed: \(error)\n".utf8)); return 1
        }
    }
    private func template(_ name:String) throws -> Data {
        if name == "solver-ccsd" { return try VivoCanonicalJSON.encode(VivoManyBodySolverRequest.ccsd(configuration:.init())) }
        if name == "solver-lih-sa-casscf" {
            return try VivoCanonicalJSON.encode(VivoAdvancedManyBodyRequest.multistateCASSCF(partition:.init(doublyOccupiedCore:[0],active:[1,2]),
                configuration:.init(weights:[0.5,0.5],rootLabels:["ground","excited"],optimization:.init(maximumMacroIterations:200,maximumEnergyEvaluations:4000,gradientTolerance:1e-6),followRoots:true)))
        }
        if name == "solver-tensor-ccsd" {return try VivoCanonicalJSON.encode(VivoAdvancedManyBodyRequest.tensorCCSD(configuration:.init()))}
        if name == "solver-direct-fci" {return try VivoCanonicalJSON.encode(VivoAdvancedManyBodyRequest.directCI(configuration:.init()))}
        if name == "solver-fci" { return try VivoCanonicalJSON.encode(VivoManyBodySolverRequest.configurationInteraction(method:.fci)) }
        let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),
            .init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
        let ecc=VivoECCDMETConfiguration(mode:.singleFragment,fragments:[.init(identifier:"hydrogen-fragment",orbitals:[0],maximumBathOrbitals:1,clusterAlphaElectrons:1,clusterBetaElectrons:1)],bathSelection:.init(minimumBathOrbitals:1))
        if name == "solver-ecc-dmet" {return try VivoCanonicalJSON.encode(VivoECCDMETSolverRequest(configuration:ecc))}
        let advanced:VivoAdvancedElectronicCalculation?
        switch name {
        case "h2-direct-fci":advanced = .correlated(reference:.init(),solver:.directCI(configuration:.init()))
        case "h2-tensor-ccsd":advanced = .correlated(reference:.init(),solver:.tensorCCSD(configuration:.init()))
        case "h2-smooth-cpcm":advanced = .smoothCPCM(configuration:.init())
        case "h2-ecc-dmet":advanced = .eccDMET(configuration:ecc)
        default:advanced=nil
        }
        if let advanced{return try VivoCanonicalJSON.encode(VivoAdvancedElectronicWorkflowRequest(system:system,basis:.hydrogenSTO3G(nucleusIndices:[0,1]),calculation:advanced))}
        let calculation:VivoElectronicWorkflowCalculation
        switch name {
        case "h2-ccsd": calculation = .correlated(reference:.init(),solver:.ccsd(configuration:.init()))
        case "h2-casci": calculation = .correlated(reference:.init(),solver:.casci(partition:.init(active:[0,1])))
        case "h2-casscf": calculation = .correlated(reference:.init(),solver:.casscf(partition:.init(active:[0,1]),configuration:.init()))
        case "h2-lda": calculation = .lda(configuration:.init())
        case "h2-lda-cpcm": calculation = .lda(configuration:.init(cpcm:.init()))
        case "h2-cpcm-rhf": calculation = .cpcmRHF(configuration:.init())
        default: throw VivoChemistryError.invalid("unknown template; see chemistry-help")
        }
        return try VivoCanonicalJSON.encode(VivoElectronicWorkflowRequest(system:system,
            basis:.hydrogenSTO3G(nucleusIndices:[0,1]),calculation:calculation))
    }
    private static let help = """
    Native electronic chemistry (FP64; no Python or CUDA runtime)

      numivivo chemistry-template h2-ccsd --output h2.json
      numivivo chemistry-run h2.json --store .numivivo/chemistry-artifacts --output h2.result.json
      numivivo chemistry-template solver-fci --output solver.json
      numivivo chemistry-solve hamiltonian.json --solver solver.json --output result.json

    Templates: h2-ccsd, h2-casci, h2-casscf, h2-lda, h2-lda-cpcm,
               h2-cpcm-rhf, h2-direct-fci, h2-tensor-ccsd, h2-smooth-cpcm,
               h2-ecc-dmet, solver-ccsd, solver-fci, solver-direct-fci,
               solver-tensor-ccsd, solver-ecc-dmet, solver-lih-sa-casscf.
    chemistry-run explicitly dispatches original and advanced request schemas.
    Advanced JSON supports densityFittedMP2 (with an explicit auxiliary basis),
    multistateCASSCF, smoothCPCM and integrated eccDMET. No method fallback.
    chemistry-solve accepts an optional --budget budget.json.
    With --output, a sibling <output>.receipt.json records input/result hashes,
    implementation identity and reused stages. Nonconverged solvers fail and
    do not publish a successful workflow receipt. All coordinates are in Bohr;
    these electronic/continuum energies are not activation Gibbs free energies.

    """
}
