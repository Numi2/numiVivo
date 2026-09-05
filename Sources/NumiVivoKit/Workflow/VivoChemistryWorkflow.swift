import Foundation

public struct VivoChemistryTaskInput: Codable, Sendable, Equatable {
    public let name: String
    public let artifact: VivoFingerprint
    public let kind: String
    public init(name: String, artifact: VivoFingerprint, kind: String) { self.name=name;self.artifact=artifact;self.kind=kind }
}
public struct VivoChemistryTaskOutput: Codable, Sendable, Equatable {
    public let name: String
    public let kind: String
    public init(name: String, kind: String) { self.name=name;self.kind=kind }
}
public struct VivoChemistryResourceContract: Codable, Sendable, Equatable {
    public let numericalBackend: String
    public let budget: VivoChemistryBudget
    public let maximumInputBytes: Int
    public let maximumOutputBytes: Int
    public init(numericalBackend: String = "cpu-fp64", budget: VivoChemistryBudget = .init(),
                maximumInputBytes: Int = 256*1024*1024, maximumOutputBytes: Int = 256*1024*1024) {
        self.numericalBackend=numericalBackend;self.budget=budget
        self.maximumInputBytes=maximumInputBytes;self.maximumOutputBytes=maximumOutputBytes
    }
}
public struct VivoChemistryTask: Codable, Sendable, Equatable {
    public let schema: String
    public let operation: String
    public let version: String
    public let implementationFingerprint: VivoFingerprint
    public let inputs: [VivoChemistryTaskInput]
    public let configuration: VivoJSONValue
    public let outputs: [VivoChemistryTaskOutput]
    public let resources: VivoChemistryResourceContract
    public init(operation: String, version: String, implementationFingerprint: VivoFingerprint,
                inputs: [VivoChemistryTaskInput], configuration: VivoJSONValue,
                outputs: [VivoChemistryTaskOutput], resources: VivoChemistryResourceContract = .init()) {
        schema="numivivo.org/chemistry-task/v1";self.operation=operation;self.version=version
        self.implementationFingerprint=implementationFingerprint;self.inputs=inputs.sorted{$0.name<$1.name}
        self.configuration=configuration;self.outputs=outputs.sorted{$0.name<$1.name};self.resources=resources
    }
    public func fingerprint() throws -> VivoFingerprint {
        try resources.budget.validate()
        guard schema=="numivivo.org/chemistry-task/v1",!operation.isEmpty,!version.isEmpty,!outputs.isEmpty,
              resources.maximumInputBytes>0,resources.maximumOutputBytes>0,!resources.numericalBackend.isEmpty,
              inputs.allSatisfy({!$0.name.isEmpty && !$0.kind.isEmpty}),outputs.allSatisfy({!$0.name.isEmpty && !$0.kind.isEmpty}),
              Set(inputs.map(\.name)).count==inputs.count,Set(outputs.map(\.name)).count==outputs.count else {
            throw VivoChemistryError.invalid("workflow task identity/contracts")
        }
        let canonical=Self(operation:operation,version:version,implementationFingerprint:implementationFingerprint,
                           inputs:inputs,configuration:configuration,outputs:outputs,resources:resources)
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(canonical))
    }
}
public struct VivoChemistryOperation: Sendable {
    public let identifier: String
    public let version: String
    public let implementationFingerprint: VivoFingerprint
    public let numericalBackend: String
    public let outputs: [VivoChemistryTaskOutput]
    public let execute: @Sendable (VivoJSONValue,[String:Data],VivoChemistryBudget) async throws -> [String:Data]
    public let validateOutputs: @Sendable (VivoJSONValue,[String:Data],[String:Data],VivoChemistryBudget) throws -> Void
    public init(identifier: String, version: String, implementationFingerprint: VivoFingerprint,
                numericalBackend: String = "cpu-fp64", outputs: [VivoChemistryTaskOutput],
                execute: @escaping @Sendable (VivoJSONValue,[String:Data],VivoChemistryBudget) async throws -> [String:Data],
                validateOutputs: @escaping @Sendable (VivoJSONValue,[String:Data],[String:Data],VivoChemistryBudget) throws -> Void) {
        self.identifier=identifier;self.version=version;self.implementationFingerprint=implementationFingerprint
        self.numericalBackend=numericalBackend;self.outputs=outputs.sorted{$0.name<$1.name}
        self.execute=execute;self.validateOutputs=validateOutputs
    }
}
private struct VivoChemistryOutputEnvelope: Codable, Sendable {
    let schema: String
    let taskFingerprint: VivoFingerprint
    let outputName: String
    let kind: String
    let payload: Data
}
public struct VivoChemistryOutputReceipt: Codable, Sendable, Equatable {
    public let name: String
    public let kind: String
    public let artifact: VivoFingerprint
}
private struct VivoChemistryTaskReceipt: Codable, Sendable {
    let schema: String
    let taskFingerprint: VivoFingerprint
    let outputs: [VivoChemistryOutputReceipt]
}
public struct VivoChemistryTaskResult: Sendable {
    public let taskFingerprint: VivoFingerprint
    public let receiptFingerprint: VivoFingerprint
    public let outputs: [VivoChemistryOutputReceipt]
    public let reused: Bool
}
public enum VivoChemistryDAGInput: Sendable {
    case artifact(name: String, fingerprint: VivoFingerprint, kind: String)
    case output(name: String, node: String, output: String, kind: String)
}
public struct VivoChemistryDAGNode: Sendable {
    public let identifier: String
    public let operation: VivoChemistryOperation
    public let inputs: [VivoChemistryDAGInput]
    public let configuration: VivoJSONValue
    public let resources: VivoChemistryResourceContract
    public init(identifier: String, operation: VivoChemistryOperation, inputs: [VivoChemistryDAGInput],
                configuration: VivoJSONValue, resources: VivoChemistryResourceContract = .init()) {
        self.identifier=identifier;self.operation=operation;self.inputs=inputs;self.configuration=configuration;self.resources=resources
    }
}
public actor VivoChemistryWorkflow {
    private let store: VivoArtifactStore
    private var inFlight: [VivoFingerprint:Task<VivoChemistryTaskResult,Error>] = [:]
    public init(store: VivoArtifactStore) { self.store=store }
    public func payload(artifact: VivoFingerprint, expectedKind: String) async throws -> Data {
        try await Self.readPayload(store:store,artifact:artifact,expectedKind:expectedKind)
    }
    private static func readPayload(store: VivoArtifactStore, artifact: VivoFingerprint, expectedKind: String) async throws -> Data {
        let descriptor=try await store.descriptor(for:artifact),data=try await store.data(for:artifact,verify:true)
        guard descriptor.byteCount==UInt64(data.count) else { throw VivoChemistryError.invalid("artifact byte count mismatch") }
        if descriptor.kind=="chemistry-output" {
            let envelope=try VivoCanonicalJSON.decode(VivoChemistryOutputEnvelope.self,from:data)
            guard envelope.schema=="numivivo.org/chemistry-output/v1",envelope.kind==expectedKind,!envelope.outputName.isEmpty else {
                throw VivoChemistryError.invalid("chemistry output kind/schema mismatch")
            }
            return envelope.payload
        }
        guard descriptor.kind==expectedKind else { throw VivoChemistryError.invalid("input artifact kind mismatch") };return data
    }
    public func run(_ task: VivoChemistryTask, using operation: VivoChemistryOperation) async throws -> VivoChemistryTaskResult {
        let id=try task.fingerprint()
        guard task.operation==operation.identifier,task.version==operation.version,
              task.implementationFingerprint==operation.implementationFingerprint,
              task.resources.numericalBackend==operation.numericalBackend,
              task.outputs.sorted(by:{$0.name<$1.name})==operation.outputs else {
            throw VivoChemistryError.invalid("operation implementation, backend or output contract mismatch")
        }
        if let running=inFlight[id] { return try await running.value }
        let store=self.store
        let job=Task<VivoChemistryTaskResult,Error> {
            let referenceName="chemistry-task-"+id.hex
            var inputs:[String:Data]=[:],inputBytes=0
            for input in task.inputs {
                let data=try await Self.readPayload(store:store,artifact:input.artifact,expectedKind:input.kind)
                guard data.count<=task.resources.maximumInputBytes-inputBytes else { throw VivoChemistryError.resourceLimit("workflow input byte budget") }
                inputBytes += data.count;inputs[input.name]=data
            }
            let cached: VivoArtifactReference?
            do { cached=try await store.reference(referenceName) }
            catch VivoArtifactStoreError.referenceMissing(_) { cached=nil }
            if let cached {
                guard cached.artifact.kind=="chemistry-task-receipt" else { throw VivoChemistryError.invalid("cache reference kind") }
                let receipt=try VivoCanonicalJSON.decode(VivoChemistryTaskReceipt.self,from:try await store.data(for:cached.artifact.fingerprint,verify:true))
                guard receipt.schema=="numivivo.org/chemistry-task-receipt/v1",receipt.taskFingerprint==id,
                      receipt.outputs.map({VivoChemistryTaskOutput(name:$0.name,kind:$0.kind)}).sorted(by:{$0.name<$1.name})==operation.outputs else {
                    throw VivoChemistryError.invalid("cache receipt identity/output mismatch")
                }
                var payloads:[String:Data]=[:],outputBytes=0
                for output in receipt.outputs {
                    let envelope=try VivoCanonicalJSON.decode(VivoChemistryOutputEnvelope.self,from:try await store.data(for:output.artifact,verify:true))
                    let descriptor=try await store.descriptor(for:output.artifact)
                    guard descriptor.kind=="chemistry-output",envelope.schema=="numivivo.org/chemistry-output/v1",
                          envelope.taskFingerprint==id,envelope.outputName==output.name,envelope.kind==output.kind,
                          envelope.payload.count<=task.resources.maximumOutputBytes-outputBytes else {
                        throw VivoChemistryError.invalid("cached output provenance or size mismatch")
                    }
                    outputBytes += envelope.payload.count;payloads[output.name]=envelope.payload
                }
                try operation.validateOutputs(task.configuration,inputs,payloads,task.resources.budget)
                return .init(taskFingerprint:id,receiptFingerprint:cached.artifact.fingerprint,outputs:receipt.outputs,reused:true)
            }
            try Task.checkCancellation()
            let payloads=try await operation.execute(task.configuration,inputs,task.resources.budget)
            guard Set(payloads.keys)==Set(operation.outputs.map(\.name)) else { throw VivoChemistryError.invalid("operation returned unexpected/missing outputs") }
            var outputBytes=0
            for data in payloads.values {
                guard data.count<=task.resources.maximumOutputBytes-outputBytes else { throw VivoChemistryError.resourceLimit("workflow output byte budget") }
                outputBytes += data.count
            }
            try operation.validateOutputs(task.configuration,inputs,payloads,task.resources.budget)
            var receipts:[VivoChemistryOutputReceipt]=[]
            for output in operation.outputs {
                let envelope=VivoChemistryOutputEnvelope(schema:"numivivo.org/chemistry-output/v1",taskFingerprint:id,
                                                         outputName:output.name,kind:output.kind,payload:payloads[output.name]!)
                let stored=try await store.put(data:VivoCanonicalJSON.encode(envelope),kind:"chemistry-output",mediaType:"application/vnd.numivivo.chemistry-output+json")
                guard try await store.verify(stored.fingerprint) else { throw VivoChemistryError.invalid("fresh output integrity failure") }
                receipts.append(.init(name:output.name,kind:output.kind,artifact:stored.fingerprint))
            }
            let receipt=VivoChemistryTaskReceipt(schema:"numivivo.org/chemistry-task-receipt/v1",taskFingerprint:id,outputs:receipts)
            let stored=try await store.put(data:VivoCanonicalJSON.encode(receipt),kind:"chemistry-task-receipt",mediaType:"application/vnd.numivivo.chemistry-task-receipt+json")
            _ = try await store.setReference(referenceName,to:stored)
            return .init(taskFingerprint:id,receiptFingerprint:stored.fingerprint,outputs:receipts,reused:false)
        }
        inFlight[id]=job
        do { let result=try await job.value;inFlight[id]=nil;return result }
        catch { inFlight[id]=nil;throw error }
    }
    public func runDAG(_ nodes: [VivoChemistryDAGNode]) async throws -> [String:VivoChemistryTaskResult] {
        guard Set(nodes.map(\.identifier)).count==nodes.count,nodes.allSatisfy({!$0.identifier.isEmpty}) else { throw VivoChemistryError.invalid("duplicate/empty DAG node") }
        let identifiers=Set(nodes.map(\.identifier));var remaining=nodes.sorted{$0.identifier<$1.identifier}
        var results:[String:VivoChemistryTaskResult]=[:]
        var dependencies:[String:Set<String>]=[:]
        for node in nodes {
            var deps=Set<String>(),inputNames=Set<String>()
            for input in node.inputs {
                let name:String
                switch input {
                case .artifact(let value,_,_): name=value
                case .output(let value,let upstream,let outputName,let kind):
                    name=value
                    guard let producer=nodes.first(where:{$0.identifier==upstream}),
                          producer.operation.outputs.contains(where:{$0.name==outputName && $0.kind==kind}) else {
                        throw VivoChemistryError.invalid("DAG edge names an unknown or incompatible output slot")
                    }
                }
                guard !name.isEmpty,inputNames.insert(name).inserted else { throw VivoChemistryError.invalid("duplicate/empty DAG input slot") }
            }
            for input in node.inputs { if case .output(_,let upstream,_,_)=input {
                guard identifiers.contains(upstream),upstream != node.identifier else { throw VivoChemistryError.invalid("unknown/self DAG dependency") };deps.insert(upstream)
            } };dependencies[node.identifier]=deps
        }
        var reached=Set<String>()
        while reached.count<nodes.count {
            let ready=nodes.filter{!reached.contains($0.identifier) && dependencies[$0.identifier]!.isSubset(of:reached)}
            guard !ready.isEmpty else { throw VivoChemistryError.invalid("cyclic chemistry DAG") };reached.formUnion(ready.map(\.identifier))
        }
        while !remaining.isEmpty {
            let ready=remaining.filter{dependencies[$0.identifier]!.isSubset(of:Set(results.keys))}
            for node in ready {
                let inputs=try node.inputs.map { input -> VivoChemistryTaskInput in
                    switch input {
                    case .artifact(let name,let fingerprint,let kind): return .init(name:name,artifact:fingerprint,kind:kind)
                    case .output(let name,let upstream,let outputName,let kind):
                        guard let output=results[upstream]?.outputs.first(where:{$0.name==outputName}),output.kind==kind else { throw VivoChemistryError.invalid("missing/incompatible DAG output slot") }
                        return .init(name:name,artifact:output.artifact,kind:kind)
                    }
                }
                let task=VivoChemistryTask(operation:node.operation.identifier,version:node.operation.version,
                    implementationFingerprint:node.operation.implementationFingerprint,inputs:inputs,configuration:node.configuration,
                    outputs:node.operation.outputs,resources:node.resources)
                results[node.identifier]=try await run(task,using:node.operation)
            }
            let completed=Set(ready.map(\.identifier));remaining.removeAll{completed.contains($0.identifier)}
        }
        return results
    }
}
