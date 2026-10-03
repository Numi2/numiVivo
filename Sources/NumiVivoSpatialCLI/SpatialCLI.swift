import Foundation
import CryptoKit
import Metal
#if canImport(NumiVivoKit)
import NumiVivoKit
#endif

struct SpatialFailure: Error, CustomStringConvertible { let description: String }
func need(_ ok: Bool, _ message: String) throws { if !ok { throw SpatialFailure(description: message) } }
struct Geometry: Codable {
    let format: String, status: String, owner: String, sourceSHA256: String
    let sourcePhysicsFingerprint: String, deviceProgramFingerprint: String, evidenceStatus: String
    let nodesMetres: [[Double]], referenceNodesMetres: [[Double]], tetrahedra: [[Int]]
    let acceptedSteps: Int
    let timeSeconds: Double
}
struct Cell: Codable { let id: String; let tetrahedron: Int; let volumeCubicMetres: Double; let membranePermeabilityMetresPerSecond: Double }
struct Pulse: Codable { let timeSeconds: Double; let extracellularIndices: [Int]; let concentrationIncrementMolPerM3: Double }
struct Plan: Codable {
    let format: String, specimenID: String
    let diffusionSquareMetresPerSecond: Double, timeStepSeconds: Double
    let extracellularInitialMolPerM3: [Double], intracellularInitialMolPerM3: [Double]
    let cells: [Cell], pulses: [Pulse], sampleTimesSeconds: [Double]
}
struct Compartment: Codable {
    let id: String, kind: String
    let positionMetres: [Double], volumeCubicMetres: Double
    let tetrahedron: Int
}
struct Edge: Codable { let source: Int, target: Int; let clearanceCubicMetresPerSecond: Double; let kind: String }
struct Sample: Codable { let timeSeconds: Double; let concentrationsMolPerM3: [Float]; let amountMol: Double; let cumulativeInjectedMol: Double }
struct Arm: Codable { let name: String; let samples: [Sample]; let acceptedSteps: Int; let maximumRelativeMassError: Double }
struct Output: Codable {
    let format: String, specimenID: String, geometrySHA256: String, planSHA256: String, modelFingerprint: String
    let device: String, backend: String
    let compartments: [Compartment], edges: [Edge], arms: [Arm]
    let evidence: [String: String]
}
func hash(_ bytes: Data) -> String { SHA256.hash(data: bytes).map{String(format:"%02x",$0)}.joined() }
func sub(_ a: [Double], _ b: [Double]) -> [Double] { zip(a,b).map(-) }
func dot(_ a: [Double], _ b: [Double]) -> Double { zip(a,b).map(*).reduce(0,+) }
func cross(_ a: [Double], _ b: [Double]) -> [Double] { [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]] }
func norm(_ a: [Double]) -> Double { sqrt(dot(a,a)) }
func average(_ a: [[Double]]) -> [Double] { (0..<3).map { j in a.map{$0[j]}.reduce(0,+)/Double(a.count) } }

@main struct SpatialCLI {
    static func main() async {
        do { try await run() }
        catch { FileHandle.standardError.write(Data("spatial-tissue: \(error)\n".utf8)); exit(1) }
    }
    static func run() async throws {
        let args=Array(CommandLine.arguments.dropFirst())
        try need(args.count==3,"usage: numivivo-spatial plan.json accepted-geometry.json new-result.json")
        let output=URL(fileURLWithPath:args[2]);try need(!FileManager.default.fileExists(atPath:output.path),"output already exists")
        let planBytes=try Data(contentsOf: URL(fileURLWithPath:args[0])),geometryBytes=try Data(contentsOf: URL(fileURLWithPath:args[1]))
        try need(planBytes.count<=4_194_304 && geometryBytes.count<=4_194_304,"input too large")
        let raw=try JSONSerialization.jsonObject(with:planBytes)
        func keys(_ value: Any, _ expected: Set<String>) throws {
            guard let object=value as? [String:Any] else { throw SpatialFailure(description:"expected plan object") }
            try need(Set(object.keys)==expected,"missing or unknown spatial plan field")
        }
        try keys(raw,["format","specimenID","diffusionSquareMetresPerSecond","timeStepSeconds","extracellularInitialMolPerM3","intracellularInitialMolPerM3","cells","pulses","sampleTimesSeconds"])
        let object=raw as! [String:Any]
        for value in object["cells"] as? [Any] ?? [] { try keys(value,["id","tetrahedron","volumeCubicMetres","membranePermeabilityMetresPerSecond"]) }
        for value in object["pulses"] as? [Any] ?? [] { try keys(value,["timeSeconds","extracellularIndices","concentrationIncrementMolPerM3"]) }
        let p=try JSONDecoder().decode(Plan.self,from:planBytes),g=try JSONDecoder().decode(Geometry.self,from:geometryBytes)
        try need(p.format=="numivivo-spatial-assay/v1" && !p.specimenID.isEmpty,"unsupported plan")
        try need(g.format=="numilab-wet-lab-geometry/v1" && g.status=="accepted" && g.owner=="NumiLab Matter" && g.acceptedSteps>0,"requires accepted Matter geometry")
        try need(!g.tetrahedra.isEmpty && g.tetrahedra.count<=8192 && g.nodesMetres.count<=4096,"geometry size")
        try need(g.nodesMetres.allSatisfy{$0.count==3 && $0.allSatisfy{$0.isFinite}},"geometry coordinates")
        try need(p.diffusionSquareMetresPerSecond.isFinite && p.diffusionSquareMetresPerSecond>=0 && p.diffusionSquareMetresPerSecond<=1e-6,"diffusivity")
        try need(p.timeStepSeconds.isFinite && p.timeStepSeconds>=1e-5 && p.timeStepSeconds<=10,"time step")
        try need(!p.sampleTimesSeconds.isEmpty && p.sampleTimesSeconds.count<=128 && p.sampleTimesSeconds[0]==0 && p.sampleTimesSeconds.last!>0 && p.sampleTimesSeconds.last!<=3600,"sample times")
        try need(zip(p.sampleTimesSeconds,p.sampleTimesSeconds.dropFirst()).allSatisfy{$0<$1},"sample times must increase")
        let end=p.sampleTimesSeconds.last!
        try need(p.sampleTimesSeconds.allSatisfy{$0.isFinite} && end/p.timeStepSeconds<=100_000,"step budget")
        try need(p.cells.count>0 && p.cells.count<=8192 && Set(p.cells.map(\.id)).count==p.cells.count,"cell identities")
        try need(Set(g.tetrahedra.map{$0.sorted().map(String.init).joined(separator:",")}).count==g.tetrahedra.count,"duplicate tetrahedron")
        var volumes:[Double]=[],centers:[[Double]]=[],faces:[String:[(Int,[Int])]]=[:]
        for (index,t) in g.tetrahedra.enumerated() {
            try need(t.count==4 && Set(t).count==4 && t.allSatisfy{g.nodesMetres.indices.contains($0)},"tetrahedron topology")
            let n=t.map{g.nodesMetres[$0]},v=dot(sub(n[1],n[0]),cross(sub(n[2],n[0]),sub(n[3],n[0])))/6
            try need(v.isFinite && v>1e-24,"inverted/degenerate accepted tetrahedron")
            volumes.append(v);centers.append(average(n))
            for omit in 0..<4 { let f=t.indices.filter{$0 != omit}.map{t[$0]}.sorted();faces[f.map(String.init).joined(separator:","),default:[]].append((index,f)) }
        }
        let count=volumes.count
        for c in p.cells {
            try need(!c.id.isEmpty && (0..<count).contains(c.tetrahedron) && c.volumeCubicMetres.isFinite && c.volumeCubicMetres>0 && c.membranePermeabilityMetresPerSecond.isFinite && c.membranePermeabilityMetresPerSecond>=0,"cell geometry or membrane law")
            volumes[c.tetrahedron]-=c.volumeCubicMetres
        }
        try need(volumes.allSatisfy{$0>0},"cells exhaust extracellular volume")
        var compartments=(0..<count).map{Compartment(id:"ecs-\($0)",kind:"extracellular",positionMetres:centers[$0],volumeCubicMetres:volumes[$0],tetrahedron:$0)}
        compartments += p.cells.map{Compartment(id:$0.id,kind:"cell",positionMetres:centers[$0.tetrahedron],volumeCubicMetres:$0.volumeCubicMetres,tetrahedron:$0.tetrahedron)}
        try need(Set(compartments.map(\.id)).count==compartments.count,"cell and space identities collide")
        var edges:[Edge]=[]
        for key in faces.keys.sorted() {
            let entries=faces[key]!;try need(entries.count<=2,"non-manifold mesh")
            if entries.count==2 {
                let a=entries[0].0,b=entries[1].0,n=entries[0].1.map{g.nodesMetres[$0]}
                let normal=cross(sub(n[1],n[0]),sub(n[2],n[0])),area=norm(normal)/2
                // Two-point finite volume transmissibility on an explicitly admitted
                // tetrahedral graph. Nonorthogonal correction is not claimed.
                let distance=norm(sub(centers[a],centers[b]));try need(distance>0,"coincident cell centers")
                edges.append(.init(source:a,target:b,clearanceCubicMetresPerSecond:p.diffusionSquareMetresPerSecond*area/distance,kind:"extracellular-diffusion"))
            }
        }
        for (i,c) in p.cells.enumerated() {
            let radius=pow(3*c.volumeCubicMetres/(4*Double.pi),1.0/3)
            edges.append(.init(source:c.tetrahedron,target:count+i,clearanceCubicMetresPerSecond:c.membranePermeabilityMetresPerSecond*4*Double.pi*radius*radius,kind:"membrane-exchange"))
        }
        try need(p.extracellularInitialMolPerM3.count==count && p.intracellularInitialMolPerM3.count==p.cells.count,"initial field dimensions")
        let initial=p.extracellularInitialMolPerM3+p.intracellularInitialMolPerM3
        try need(initial.allSatisfy{$0.isFinite && $0>=0 && $0<=1000},"initial concentration bounds")
        try need(p.pulses.count<=64,"pulse budget")
        for pulse in p.pulses {
            try need(pulse.timeSeconds.isFinite && pulse.timeSeconds>=0 && pulse.timeSeconds<=end && !pulse.extracellularIndices.isEmpty && Set(pulse.extracellularIndices).count==pulse.extracellularIndices.count && pulse.extracellularIndices.allSatisfy{(0..<count).contains($0)} && pulse.concentrationIncrementMolPerM3.isFinite && pulse.concentrationIncrementMolPerM3>=0 && pulse.concentrationIncrementMolPerM3<=1000,"intervention outside support")
        }
        let model=try VivoPhysiologicalPartitionModel(name:p.specimenID,compartments:compartments.map{.init(id:$0.id,volumeCubicMetres:Float($0.volumeCubicMetres))},analytes:[.init(id:"abstract-tracer",unit:"mol/m3",maximumConcentration:1e6)],edges:edges.enumerated().map{i,e in .init(id:"edge-\(i)",analyte:"abstract-tracer",sourceCompartment:compartments[e.source].id,targetCompartment:compartments[e.target].id,partitionCoefficient:1,clearanceCubicMetresPerSecond:Float(e.clearanceCubicMetresPerSecond),evidenceClass:"hypothetical")})
        guard let device=MTLCreateSystemDefaultDevice() else { throw SpatialFailure(description:"Metal unavailable") }
        var arms:[Arm]=[]
        for treated in [false,true] {
            var values=initial.map(Float.init),samples:[Sample]=[],steps=0,clock=0.0,injected=0.0,maxError=0.0
            let initialAmount=zip(values,model.compartments).map{Double($0)*Double($1.volumeCubicMetres)}.reduce(0,+)
            var runtime=try await VivoPhysiologicalPartitionRuntime.make(model:model,configuration:.init(timeStep:Float(p.timeStepSeconds),minimumTimeStep:1e-8,maximumTimeStep:10),initialConcentrations:values,device:device)
            let events=Set(p.sampleTimesSeconds+p.pulses.map(\.timeSeconds)).sorted()
            for event in events {
                while clock < event-1e-8 {
                    let dt=min(p.timeStepSeconds,event-clock)
                    let certificate=try await runtime.step(deltaTime:Float(dt),permitAdaptiveReduction:true,certifyAmountConservation:true)
                    try need(certificate.disposition != .rejected && certificate.acceptedStep != nil,"transport transaction rejected")
                    clock += Double(certificate.acceptedStep!);steps += 1;try need(steps<=200_000,"accepted step budget exceeded")
                }
                try need(abs(clock-event)<1e-5,"sampling clock mismatch")
                values=try await runtime.snapshot().concentrations
                let here=p.pulses.filter{$0.timeSeconds==event}
                if !here.isEmpty {
                    for pulse in here where treated {
                        for i in pulse.extracellularIndices { let before=values[i];values[i]+=Float(pulse.concentrationIncrementMolPerM3);injected+=Double(values[i]-before)*Double(model.compartments[i].volumeCubicMetres) }
                    }
                    runtime=try await VivoPhysiologicalPartitionRuntime.make(model:model,configuration:.init(timeStep:Float(p.timeStepSeconds),minimumTimeStep:1e-8,maximumTimeStep:10),initialConcentrations:values,device:device)
                }
                let total=zip(values,model.compartments).map{Double($0)*Double($1.volumeCubicMetres)}.reduce(0,+)
                let error=abs(total-initialAmount-injected)/max(1e-30,initialAmount+injected);maxError=max(maxError,error)
                try need(error<=2e-5,"amount balance failed")
                if p.sampleTimesSeconds.contains(event) { samples.append(.init(timeSeconds:event,concentrationsMolPerM3:values,amountMol:total,cumulativeInjectedMol:injected)) }
            }
            arms.append(.init(name:treated ? "intervention":"control",samples:samples,acceptedSteps:steps,maximumRelativeMassError:maxError))
        }
        let result=Output(format:"numivivo-spatial-result/v1",specimenID:p.specimenID,geometrySHA256:hash(geometryBytes),planSHA256:hash(planBytes),modelFingerprint:model.fingerprint,device:device.name,backend:"NumiVivo Metal reversible partition RK2",compartments:compartments,edges:edges,arms:arms,evidence:["mechanics":"accepted-numerical-state","geometryToTransport":"one-way-geometric-transfer","transport":"conservative-discrete-graph","cellResponse":"hypothetical-passive-membrane-exchange","biologicalValidation":"not-established","feedbackToMechanics":"not-implemented","continuumDiffusion":"nonorthogonal-mesh-accuracy-unqualified"])
        let encoder=JSONEncoder();encoder.outputFormatting=[.prettyPrinted,.sortedKeys,.withoutEscapingSlashes]
        try encoder.encode(result).write(to:output,options:.withoutOverwriting)
        print("Spatial assay completed: \(compartments.count) compartments, \(edges.count) edges, \(arms.count) arms")
    }
}
