import Foundation
import NumiVivoKit

struct VivoMDCLICommands {
    static func handles(_ name:String?)->Bool{["md-capabilities","md-minimize","md-run","md-help"].contains(name ?? "")}

    func run(arguments:[String])async->Int32{
        do{
            guard let command=arguments.first else{throw MDCLIError.usage("missing MD command")}
            if command=="md-help"{FileHandle.standardOutput.write(Data(Self.help.utf8));return 0}
            let parsed=try Arguments(Array(arguments.dropFirst()))
            switch command{
            case "md-capabilities":return try capabilities(parsed)
            case "md-minimize":return try await minimize(parsed)
            case "md-run":return try await runMD(parsed)
            default:throw MDCLIError.usage("unknown MD command")
            }
        }catch let error as MDCLIError{stderr("numivivo md: \(error)\n");return 64}
        catch let error as VivoMDRuntimeError{stderr("numivivo md: \(error.description)\n");return 75}
        catch{stderr("numivivo md: \(error.localizedDescription)\n");return 65}
    }

    private func capabilities(_ a:Arguments)throws->Int32{
        try a.requirePositionals(1,"md-capabilities <system.json> --state <initial-state.json> --config <md-config.json>")
        try a.allow(["state","config","output","force"])
        let system:VivoClassicalSystem=try load(a.positionals[0]),state:VivoClassicalInitialState=try load(a.required("state")),config:VivoMDConfiguration=try load(a.required("config"))
        let report=try VivoMDCapabilityAnalyzer.analyze(system:system,initialState:state,configuration:config)
        try write(VivoCanonicalJSON.encode(report),to:a.options["output"] ?? "-",force:a.force)
        return report.executable ? 0:65
    }

    private func minimize(_ a:Arguments)async throws->Int32{
        try a.requirePositionals(1,"md-minimize <system.json> --config <md-config.json> (--state <initial-state.json> | --restore <checkpoint.json>)")
        try a.allow(["state","restore","config","min-config","output","checkpoint","force"])
        let system:VivoClassicalSystem=try load(a.positionals[0]),config:VivoMDConfiguration=try load(a.required("config"))
        let hasState=a.options["state"] != nil,hasRestore=a.options["restore"] != nil
        guard hasState != hasRestore else{throw MDCLIError.usage("provide exactly one of --state or --restore")}
        let runtime:VivoMDMetalRuntime
        if let path=a.options["restore"]{let checkpoint:VivoMDCheckpoint=try load(path);runtime=try await VivoMDMetalRuntime.restore(system:system,configuration:config,checkpoint:checkpoint)}
        else{let state:VivoClassicalInitialState=try load(a.required("state"));runtime=try await VivoMDMetalRuntime.make(system:system,initialState:state,configuration:config)}
        let minConfig:VivoMDMinimizationConfiguration
        if let path=a.options["min-config"]{minConfig=try load(path)}else{minConfig=.init()}
        let certificate=try await runtime.minimize(minConfig)
        if let destination=a.options["checkpoint"]{let checkpoint=try await runtime.checkpoint();try write(VivoCanonicalJSON.encode(checkpoint),to:destination,force:a.force)}
        try write(VivoCanonicalJSON.encode(certificate),to:a.options["output"] ?? "-",force:a.force)
        return certificate.converged ? 0:75
    }

    private func runMD(_ a:Arguments)async throws->Int32{
        try a.requirePositionals(1,"md-run <system.json> --config <md-config.json> (--state <initial-state.json> | --restore <checkpoint.json>) --steps N")
        try a.allow(["state","restore","config","steps","output","checkpoint","sample-every","observables-every","max-samples","force"])
        let system:VivoClassicalSystem=try load(a.positionals[0]),config:VivoMDConfiguration=try load(a.required("config"))
        let hasState=a.options["state"] != nil,hasRestore=a.options["restore"] != nil
        guard hasState != hasRestore else{throw MDCLIError.usage("provide exactly one of --state or --restore")}
        let requested=try a.uint64("steps",required:true)
        guard requested<=10_000_000_000 else{throw MDCLIError.usage("--steps exceeds per-invocation safety bound")}
        let sampleEvery=try a.optionalUInt64("sample-every"),observablesEvery=try a.optionalUInt64("observables-every")
        if sampleEvery==0 || observablesEvery==0{throw MDCLIError.usage("sampling intervals must be positive")}
        let maxSamples=try a.uint64("max-samples",default:10_000)
        if let sampleEvery, requested/sampleEvery>maxSamples{throw MDCLIError.usage("requested trajectory samples exceed --max-samples")}
        let runtime:VivoMDMetalRuntime
        if let path=a.options["restore"]{let checkpoint:VivoMDCheckpoint=try load(path);runtime=try await VivoMDMetalRuntime.restore(system:system,configuration:config,checkpoint:checkpoint)}
        else{let state:VivoClassicalInitialState=try load(a.required("state"));runtime=try await VivoMDMetalRuntime.make(system:system,initialState:state,configuration:config)}
        let start=try await runtime.checkpoint();var committed:UInt64=0,rejected:VivoMDStepCertificate?,samples:[VivoMDTrajectorySample]=[]
        if requested>0{for ordinal in 1...requested{
            try Task.checkCancellation();let certificate=try await runtime.step();if !certificate.committed{rejected=certificate;break};committed+=1
            let sampleDue=sampleEvery.map{ordinal % $0 == 0} ?? false,obsDue=observablesEvery.map{ordinal % $0 == 0} ?? false
            if sampleDue{let state=try await runtime.snapshot();let obs=obsDue ? try await runtime.observables():nil;samples.append(.init(stepIndex:state.stepIndex,timePS:state.timePS,positionsNM:state.positionsNM,periodicCell:state.periodicCell,observables:obs))}
            else if obsDue{_ = try await runtime.observables()}
        }}
        let final=try await runtime.checkpoint()
        if let destination=a.options["checkpoint"]{try write(VivoCanonicalJSON.encode(final),to:destination,force:a.force)}
        let report=VivoMDRunReport(systemFingerprint:runtime.systemFingerprint,configurationFingerprint:runtime.configurationFingerprint,deviceName:runtime.deviceName,deviceRegistryID:runtime.deviceRegistryID,requestedSteps:requested,committedSteps:committed,startStep:start.acceptedStep,endStep:final.acceptedStep,startTimePS:start.timePS,endTimePS:final.timePS,rejected:rejected,samples:samples,finalCheckpoint:final)
        try write(VivoCanonicalJSON.encode(report),to:a.options["output"] ?? "-",force:a.force)
        return rejected == nil ? 0:75
    }

    private func load<T:Decodable>(_ path:String)throws->T{let data=try Data(contentsOf:URL(fileURLWithPath:path));return try VivoCanonicalJSON.decode(T.self,from:data)}
    private func write(_ data:Data,to path:String,force:Bool)throws{
        if path=="-"{FileHandle.standardOutput.write(data);FileHandle.standardOutput.write(Data("\n".utf8));return}
        let url=URL(fileURLWithPath:path);if FileManager.default.fileExists(atPath:url.path),!force{throw MDCLIError.usage("output exists: \(path); use --force")}
        try FileManager.default.createDirectory(at:url.deletingLastPathComponent(),withIntermediateDirectories:true);try data.write(to:url,options:.atomic)
    }
    private func stderr(_ value:String){FileHandle.standardError.write(Data(value.utf8))}

    static let help="""
    NumiVivo Wave B Apple-native molecular dynamics

      numivivo md-capabilities <system.json> --state <initial-state.json>
          --config <md-config.json> [--output report.json]

      numivivo md-minimize <system.json> --state <initial-state.json>
          --config <md-config.json> [--min-config minimization.json]
          [--checkpoint minimized.json] [--output certificate.json] [--force]

      numivivo md-run <system.json> --state <initial-state.json>
          --config <md-config.json> --steps 10000
          [--sample-every 1000] [--observables-every 1000]
          [--checkpoint checkpoint.json] [--output report.json] [--force]

      numivivo md-run <system.json> --restore <checkpoint.json>
          --config <md-config.json> --steps 10000 [same output options]

    Minimization and MD share the same force authority, including PME, virtual
    sites and force-field exceptions. Accepted dynamics remain GPU-resident; full
    positions are read only for requested samples/checkpoints.
    """+"\n"
}

private enum MDCLIError:Error,CustomStringConvertible{case usage(String);var description:String{switch self{case .usage(let v):return v}}}
private struct Arguments{
    let positionals:[String],options:[String:String],force:Bool
    init(_ raw:[String])throws{var p:[String]=[],o:[String:String]=[:],f=false,i=0;while i<raw.count{let token=raw[i];if token=="--force"{guard !f else{throw MDCLIError.usage("duplicate --force")};f=true;i+=1;continue};if token.hasPrefix("--"){let key=String(token.dropFirst(2));guard o[key]==nil,i+1<raw.count,!raw[i+1].hasPrefix("--")else{throw MDCLIError.usage("missing value or duplicate option \(token)")};o[key]=raw[i+1];i+=2}else{p.append(token);i+=1}};positionals=p;options=o;force=f}
    func requirePositionals(_ n:Int,_ usage:String)throws{guard positionals.count==n else{throw MDCLIError.usage(usage)}}
    func allow(_ values:Set<String>)throws{guard options.keys.allSatisfy(values.contains)else{throw MDCLIError.usage("unknown option")}}
    func required(_ key:String)throws->String{guard let v=options[key],!v.isEmpty else{throw MDCLIError.usage("missing --\(key)")};return v}
    func uint64(_ key:String,required:Bool=false,default fallback:UInt64=0)throws->UInt64{guard let raw=options[key]else{if required{throw MDCLIError.usage("missing --\(key)")};return fallback};guard let v=UInt64(raw)else{throw MDCLIError.usage("--\(key) requires UInt64")};return v}
    func optionalUInt64(_ key:String)throws->UInt64?{guard let raw=options[key]else{return nil};guard let v=UInt64(raw)else{throw MDCLIError.usage("--\(key) requires UInt64")};return v}
}
