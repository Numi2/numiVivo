import Foundation

let arguments = Array(CommandLine.arguments.dropFirst())
let status: Int32
if VivoReactionCLICommands.handles(arguments.first) {
    status = await VivoReactionCLICommands().run(arguments: arguments)
} else if VivoECCPathCLICommands.handles(arguments.first) {
    status = await VivoECCPathCLICommands().run(arguments: arguments)
} else if VivoScreenedInferenceCLICommands.handles(arguments.first) {
    status = await VivoScreenedInferenceCLICommands().run(arguments: arguments)
} else if VivoResearchCLICommands.handles(arguments.first) {
    status = await VivoResearchCLICommands().run(arguments: arguments)
} else if VivoElectronicCLICommands.handles(arguments.first) {
    status = await VivoElectronicCLICommands().run(arguments: arguments)
} else if VivoPosteriorCLICommands.handles(arguments.first) {
    status = await VivoPosteriorCLICommands().run(arguments: arguments)
} else if VivoTargetEngagementCLICommands.handles(arguments.first) {
    status = await VivoTargetEngagementCLICommands().run(arguments: arguments)
} else if VivoMDProtocolCLICommands.handles(arguments.first) {
    status = await VivoMDProtocolCLICommands().run(arguments: arguments)
} else if VivoHybridCLICommands.handles(arguments.first) {
    status = await VivoHybridCLICommands().run(arguments: arguments)
} else if VivoMDCLICommands.handles(arguments.first) {
    status = await VivoMDCLICommands().run(arguments: arguments)
} else if VivoStructureCLICommands.handles(arguments.first) {
    status = VivoStructureCLICommands().run(arguments: arguments)
} else if VivoStructurePrepCLICommands.handles(arguments.first) {
    status = VivoStructurePrepCLICommands().run(arguments: arguments)
} else if VivoForceFieldCLICommands.handles(arguments.first) {
    status = VivoForceFieldCLICommands().run(arguments: arguments)
} else if VivoAmberExecutableCLICommands.handles(arguments.first) {
    status = VivoAmberExecutableCLICommands().run(arguments: arguments)
} else if VivoAmberCLICommands.handles(arguments.first) {
    status = VivoAmberCLICommands().run(arguments: arguments)
} else {
    status = VivoCLICommandRouter().run(arguments: arguments)
    if arguments.isEmpty || ["help", "--help", "-h"].contains(arguments.first ?? "") {
        FileHandle.standardOutput.write(Data("\nReaction qualification: reaction-template, reaction-run, reaction-help.\n".utf8))
        FileHandle.standardOutput.write(Data("\nShared-orbital ECC paths: chemistry-path-template, chemistry-path-run, chemistry-path-help.\n".utf8))
        FileHandle.standardOutput.write(Data("\nMetal-screened inference: engagement-screen-check, engagement-fit-screened, screening-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Research workflows: finite-drug-run, engagement-design, occupancy-benchmark-inspect, occupancy-benchmark-compare, research-help.\n".utf8))
        FileHandle.standardOutput.write(Data("\nElectronic chemistry: chemistry-template, chemistry-run, chemistry-solve, chemistry-help.\n".utf8))
        FileHandle.standardOutput.write(Data("\nKinetic inference: engagement-fit, engagement-predict, engagement-sensitivity, posterior-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Target engagement: engagement-validate, engagement-source, engagement-compile, engagement-run, engagement-batch, engagement-study, engagement-rate, engagement-apply-rate, engagement-help.\n".utf8))
        FileHandle.standardOutput.write(Data("\nMD workflows: md-protocol-template, md-protocol-validate, md-protocol-run, md-protocol-resume, md-protocol-inspect, md-trajectory-inspect, md-protocol-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Hybrid GPU execution: plan-hybrid, run-hybrid, hybrid-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Apple-native MD: md-capabilities, md-run, md-minimize, md-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Molecular structures: structure-import, structure-export, structure-inspect, structure-select, structure-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Structure preparation: structure-resolve-altloc, structure-build-topology, structure-slice, structure-prep-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Force fields: forcefield-validate, forcefield-assign, forcefield-compile, forcefield-help.\n".utf8))
        FileHandle.standardOutput.write(Data("AMBER bridge: amber-import, amber-import-md, amber-help, amber-md-help.\n".utf8))
    }
}
exit(status)
