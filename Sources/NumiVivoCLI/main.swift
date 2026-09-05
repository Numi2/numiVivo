import Foundation

let arguments = Array(CommandLine.arguments.dropFirst())
let status: Int32
if VivoHybridCLICommands.handles(arguments.first) {
    status = await VivoHybridCLICommands().run(arguments: arguments)
} else if VivoStructureCLICommands.handles(arguments.first) {
    status = VivoStructureCLICommands().run(arguments: arguments)
} else if VivoStructurePrepCLICommands.handles(arguments.first) {
    status = VivoStructurePrepCLICommands().run(arguments: arguments)
} else if VivoForceFieldCLICommands.handles(arguments.first) {
    status = VivoForceFieldCLICommands().run(arguments: arguments)
} else {
    status = VivoCLICommandRouter().run(arguments: arguments)
    if arguments.isEmpty || ["help", "--help", "-h"].contains(arguments.first ?? "") {
        FileHandle.standardOutput.write(Data("\nHybrid GPU execution: plan-hybrid, run-hybrid, hybrid-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Molecular structures: structure-import, structure-export, structure-inspect, structure-select, structure-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Structure preparation: structure-resolve-altloc, structure-build-topology, structure-slice, structure-prep-help.\n".utf8))
        FileHandle.standardOutput.write(Data("Force fields: forcefield-validate, forcefield-assign, forcefield-compile, forcefield-help.\n".utf8))
    }
}
exit(status)
