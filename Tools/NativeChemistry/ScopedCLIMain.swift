import Foundation

/// Test executable entrypoint only. Dispatch, parsing, storage and numerics use
/// the real production electronic and shared-path command classes.
@main struct ScopedChemistryCLI {
    static func main() async {
        let arguments=Array(CommandLine.arguments.dropFirst())
        if VivoECCPathCLICommands.handles(arguments.first) {
            exit(await VivoECCPathCLICommands().run(arguments:arguments))
        }
        exit(await VivoElectronicCLICommands().run(arguments:arguments))
    }
}
