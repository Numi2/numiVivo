import Foundation

/// Test executable entrypoint only. All dispatch, parsing, storage and numerical
/// work is performed by the unchanged production electronic-command class.
@main struct ScopedChemistryCLI {
    static func main() async {
        exit(await VivoElectronicCLICommands().run(arguments:Array(CommandLine.arguments.dropFirst())))
    }
}
