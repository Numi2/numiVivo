import Foundation

/// Entry point for the actual product single-cell router, linked to real owners.
/// This scoped executable does not claim the rest of the NumiVivo application.
@main struct OmicsCLICheck {
    static func main() async {
        exit(await VivoSingleCellCLICommands().run(arguments: Array(CommandLine.arguments.dropFirst())))
    }
}
