import Foundation

let status = VivoCLICommandRouter().run(
    arguments: Array(CommandLine.arguments.dropFirst())
)
exit(status)
