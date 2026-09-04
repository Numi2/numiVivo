import Foundation

// Async execution commands install their own dispatcher; artifact-only commands
// remain synchronous and retain their existing exit-code contract.
let status = VivoCLICommandRouter().run(
    arguments: Array(CommandLine.arguments.dropFirst())
)
exit(status)
