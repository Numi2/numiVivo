import Foundation
import NumiVivoKit

/// Link this unchanged driver to the retained pre-change owner and new owner.
@main struct CPURegression {
    static func main() throws {
        let args = CommandLine.arguments
        let implementation = try VivoFingerprint(bytes: Array(repeating: 0, count: 32))
        let result = try VivoH5ADCountStore.normalize(URL(fileURLWithPath: args[1]), target: 10_000,
            implementation: implementation, to: URL(fileURLWithPath: args[2]))
        FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(result))
    }
}
