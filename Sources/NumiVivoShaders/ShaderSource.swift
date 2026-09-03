import Foundation

public enum NumiVivoShaderSourceError: Error, Sendable {
    case missingResource(String)
    case unreadableResource(String)
}

public extension NumiVivoShaderResources {
    static func completeMetalSource() throws -> String {
        guard let abiURL = bundle.url(forResource: "NumiVivoMetalABI", withExtension: "h") else {
            throw NumiVivoShaderSourceError.missingResource("NumiVivoMetalABI.h")
        }
        guard let kernelURL = bundle.url(forResource: "NumiVivoKernels", withExtension: "metal") else {
            throw NumiVivoShaderSourceError.missingResource("NumiVivoKernels.metal")
        }
        guard let abi = try? String(contentsOf: abiURL, encoding: .utf8) else {
            throw NumiVivoShaderSourceError.unreadableResource("NumiVivoMetalABI.h")
        }
        guard var kernels = try? String(contentsOf: kernelURL, encoding: .utf8) else {
            throw NumiVivoShaderSourceError.unreadableResource("NumiVivoKernels.metal")
        }

        // Runtime source compilation receives the ABI and kernels as one
        // translation unit. The guarded local include remains useful when the
        // .metal file is compiled directly by an Xcode target.
        kernels = kernels.replacingOccurrences(
            of: "#ifndef NUMIVIVO_METAL_ABI_H\n#include \"NumiVivoMetalABI.h\"\n#endif\n",
            with: ""
        )
        return abi + "\n" + kernels
    }
}
