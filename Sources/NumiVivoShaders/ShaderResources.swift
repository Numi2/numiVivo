import Foundation

public enum NumiVivoShaderResources: Sendable {
    public static var bundle: Bundle { .module }

    public static func sourceURL(named name: String = "NumiVivoKernels") -> URL? {
        bundle.url(forResource: name, withExtension: "metal")
    }

    public static func libraryURL(named name: String = "NumiVivo") -> URL? {
        bundle.url(forResource: name, withExtension: "metallib")
    }
}
