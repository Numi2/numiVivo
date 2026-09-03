import Foundation
import NumiVivoCore

public enum NumiVivo: Sendable {
    public static var version: String {
        String(cString: nvivo_version_string())
    }

    public static var programPackVersion: UInt32 {
        nvivo_program_pack_version()
    }
}
