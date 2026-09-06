import Foundation

extension VivoChemistryError {
    /// Local source-compatibility helper for mapped free-energy preparation.
    /// Missing atom-to-particle mappings are invalid prepared inputs, not a new
    /// public error category.
    static func unresolved(_ message: String) -> VivoChemistryError { .invalid(message) }
}
