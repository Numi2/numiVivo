import Foundation
import NumiVivoKit

/// Shared path preflight before capturing existing pinned-root output handles.
enum VivoWorkflowCLIDocumentPaths {
    /// Resolve existing ancestors without treating an unresolved symbolic link
    /// as a safe absent leaf. OutputPlan then captures the existing pinned-root,
    /// immutable-link publication boundary before asynchronous work begins.
    static func canonicalURL(_ url: URL) throws -> URL {
        let manager = FileManager.default
        var ancestor = url.absoluteURL, suffix: [String] = []
        guard ancestor.isFileURL, ancestor.path.utf8.count <= 8192, !ancestor.path.contains("\0") else {
            throw VivoChemistryError.invalid("workflow document needs a bounded local path")
        }
        while !manager.fileExists(atPath: ancestor.path) {
            if let attributes = try? manager.attributesOfItem(atPath: ancestor.path),
               attributes[.type] as? FileAttributeType == .typeSymbolicLink {
                throw VivoChemistryError.invalid("workflow document path contains a dangling symbolic link")
            }
            let parent = ancestor.deletingLastPathComponent()
            guard parent.path != ancestor.path, !ancestor.lastPathComponent.isEmpty, suffix.count < 4096 else {
                throw VivoChemistryError.invalid("workflow document path has no existing ancestor")
            }
            suffix.append(ancestor.lastPathComponent); ancestor = parent
        }
        var result = ancestor.resolvingSymlinksInPath().standardizedFileURL
        if !suffix.isEmpty {
            guard try result.resourceValues(forKeys: [.isDirectoryKey]).isDirectory == true else {
                throw VivoChemistryError.invalid("workflow document path parent is not a directory")
            }
        }
        for component in suffix.reversed() { result.appendPathComponent(component) }
        return result.standardizedFileURL
    }

    static func aliases(_ left: URL, _ right: URL) throws -> Bool {
        let a = try canonicalURL(left), b = try canonicalURL(right)
        if a == b { return true }
        let manager = FileManager.default
        guard manager.fileExists(atPath: a.path), manager.fileExists(atPath: b.path) else { return false }
        let x = try manager.attributesOfItem(atPath: a.path), y = try manager.attributesOfItem(atPath: b.path)
        guard let xi = x[.systemFileNumber] as? NSNumber, let yi = y[.systemFileNumber] as? NSNumber,
              let xd = x[.systemNumber] as? NSNumber, let yd = y[.systemNumber] as? NSNumber else {
            throw VivoChemistryError.invalid("cannot establish workflow output file identity")
        }
        return xi == yi && xd == yd
    }

}
