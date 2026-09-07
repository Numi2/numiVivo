import Foundation

/// Public kinetic-document boundary over the EXISTING pinned-root filesystem
/// implementation. No separate artifact store, parser or reference namespace.
public enum VivoKineticsDocumentIO {
    /// A no-clobber destination whose parent directory is captured before any
    /// later asynchronous work. Renaming or replacing its path cannot redirect
    /// publication away from that already opened directory.
    public struct PreparedOutput: Sendable {
        private let files: VivoRootedFileStore
        private let leaf: String
        private let maximumBytes: Int

        fileprivate init(files: VivoRootedFileStore, leaf: String, maximumBytes: Int) {
            self.files = files; self.leaf = leaf; self.maximumBytes = maximumBytes
        }

        public func write(_ data: Data) throws {
            guard data.count <= maximumBytes else { throw VivoKineticsError.capacity("document output") }
            guard try files.writeFile(data, relative: leaf, immutable: true) else {
                throw VivoKineticsError.invalid("output already exists; explicitly request overwrite")
            }
        }
    }

    /// Validate all destinations before preparing them: preparation can create
    /// missing parent directories, but never publishes the document itself.
    /// The caller owns any input/store alias policy before capturing this handle.
    public static func prepareNoClobberOutput(to url: URL, maximumBytes: Int = 536_870_912) throws -> PreparedOutput {
        guard url.isFileURL, (1...536_870_912).contains(maximumBytes) else {
            throw VivoKineticsError.invalid("document URL or byte limit")
        }
        let leaf = url.lastPathComponent
        guard !leaf.isEmpty, leaf != ".", leaf != "..", leaf.utf8.count <= 240,
              !leaf.contains("/"), !leaf.contains("\\"), !leaf.contains("\0") else {
            throw VivoKineticsError.invalid("document output filename")
        }
        let files = try VivoRootedFileStore(rootURL: url.deletingLastPathComponent(), createIfNeeded: true)
        return PreparedOutput(files: files, leaf: leaf, maximumBytes: maximumBytes)
    }

    public static func read<T: Decodable>(_ type: T.Type, from url: URL,
                                          maximumBytes: Int = 64 * 1024 * 1024) throws -> T {
        guard url.isFileURL, (1...536_870_912).contains(maximumBytes) else {
            throw VivoKineticsError.invalid("document URL or byte limit")
        }
        let files = try VivoRootedFileStore(rootURL: url.deletingLastPathComponent(), createIfNeeded: false)
        let data = try files.readFile(url.lastPathComponent, maximumBytes: maximumBytes)
        return try VivoCanonicalJSON.decode(type, from: data)
    }
    public static func write(_ data: Data, to url: URL, overwrite: Bool = false) throws {
        guard url.isFileURL, data.count <= 536_870_912 else { throw VivoKineticsError.capacity("document output") }
        let files = try VivoRootedFileStore(rootURL: url.deletingLastPathComponent(), createIfNeeded: true)
        guard try files.writeFile(data, relative: url.lastPathComponent, immutable: !overwrite) else {
            throw VivoKineticsError.invalid("output already exists; explicitly request overwrite")
        }
    }
}
