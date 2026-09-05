import Foundation

/// Public kinetic-document boundary over the EXISTING pinned-root filesystem
/// implementation. No separate artifact store, parser or reference namespace.
public enum VivoKineticsDocumentIO {
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
