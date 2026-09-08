import Foundation

public enum VivoSingleCellCampaignIO {
    /// Reuses descriptor-relative O_NOFOLLOW I/O from the artifact store. Sources
    /// cannot escape the selected directory through ../, symlinks or named pipes.
    public static func snapshot(manifestURL: URL, ceiling: VivoOmicsLimits = .init()) throws -> VivoSingleCellInputBundle {
        guard manifestURL.isFileURL else { throw VivoOmicsError.invalid("manifest must be a local file") }
        let url = manifestURL.standardizedFileURL
        let files = try VivoRootedFileStore(rootURL: url.deletingLastPathComponent(), createIfNeeded: false)
        let bytes = try files.readFile(url.lastPathComponent, maximumBytes: VivoSingleCellCampaign.maximumManifestBytes)
        let manifest = try VivoSingleCellCampaign.manifest(from: bytes, ceiling: ceiling)
        let limits = try manifest.admittedLimits(ceiling: ceiling)
        guard bytes.count <= limits.maximumInputBytes else { throw VivoOmicsError.limit("manifest exceeds total input budget") }
        var remaining = limits.maximumInputBytes - bytes.count
        func read(_ path: String) throws -> Data {
            let data = try files.readFile(path, maximumBytes: remaining)
            remaining -= data.count; return data
        }
        var libraries: [VivoSingleCellLibraryBytes] = []
        for source in manifest.libraries {
            try Task.checkCancellation()
            libraries.append(try .init(matrix: read(source.matrix), features: read(source.features), barcodes: read(source.barcodes)))
        }
        return .init(manifest: bytes, libraries: libraries)
    }
}
