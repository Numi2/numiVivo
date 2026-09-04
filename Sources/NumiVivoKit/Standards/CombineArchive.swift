import Foundation

public struct VivoCombineFile: Sendable, Codable, Equatable {
    public let location: String
    public let formatURI: String
    public let master: Bool
    public let data: Data
    public let fingerprint: VivoFingerprint

    public init(
        location: String,
        formatURI: String,
        master: Bool = false,
        data: Data
    ) throws {
        guard !location.isEmpty,
              !location.hasPrefix("/"),
              !location.split(separator: "/").contains(".."),
              URL(string: formatURI)?.scheme != nil else {
            throw VivoArtifactValidationError.invalid("COMBINE file location or format URI is invalid")
        }
        self.location = location
        self.formatURI = formatURI
        self.master = master
        self.data = data
        self.fingerprint = try VivoCanonicalJSON.fingerprint(data)
    }
}

public struct VivoCombineArchivePlan: Sendable {
    public let files: [VivoCombineFile]
    public let fingerprint: VivoFingerprint

    public init(files: [VivoCombineFile]) throws {
        var locations = Set<String>()
        for file in files {
            guard locations.insert(file.location).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate COMBINE location \(file.location)")
            }
        }
        struct DigestRecord: Codable {
            let location: String
            let formatURI: String
            let master: Bool
            let fingerprint: VivoFingerprint
        }
        self.files = files.sorted { $0.location < $1.location }
        self.fingerprint = try VivoCanonicalJSON.fingerprint(
            VivoCanonicalJSON.encode(self.files.map {
                DigestRecord(
                    location: $0.location,
                    formatURI: $0.formatURI,
                    master: $0.master,
                    fingerprint: $0.fingerprint
                )
            })
        )
    }

    public func writeDirectory(to root: URL) throws {
        let manager = FileManager.default
        let parent = root.deletingLastPathComponent()
        try manager.createDirectory(at: parent, withIntermediateDirectories: true)
        let staging = parent.appendingPathComponent(".\(root.lastPathComponent).\(UUID().uuidString).staging", isDirectory: true)
        try manager.createDirectory(at: staging, withIntermediateDirectories: true)
        do {
            for file in files {
                let destination = staging.appendingPathComponent(file.location)
                try manager.createDirectory(
                    at: destination.deletingLastPathComponent(),
                    withIntermediateDirectories: true
                )
                try file.data.write(to: destination, options: [.atomic])
            }
            if manager.fileExists(atPath: root.path) {
                let backup = parent.appendingPathComponent(".\(root.lastPathComponent).\(UUID().uuidString).backup")
                try manager.moveItem(at: root, to: backup)
                do {
                    try manager.moveItem(at: staging, to: root)
                    try? manager.removeItem(at: backup)
                } catch {
                    try? manager.moveItem(at: backup, to: root)
                    throw error
                }
            } else {
                try manager.moveItem(at: staging, to: root)
            }
        } catch {
            try? manager.removeItem(at: staging)
            throw error
        }
    }
}

public struct VivoCombineArchiveBuilder: Sendable {
    public init() {}

    public func build(
        programPack: VivoProgramPack,
        sbml: VivoSBML3Export,
        sbol: VivoSBOL3Export? = nil,
        sedml: VivoSEDMLExport? = nil,
        additionalFiles: [VivoCombineFile] = []
    ) throws -> VivoCombineArchivePlan {
        var scientificFiles: [VivoCombineFile] = [
            try .init(
                location: "program.nvivo",
                formatURI: "https://numivivo.org/formats/program-pack/v1",
                data: programPack.data
            ),
            try .init(
                location: "model.xml",
                formatURI: "http://identifiers.org/combine.specifications/sbml.level-3.version-2",
                master: sedml == nil,
                data: sbml.xml
            )
        ]
        if let sbol {
            scientificFiles.append(try .init(
                location: "design.ttl",
                formatURI: "http://identifiers.org/combine.specifications/sbol",
                data: sbol.turtle
            ))
        }
        if let sedml {
            scientificFiles.append(try .init(
                location: "simulation.sedml",
                formatURI: "http://identifiers.org/combine.specifications/sed-ml.level-1.version-4",
                master: true,
                data: sedml.xml
            ))
        }
        scientificFiles.append(contentsOf: additionalFiles)

        let metadata = try metadataRDF(
            programPack: programPack,
            scientificFiles: scientificFiles
        )
        scientificFiles.append(try .init(
            location: "metadata.rdf",
            formatURI: "http://identifiers.org/combine.specifications/omex-metadata",
            data: metadata
        ))
        let manifest = try manifestXML(scientificFiles)
        scientificFiles.append(try .init(
            location: "manifest.xml",
            formatURI: "http://identifiers.org/combine.specifications/omex-manifest",
            data: manifest
        ))
        return try VivoCombineArchivePlan(files: scientificFiles)
    }

    private func manifestXML(_ files: [VivoCombineFile]) throws -> Data {
        var lines = [
            #"<?xml version="1.0" encoding="UTF-8"?>"#,
            #"<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">"#,
            #"  <content location="." format="http://identifiers.org/combine.specifications/omex"/>"#
        ]
        for file in files.sorted(by: { $0.location < $1.location }) {
            let master = file.master ? #" master="true""# : ""
            lines.append(#"  <content location="\#(VivoXMLCodec.escapedAttribute(file.location))" format="\#(VivoXMLCodec.escapedAttribute(file.formatURI))"\#(master)/>"#)
        }
        lines.append("</omexManifest>")
        guard let data = (lines.joined(separator: "\n") + "\n").data(using: .utf8) else {
            throw VivoArtifactValidationError.invalid("COMBINE manifest did not encode as UTF-8")
        }
        return data
    }

    private func metadataRDF(
        programPack: VivoProgramPack,
        scientificFiles: [VivoCombineFile]
    ) throws -> Data {
        var lines = [
            #"<?xml version="1.0" encoding="UTF-8"?>"#,
            #"<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:prov="http://www.w3.org/ns/prov#" xmlns:numivivo="https://numivivo.org/ns#">"#,
            #"  <rdf:Description rdf:about=".">"#,
            #"    <dcterms:title>NumiVivo computational molecular biotechnology archive</dcterms:title>"#,
            #"    <numivivo:programPackFingerprint>\#(programPack.header.contentFingerprint.hex)</numivivo:programPackFingerprint>"#,
            #"    <numivivo:sourceFingerprint>\#(programPack.header.sourceFingerprint.hex)</numivivo:sourceFingerprint>"#,
            #"    <numivivo:fidelity>\#(programPack.header.fidelity.label)</numivivo:fidelity>"#,
            #"  </rdf:Description>"#
        ]
        for file in scientificFiles.sorted(by: { $0.location < $1.location }) {
            lines.append(#"  <rdf:Description rdf:about="\#(VivoXMLCodec.escapedAttribute(file.location))">"#)
            lines.append(#"    <dcterms:format>\#(VivoXMLCodec.escapedText(file.formatURI))</dcterms:format>"#)
            lines.append(#"    <numivivo:sha256>\#(file.fingerprint.hex)</numivivo:sha256>"#)
            lines.append(#"    <prov:wasDerivedFrom rdf:resource="urn:sha256:\#(programPack.header.contentFingerprint.hex)"/>"#)
            lines.append("  </rdf:Description>")
        }
        lines.append("</rdf:RDF>")
        guard let data = (lines.joined(separator: "\n") + "\n").data(using: .utf8) else {
            throw VivoArtifactValidationError.invalid("COMBINE metadata did not encode as UTF-8")
        }
        return data
    }
}
