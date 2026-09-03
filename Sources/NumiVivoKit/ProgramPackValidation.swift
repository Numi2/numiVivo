import Foundation
import NumiVivoCore

public extension VivoProgramPack {
    func semanticValidationReport() throws -> Data {
        var report = NVivoByteBuffer(data: nil, size: 0)
        defer { nvivo_buffer_release(&report) }
        let status = data.withUnsafeBytes { rawBuffer in
            nvivo_validate_program_pack_semantics(
                rawBuffer.bindMemory(to: UInt8.self).baseAddress,
                rawBuffer.count,
                &report
            )
        }
        let reportData: Data
        if let pointer = report.data, report.size > 0 {
            reportData = Data(bytes: pointer, count: report.size)
        } else {
            reportData = Data()
        }
        guard status.rawValue == 0 else {
            throw VivoRuntimeError.invalidProgramPack(
                String(decoding: reportData, as: UTF8.self)
            )
        }
        return reportData
    }
}
