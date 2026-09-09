import Foundation
import NumiVivoKit

struct Request: Decodable {
    let design: [[Double]]
    let offsets: [Double]
    let contrast: [Double]
    let counts: [[UInt64]]
    let dispersion: Double
    let profileCount: Int
}
struct Response: Encodable {
    let fits: [VivoOmicsNBFit]
    let profiles: [ProfileAttempt]
}
struct ProfileAttempt: Encodable {
    let fit: VivoOmicsNBDispersionFit?
    let error: String?
}
@main struct Main {
    static func main() throws {
        let args = CommandLine.arguments
        guard args.count == 3 else { throw NSError(domain: "usage: nb-check input.json output.json",code: 1) }
        let input = try JSONDecoder().decode(Request.self,from: Data(contentsOf: URL(fileURLWithPath: args[1])))
        guard (0...input.counts.count).contains(input.profileCount) else { throw NSError(domain: "profile count",code: 1) }
        let fits = try input.counts.map { try VivoOmicsNegativeBinomial.fit(counts: $0, design: input.design,
            offsets: input.offsets, contrast: input.contrast, dispersion: input.dispersion) }
        let profiles = input.counts.prefix(input.profileCount).map { counts in
            do { return ProfileAttempt(fit: try VivoOmicsNegativeBinomial.estimateDispersion(
                counts: counts, design: input.design, offsets: input.offsets, contrast: input.contrast),error: nil)
            } catch { return ProfileAttempt(fit: nil,error: error.localizedDescription) }
        }
        let encoder = JSONEncoder(); encoder.outputFormatting = [.sortedKeys]
        try encoder.encode(Response(fits: fits,profiles: profiles)).write(to: URL(fileURLWithPath: args[2]),options: .atomic)
    }
}
