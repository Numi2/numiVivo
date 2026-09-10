import Foundation

struct Request: Codable {
    let id: String
    let design: [[Double]]
    let offsets: [Double]
    let contrast: [Double]
    let counts: [UInt64]
    let dispersions: [Double]
}
struct Attempt: Encodable {
    let dispersion: Double
    let fit: VivoOmicsNBFit?
    let error: String?
}
struct Result: Encodable {
    let id: String
    let profile: VivoOmicsNBDispersionFit?
    let profileError: String?
    let points: [Attempt]
}
@main struct Main {
    static func main() throws {
        let args=CommandLine.arguments
        guard args.count==3 else { throw NSError(domain:"usage: profile-audit input.json output.json",code:1) }
        let requests=try JSONDecoder().decode([Request].self,from:Data(contentsOf:URL(fileURLWithPath:args[1])))
        var results:[Result]=[]
        for r in requests {
            var profile:VivoOmicsNBDispersionFit?,failure:String?
            do { profile=try VivoOmicsNegativeBinomial.estimateDispersion(counts:r.counts,design:r.design,offsets:r.offsets,contrast:r.contrast) }
            catch { failure=error.localizedDescription }
            let points=r.dispersions.map { alpha -> Attempt in
                do { return .init(dispersion:alpha,fit:try VivoOmicsNegativeBinomial.fit(counts:r.counts,design:r.design,offsets:r.offsets,contrast:r.contrast,dispersion:alpha),error:nil) }
                catch { return .init(dispersion:alpha,fit:nil,error:error.localizedDescription) }
            }
            results.append(.init(id:r.id,profile:profile,profileError:failure,points:points))
        }
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys]
        try encoder.encode(results).write(to:URL(fileURLWithPath:args[2]),options:.atomic)
    }
}
