import Foundation
#if canImport(NumiVivoKit)
import NumiVivoKit
#endif
@main struct TissueCLI {
    static func main() {
        do {
            let a=CommandLine.arguments
            guard a.count==2 else { throw NSError(domain:"tissue",code:1,userInfo:[NSLocalizedDescriptionKey:"usage: numivivo-tissue specimen.json"]) }
            let data=try Data(contentsOf:URL(fileURLWithPath:a[1]))
            guard data.count<=536_870_912 else { throw NSError(domain:"tissue",code:2) }
            let s=try JSONDecoder().decode(VivoTissueSpecimen.self,from:data);try s.validate()
            let result:[String:Any]=["status":"valid-data-contract","specimen":s.id,"entities":s.entities.count,"measurements":s.measurements.count,"biologicalValidation":"not-conferred-by-schema"]
            FileHandle.standardOutput.write(try JSONSerialization.data(withJSONObject:result,options:[.sortedKeys]));print("")
        } catch { FileHandle.standardError.write(Data("tissue: \(error.localizedDescription)\n".utf8));exit(1) }
    }
}
