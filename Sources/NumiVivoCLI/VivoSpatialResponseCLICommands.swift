import Foundation
import NumiVivoLearning
struct VivoSpatialResponseCLICommands {
    func run(arguments: [String]) -> Int32 {
        do { try VivoSpatialResponseLearning.run(Array(arguments.dropFirst())); return 0 }
        catch { FileHandle.standardError.write(Data("spatial-response: \(error)\n".utf8)); return 1 }
    }
}
