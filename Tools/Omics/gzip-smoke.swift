import Foundation
@main struct GzipSmoke {
    static func main() throws {
        let original = Data("native compressed count fixture\n".utf8)
        // Fixed RFC 1952 member generated from the exact UTF-8 fixture below.
        let gzip = Data(base64Encoded: "H4sIAAAAAAAC/8tLLMksS1VIzs8tKEotLk5NATJL80oU0jIrSkqLUrkA65t8niAAAAA=")!
        let concat = gzip + gzip
        let one = try VivoOmicsSourceDecoder.decode(gzip, maximumExpandedBytes: original.count)
        precondition(one == original)
        let two = try VivoOmicsSourceDecoder.decode(concat, maximumExpandedBytes: 2 * original.count)
        precondition(two == original + original)
        var rejected = 0
        for bad in [gzip.dropLast(), gzip + Data([0]), Data([0x1f, 0x8b])] {
            do { _ = try VivoOmicsSourceDecoder.decode(Data(bad), maximumExpandedBytes: 1000) }
            catch { rejected += 1 }
        }
        do { _ = try VivoOmicsSourceDecoder.decode(gzip, maximumExpandedBytes: original.count - 1) }
        catch { rejected += 1 }
        precondition(rejected == 4)
        let plain = try VivoOmicsSourceDecoder.decode(original, maximumExpandedBytes: original.count)
        precondition(plain == original)
        print("Native gzip smoke passed")
    }
}
