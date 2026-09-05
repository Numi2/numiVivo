// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "NumiVivoNativeChemistryExample",
    platforms: [.macOS(.v15)],
    products: [.executable(name: "native-chemistry", targets: ["NativeChemistry"])],
    dependencies: [.package(name: "NumiVivo", path: "../..")],
    targets: [.executableTarget(name: "NativeChemistry", dependencies: [.product(name: "NumiVivoKit", package: "NumiVivo")])]
)
