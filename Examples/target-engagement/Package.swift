// swift-tools-version: 6.0
import PackageDescription
let package = Package(
    name: "NumiVivoTargetEngagementChecks",
    platforms: [.macOS(.v15)],
    products: [.executable(name: "target-engagement-checks", targets: ["TargetEngagementChecks"])],
    dependencies: [.package(path: "../..")],
    targets: [.executableTarget(name: "TargetEngagementChecks", dependencies: [
        .product(name: "NumiVivoKit", package: "numivivo")
    ])],
    swiftLanguageModes: [.v6]
)
