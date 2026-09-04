// swift-tools-version: 6.0

import PackageDescription

let package = Package(
    name: "DigitalTissueHomeostasis",
    platforms: [.macOS(.v15)],
    dependencies: [
        .package(path: "../..")
    ],
    targets: [
        .executableTarget(
            name: "DigitalTissueHomeostasis",
            dependencies: [
                .product(name: "NumiVivoKit", package: "NumiVivo")
            ]
        )
    ],
    swiftLanguageModes: [.v6]
)
