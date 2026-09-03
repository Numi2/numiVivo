// swift-tools-version: 6.0

import PackageDescription

let package = Package(
    name: "NumiVivo",
    platforms: [
        .macOS(.v15),
        .iOS(.v18)
    ],
    products: [
        .library(name: "NumiVivoKit", targets: ["NumiVivoKit"]),
        .executable(name: "numivivo", targets: ["NumiVivoCLI"])
    ],
    targets: [
        .target(
            name: "NumiVivoCore",
            path: "Sources/NumiVivoCore",
            publicHeadersPath: "include",
            cxxSettings: [
                .headerSearchPath("include"),
                .define("NVIVO_BUILDING_CORE")
            ]
        ),
        .target(
            name: "NumiVivoShaders",
            path: "Sources/NumiVivoShaders",
            resources: [
                .process("Resources")
            ]
        ),
        .target(
            name: "NumiVivoKit",
            dependencies: [
                "NumiVivoCore",
                "NumiVivoShaders"
            ],
            path: "Sources/NumiVivoKit",
            linkerSettings: [
                .linkedFramework("Metal"),
                .linkedFramework("Accelerate")
            ]
        ),
        .executableTarget(
            name: "NumiVivoCLI",
            dependencies: ["NumiVivoKit"],
            path: "Sources/NumiVivoCLI"
        )
    ],
    cxxLanguageStandard: .cxx23,
    swiftLanguageModes: [.v6]
)
