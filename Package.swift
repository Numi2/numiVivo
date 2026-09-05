// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "NumiVivo",
    platforms: [.macOS(.v15), .iOS(.v18)],
    products: [
        .library(name: "NumiVivoKit", targets: ["NumiVivoKit"]),
        .executable(name: "numivivo", targets: ["NumiVivoCLI"])
    ],
    targets: [
        .target(name: "NumiVivoCore", path: "Sources/NumiVivoCore", publicHeadersPath: "include",
                cxxSettings: [.headerSearchPath("include"), .define("NVIVO_BUILDING_CORE")]),
        .target(name: "NumiVivoShaders", path: "Sources/NumiVivoShaders", resources: [
            .copy("Resources/NumiVivoProgramPackRuntime.metal"),
            .copy("Resources/NumiVivoTargetLikelihood.metal"),
            .copy("Resources/NumiVivoHybridExecution.metal"),
            .copy("Resources/NumiVivoPhysiologyKernels.metal"),
            .copy("Resources/NumiVivoMDKernels.metal"),
            .copy("Resources/NumiVivoMDVirtualSites.metal"),
            .copy("Resources/NumiVivoMDNeighborGrid.metal"),
            .copy("Resources/NumiVivoMDPME.metal"),
            .copy("Resources/NumiVivoMDPMECorrections.metal"),
            .copy("Resources/NumiVivoMDMinimization.metal"),
            .copy("Resources/NumiVivoMDBarostat.metal"),
            .copy("Resources/NumiVivoMDExecutionSupport.metal"),
            .copy("Resources/NumiVivoExactSSAKernels.metal"),
            .copy("Resources/NumiVivoMigrationKernels.metal"),
            .copy("Resources/NumiVivoPartitionKernels.metal"),
            .copy("Resources/NumiVivoPopulationKernels.metal"),
            .copy("Resources/NumiVivoTransactionCommitKernels.metal"),
            .copy("Resources/NumiVivoTransactionKernels.metal"),
            .copy("Resources/NumiVivoKernels.metal"),
            .copy("Resources/NumiVivoMetalABI.h"),
            .copy("Resources/README.txt")
        ]),
        .target(name: "NumiVivoKit", dependencies: ["NumiVivoCore", "NumiVivoShaders"],
                path: "Sources/NumiVivoKit", linkerSettings: [.linkedFramework("Metal"), .linkedFramework("Accelerate")]),
        .executableTarget(name: "NumiVivoCLI", dependencies: ["NumiVivoKit"], path: "Sources/NumiVivoCLI")
    ],
    swiftLanguageModes: [.v6],
    // SwiftPM's manifest API calls its C++23 compiler mode cxx2b.
    cxxLanguageStandard: .cxx2b
)
