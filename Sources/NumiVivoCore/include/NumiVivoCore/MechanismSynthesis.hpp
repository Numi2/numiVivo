#pragma once

#include "NumiVivoCore/Core.hpp"

#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace nvivo::mechanism {

enum class Role : std::uint8_t {
    sensor = 0,
    transducer = 1,
    logic = 2,
    temporal = 3,
    memory = 4,
    effector = 5,
    communication = 6,
    containment = 7,
    shutdown = 8,
    monitor = 9
};

enum class EvidenceTier : std::uint8_t {
    observedInTargetContext = 0,
    observedInRelatedContext = 1,
    calibrated = 2,
    inferred = 3,
    assumed = 4,
    hypothetical = 5
};

enum class Reversibility : std::uint8_t {
    reversible = 0,
    resettable = 1,
    conditionallyIrreversible = 2,
    irreversible = 3
};

struct PerformanceEnvelope {
    double payloadBytes = 0.0;
    double cellularBurden = 0.0;
    double latencySeconds = 0.0;
    double leakProbability = 0.0;
    double dynamicRange = 1.0;
    double specificity = 1.0;
    double robustness = 1.0;
    double relativeUncertainty = 0.0;
};

struct Part {
    std::string id;
    Role role = Role::sensor;
    std::set<std::string, std::less<>> acceptedInputs;
    std::set<std::string, std::less<>> producedOutputs;
    std::set<std::string, std::less<>> supportedHosts;
    std::set<std::string, std::less<>> supportedDeliveryModes;
    std::set<std::string, std::less<>> dependencies;
    std::set<std::string, std::less<>> incompatibilities;
    std::string orthogonalityGroup;
    EvidenceTier evidence = EvidenceTier::assumed;
    Reversibility reversibility = Reversibility::reversible;
    bool independentlyControlled = false;
    bool contextInsulated = false;
    bool resourceBuffered = false;
    PerformanceEnvelope performance;
};

struct Slot {
    std::string id;
    Role role = Role::sensor;
    std::string requiredInput;
    std::string requiredOutput;
    bool optional = false;
    bool requireIndependentControl = false;
    Reversibility maximumReversibility = Reversibility::irreversible;
    std::set<std::string, std::less<>> requiredTags;
};

struct ObjectiveWeights {
    double payloadBytes = 1.0;
    double cellularBurden = 1.0;
    double latencySeconds = 1.0;
    double leakProbability = 1.0;
    double uncertainty = 1.0;
    double specificityPenalty = 1.0;
    double robustnessPenalty = 1.0;
    double weakEvidencePenalty = 1.0;
    double irreversiblePenalty = 1.0;
};

struct Constraints {
    std::string host;
    std::string deliveryMode;
    double maximumPayloadBytes = 1.0e12;
    double maximumCellularBurden = 1.0e12;
    double maximumLatencySeconds = 1.0e12;
    double maximumLeakProbability = 1.0;
    double minimumDynamicRange = 1.0;
    double minimumSpecificity = 0.0;
    double minimumRobustness = 0.0;
    double maximumRelativeUncertainty = 1.0e12;
    EvidenceTier weakestPermittedEvidence = EvidenceTier::hypothetical;
    bool requireIndependentShutdown = true;
    bool requireMonitor = true;
    bool requireContextInsulation = false;
    bool requireResourceBuffering = false;
    bool requireDistinctOrthogonalityGroups = false;
    std::uint32_t maximumSolutions = 64;
    std::uint64_t maximumVisitedNodes = 10'000'000;
    ObjectiveWeights weights;
};

struct Problem {
    std::string id;
    std::vector<Slot> slots;
    Constraints constraints;
};

struct CandidateMetrics {
    PerformanceEnvelope aggregate;
    double objective = 0.0;
    EvidenceTier weakestEvidence = EvidenceTier::observedInTargetContext;
    Reversibility worstReversibility = Reversibility::reversible;
};

struct Candidate {
    std::vector<std::string> selectedPartIDs;
    std::map<std::string, std::string, std::less<>> slotAssignments;
    CandidateMetrics metrics;
    std::array<std::uint8_t, kFingerprintBytes> fingerprint{};
};

enum class RejectionReason : std::uint8_t {
    hostIncompatible = 0,
    deliveryIncompatible = 1,
    signalIncompatible = 2,
    evidenceTooWeak = 3,
    reversibilityTooWeak = 4,
    missingIndependentControl = 5,
    incompatibleParts = 6,
    duplicateOrthogonalityGroup = 7,
    missingDependency = 8,
    payloadExceeded = 9,
    burdenExceeded = 10,
    latencyExceeded = 11,
    leakageExceeded = 12,
    uncertaintyExceeded = 13,
    dynamicRangeInsufficient = 14,
    specificityInsufficient = 15,
    robustnessInsufficient = 16,
    requiredSafetyRoleMissing = 17,
    searchBudgetExceeded = 18
};

struct RejectionCount {
    RejectionReason reason = RejectionReason::hostIncompatible;
    std::uint64_t count = 0;
};

struct Result {
    std::vector<Candidate> candidates;
    std::vector<RejectionCount> rejectionCounts;
    std::uint64_t visitedNodes = 0;
    std::uint64_t completedAssignments = 0;
    bool searchBudgetExhausted = false;
    Diagnostics diagnostics;
};

class Synthesizer {
public:
    [[nodiscard]] Result synthesize(const Problem& problem,
                                    std::span<const Part> library) const;
};

[[nodiscard]] bool validatePart(const Part& part, Diagnostics& diagnostics,
                                std::string_view path = {});
[[nodiscard]] bool validateProblem(const Problem& problem, Diagnostics& diagnostics);
[[nodiscard]] std::string_view roleName(Role role) noexcept;
[[nodiscard]] std::string_view evidenceTierName(EvidenceTier tier) noexcept;
[[nodiscard]] std::string_view reversibilityName(Reversibility value) noexcept;
[[nodiscard]] std::string_view rejectionReasonName(RejectionReason reason) noexcept;

} // namespace nvivo::mechanism
