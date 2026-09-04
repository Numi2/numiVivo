#include "NumiVivoCore/MechanismSynthesis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace nvivo::mechanism {

namespace {

constexpr std::size_t kRejectionReasonCount = 19;

bool finite(double value) noexcept {
    return std::isfinite(value);
}

std::size_t reasonIndex(RejectionReason reason) noexcept {
    return static_cast<std::size_t>(reason);
}

void reject(std::array<std::uint64_t, kRejectionReasonCount>& counts,
            RejectionReason reason) {
    ++counts[reasonIndex(reason)];
}

bool supports(const std::set<std::string, std::less<>>& values,
              std::string_view required) {
    return required.empty() || values.empty() || values.contains(required);
}

bool acceptsSignal(const std::set<std::string, std::less<>>& values,
                   std::string_view required) {
    return required.empty() || values.contains(required) || values.contains("*");
}

bool hasRequiredTag(const Part& part, std::string_view tag) {
    if (tag == "independently-controlled") return part.independentlyControlled;
    if (tag == "context-insulated") return part.contextInsulated;
    if (tag == "resource-buffered") return part.resourceBuffered;
    if (tag == "reversible") return part.reversibility == Reversibility::reversible;
    if (tag == "resettable") {
        return static_cast<std::uint8_t>(part.reversibility) <=
               static_cast<std::uint8_t>(Reversibility::resettable);
    }
    return false;
}

double combinedLeak(double current, double next) {
    return 1.0 - (1.0 - std::clamp(current, 0.0, 1.0)) *
                 (1.0 - std::clamp(next, 0.0, 1.0));
}

CandidateMetrics addPart(const CandidateMetrics& current,
                         const Part& part,
                         const ObjectiveWeights& weights) {
    CandidateMetrics next = current;
    next.aggregate.payloadBytes += part.performance.payloadBytes;
    next.aggregate.cellularBurden += part.performance.cellularBurden;
    next.aggregate.latencySeconds += part.performance.latencySeconds;
    next.aggregate.leakProbability = combinedLeak(
        current.aggregate.leakProbability,
        part.performance.leakProbability
    );
    next.aggregate.dynamicRange = std::min(
        current.aggregate.dynamicRange,
        part.performance.dynamicRange
    );
    next.aggregate.specificity = std::min(
        current.aggregate.specificity,
        part.performance.specificity
    );
    next.aggregate.robustness = std::min(
        current.aggregate.robustness,
        part.performance.robustness
    );
    next.aggregate.relativeUncertainty = std::hypot(
        current.aggregate.relativeUncertainty,
        part.performance.relativeUncertainty
    );
    next.weakestEvidence = std::max(current.weakestEvidence, part.evidence);
    next.worstReversibility = std::max(current.worstReversibility, part.reversibility);

    next.objective =
        weights.payloadBytes * next.aggregate.payloadBytes +
        weights.cellularBurden * next.aggregate.cellularBurden +
        weights.latencySeconds * next.aggregate.latencySeconds +
        weights.leakProbability * next.aggregate.leakProbability +
        weights.uncertainty * next.aggregate.relativeUncertainty +
        weights.specificityPenalty * (1.0 - next.aggregate.specificity) +
        weights.robustnessPenalty * (1.0 - next.aggregate.robustness) +
        weights.weakEvidencePenalty * static_cast<double>(next.weakestEvidence) +
        weights.irreversiblePenalty * static_cast<double>(next.worstReversibility);
    return next;
}

CandidateMetrics initialMetrics() {
    CandidateMetrics metrics;
    metrics.aggregate.dynamicRange = std::numeric_limits<double>::infinity();
    metrics.aggregate.specificity = 1.0;
    metrics.aggregate.robustness = 1.0;
    metrics.weakestEvidence = EvidenceTier::observedInTargetContext;
    metrics.worstReversibility = Reversibility::reversible;
    return metrics;
}

std::optional<RejectionReason> violatesUpperBounds(
    const CandidateMetrics& metrics,
    const Constraints& constraints) {
    if (metrics.aggregate.payloadBytes > constraints.maximumPayloadBytes) {
        return RejectionReason::payloadExceeded;
    }
    if (metrics.aggregate.cellularBurden > constraints.maximumCellularBurden) {
        return RejectionReason::burdenExceeded;
    }
    if (metrics.aggregate.latencySeconds > constraints.maximumLatencySeconds) {
        return RejectionReason::latencyExceeded;
    }
    if (metrics.aggregate.leakProbability > constraints.maximumLeakProbability) {
        return RejectionReason::leakageExceeded;
    }
    if (metrics.aggregate.relativeUncertainty > constraints.maximumRelativeUncertainty) {
        return RejectionReason::uncertaintyExceeded;
    }
    if (metrics.aggregate.specificity < constraints.minimumSpecificity) {
        return RejectionReason::specificityInsufficient;
    }
    if (metrics.aggregate.robustness < constraints.minimumRobustness) {
        return RejectionReason::robustnessInsufficient;
    }
    return std::nullopt;
}

std::optional<RejectionReason> violatesFinalBounds(
    const CandidateMetrics& metrics,
    const Constraints& constraints) {
    if (const auto upper = violatesUpperBounds(metrics, constraints)) {
        return upper;
    }
    if (metrics.aggregate.dynamicRange < constraints.minimumDynamicRange) {
        return RejectionReason::dynamicRangeInsufficient;
    }
    return std::nullopt;
}

bool incompatible(const Part& candidate,
                  std::span<const Part* const> selected,
                  bool distinctOrthogonalityGroups,
                  RejectionReason& reason) {
    for (const Part* existing : selected) {
        if (existing == nullptr) continue;
        if (candidate.incompatibilities.contains(existing->id) ||
            existing->incompatibilities.contains(candidate.id)) {
            reason = RejectionReason::incompatibleParts;
            return true;
        }
        if (distinctOrthogonalityGroups &&
            !candidate.orthogonalityGroup.empty() &&
            candidate.orthogonalityGroup == existing->orthogonalityGroup) {
            reason = RejectionReason::duplicateOrthogonalityGroup;
            return true;
        }
    }
    return false;
}

bool allDependenciesSelected(std::span<const Part* const> selected) {
    std::unordered_set<std::string_view> ids;
    ids.reserve(selected.size());
    for (const Part* part : selected) {
        if (part != nullptr) ids.insert(part->id);
    }
    for (const Part* part : selected) {
        if (part == nullptr) continue;
        for (const auto& dependency : part->dependencies) {
            if (!ids.contains(dependency)) return false;
        }
    }
    return true;
}

bool hasRole(std::span<const Part* const> selected, Role role) {
    return std::any_of(selected.begin(), selected.end(), [role](const Part* part) {
        return part != nullptr && part->role == role;
    });
}

bool hasIndependentShutdown(std::span<const Part* const> selected) {
    return std::any_of(selected.begin(), selected.end(), [](const Part* part) {
        return part != nullptr &&
               part->role == Role::shutdown &&
               part->independentlyControlled &&
               static_cast<std::uint8_t>(part->reversibility) <=
                   static_cast<std::uint8_t>(Reversibility::resettable);
    });
}

bool dominates(const Candidate& left, const Candidate& right) {
    const auto& a = left.metrics;
    const auto& b = right.metrics;
    const bool noWorse =
        a.aggregate.payloadBytes <= b.aggregate.payloadBytes &&
        a.aggregate.cellularBurden <= b.aggregate.cellularBurden &&
        a.aggregate.latencySeconds <= b.aggregate.latencySeconds &&
        a.aggregate.leakProbability <= b.aggregate.leakProbability &&
        a.aggregate.relativeUncertainty <= b.aggregate.relativeUncertainty &&
        a.aggregate.dynamicRange >= b.aggregate.dynamicRange &&
        a.aggregate.specificity >= b.aggregate.specificity &&
        a.aggregate.robustness >= b.aggregate.robustness &&
        a.weakestEvidence <= b.weakestEvidence &&
        a.worstReversibility <= b.worstReversibility;
    const bool strictlyBetter =
        a.aggregate.payloadBytes < b.aggregate.payloadBytes ||
        a.aggregate.cellularBurden < b.aggregate.cellularBurden ||
        a.aggregate.latencySeconds < b.aggregate.latencySeconds ||
        a.aggregate.leakProbability < b.aggregate.leakProbability ||
        a.aggregate.relativeUncertainty < b.aggregate.relativeUncertainty ||
        a.aggregate.dynamicRange > b.aggregate.dynamicRange ||
        a.aggregate.specificity > b.aggregate.specificity ||
        a.aggregate.robustness > b.aggregate.robustness ||
        a.weakestEvidence < b.weakestEvidence ||
        a.worstReversibility < b.worstReversibility;
    return noWorse && strictlyBetter;
}

std::array<std::uint8_t, kFingerprintBytes> candidateFingerprint(
    const Problem& problem,
    const std::map<std::string, std::string, std::less<>>& assignments) {
    std::ostringstream stream;
    stream << "numivivo-mechanism-candidate-v1\n" << problem.id << '\n'
           << problem.constraints.host << '\n'
           << problem.constraints.deliveryMode << '\n';
    for (const auto& [slot, part] : assignments) {
        stream << slot << '=' << part << '\n';
    }
    const std::string canonical = stream.str();
    return sha256({
        reinterpret_cast<const std::byte*>(canonical.data()),
        canonical.size()
    });
}

void insertCandidate(std::vector<Candidate>& candidates,
                     Candidate candidate,
                     std::uint32_t maximumSolutions) {
    if (std::any_of(candidates.begin(), candidates.end(), [&](const Candidate& existing) {
            return existing.fingerprint == candidate.fingerprint || dominates(existing, candidate);
        })) {
        return;
    }
    std::erase_if(candidates, [&](const Candidate& existing) {
        return dominates(candidate, existing);
    });
    candidates.push_back(std::move(candidate));
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& left, const Candidate& right) {
        if (left.metrics.objective != right.metrics.objective) {
            return left.metrics.objective < right.metrics.objective;
        }
        return left.fingerprint < right.fingerprint;
    });
    if (candidates.size() > maximumSolutions) {
        candidates.resize(maximumSolutions);
    }
}

} // namespace

bool validatePart(const Part& part,
                  Diagnostics& diagnostics,
                  std::string_view path) {
    bool valid = true;
    const std::string base(path);
    if (!isValidIdentifier(part.id)) {
        diagnostics.error("NVM001", "Mechanism part has an invalid identifier.", base + ".id");
        valid = false;
    }
    const auto& metrics = part.performance;
    const std::array values = {
        metrics.payloadBytes,
        metrics.cellularBurden,
        metrics.latencySeconds,
        metrics.leakProbability,
        metrics.dynamicRange,
        metrics.specificity,
        metrics.robustness,
        metrics.relativeUncertainty
    };
    if (!std::all_of(values.begin(), values.end(), finite)) {
        diagnostics.error("NVM002", "Mechanism performance envelope contains a non-finite value.", base + ".performance");
        valid = false;
    }
    if (metrics.payloadBytes < 0.0 || metrics.cellularBurden < 0.0 ||
        metrics.latencySeconds < 0.0 || metrics.leakProbability < 0.0 ||
        metrics.leakProbability > 1.0 || metrics.dynamicRange < 1.0 ||
        metrics.specificity < 0.0 || metrics.specificity > 1.0 ||
        metrics.robustness < 0.0 || metrics.robustness > 1.0 ||
        metrics.relativeUncertainty < 0.0) {
        diagnostics.error("NVM003", "Mechanism performance envelope is outside valid bounds.", base + ".performance");
        valid = false;
    }
    if (part.role != Role::sensor && part.producedOutputs.empty() &&
        part.role != Role::containment && part.role != Role::shutdown &&
        part.role != Role::monitor) {
        diagnostics.warning("NVM004", "Mechanism part does not declare an output signal.", base + ".producedOutputs");
    }
    if (part.dependencies.contains(part.id) || part.incompatibilities.contains(part.id)) {
        diagnostics.error("NVM005", "Mechanism part cannot depend on or conflict with itself.", base);
        valid = false;
    }
    return valid;
}

bool validateProblem(const Problem& problem, Diagnostics& diagnostics) {
    bool valid = true;
    if (!isValidIdentifier(problem.id)) {
        diagnostics.error("NVM006", "Mechanism synthesis problem has an invalid identifier.", "$.id");
        valid = false;
    }
    if (problem.slots.empty()) {
        diagnostics.error("NVM007", "Mechanism synthesis problem has no slots.", "$.slots");
        valid = false;
    }
    std::set<std::string, std::less<>> slotIDs;
    for (std::size_t index = 0; index < problem.slots.size(); ++index) {
        const auto& slot = problem.slots[index];
        const std::string path = "$.slots[" + std::to_string(index) + "]";
        if (!isValidIdentifier(slot.id)) {
            diagnostics.error("NVM008", "Mechanism slot has an invalid identifier.", path + ".id");
            valid = false;
        }
        if (!slotIDs.insert(slot.id).second) {
            diagnostics.error("NVM009", "Mechanism slot identifier is duplicated.", path + ".id");
            valid = false;
        }
        for (const auto& tag : slot.requiredTags) {
            if (tag != "independently-controlled" && tag != "context-insulated" &&
                tag != "resource-buffered" && tag != "reversible" && tag != "resettable") {
                diagnostics.error("NVM010", "Mechanism slot contains an unsupported required tag '" + tag + "'.", path + ".requiredTags");
                valid = false;
            }
        }
    }

    const auto& constraints = problem.constraints;
    const std::array nonnegative = {
        constraints.maximumPayloadBytes,
        constraints.maximumCellularBurden,
        constraints.maximumLatencySeconds,
        constraints.maximumRelativeUncertainty
    };
    if (!std::all_of(nonnegative.begin(), nonnegative.end(), [](double value) {
            return finite(value) && value >= 0.0;
        }) ||
        !finite(constraints.maximumLeakProbability) ||
        constraints.maximumLeakProbability < 0.0 ||
        constraints.maximumLeakProbability > 1.0 ||
        !finite(constraints.minimumDynamicRange) ||
        constraints.minimumDynamicRange < 1.0 ||
        !finite(constraints.minimumSpecificity) ||
        constraints.minimumSpecificity < 0.0 ||
        constraints.minimumSpecificity > 1.0 ||
        !finite(constraints.minimumRobustness) ||
        constraints.minimumRobustness < 0.0 ||
        constraints.minimumRobustness > 1.0 ||
        constraints.maximumSolutions == 0 ||
        constraints.maximumVisitedNodes == 0) {
        diagnostics.error("NVM011", "Mechanism synthesis constraints are invalid.", "$.constraints");
        valid = false;
    }
    const auto& weights = constraints.weights;
    const std::array weightValues = {
        weights.payloadBytes,
        weights.cellularBurden,
        weights.latencySeconds,
        weights.leakProbability,
        weights.uncertainty,
        weights.specificityPenalty,
        weights.robustnessPenalty,
        weights.weakEvidencePenalty,
        weights.irreversiblePenalty
    };
    if (!std::all_of(weightValues.begin(), weightValues.end(), [](double value) {
            return finite(value) && value >= 0.0;
        })) {
        diagnostics.error("NVM012", "Mechanism objective weights must be finite and non-negative.", "$.constraints.weights");
        valid = false;
    }
    return valid;
}

Result Synthesizer::synthesize(const Problem& problem,
                               std::span<const Part> library) const {
    Result result;
    std::array<std::uint64_t, kRejectionReasonCount> rejectionCounts{};
    if (!validateProblem(problem, result.diagnostics)) return result;
    if (library.empty()) {
        result.diagnostics.error("NVM013", "Mechanism library is empty.", "$.library");
        return result;
    }

    std::unordered_map<std::string_view, const Part*> byID;
    byID.reserve(library.size());
    for (std::size_t index = 0; index < library.size(); ++index) {
        const Part& part = library[index];
        if (!validatePart(part, result.diagnostics, "$.library[" + std::to_string(index) + "]")) {
            continue;
        }
        if (!byID.emplace(part.id, &part).second) {
            result.diagnostics.error("NVM014", "Mechanism part identifier is duplicated.", "$.library[" + std::to_string(index) + "].id");
        }
    }
    if (result.diagnostics.hasErrors()) return result;

    struct SlotCandidates {
        const Slot* slot = nullptr;
        std::vector<const Part*> parts;
        std::size_t originalIndex = 0;
    };
    std::vector<SlotCandidates> choices;
    choices.reserve(problem.slots.size());

    for (std::size_t slotIndex = 0; slotIndex < problem.slots.size(); ++slotIndex) {
        const Slot& slot = problem.slots[slotIndex];
        SlotCandidates candidates{&slot, {}, slotIndex};
        for (const Part& part : library) {
            if (part.role != slot.role) continue;
            if (!supports(part.supportedHosts, problem.constraints.host)) {
                reject(rejectionCounts, RejectionReason::hostIncompatible);
                continue;
            }
            if (!supports(part.supportedDeliveryModes, problem.constraints.deliveryMode)) {
                reject(rejectionCounts, RejectionReason::deliveryIncompatible);
                continue;
            }
            if (!acceptsSignal(part.acceptedInputs, slot.requiredInput) ||
                !acceptsSignal(part.producedOutputs, slot.requiredOutput)) {
                reject(rejectionCounts, RejectionReason::signalIncompatible);
                continue;
            }
            if (part.evidence > problem.constraints.weakestPermittedEvidence) {
                reject(rejectionCounts, RejectionReason::evidenceTooWeak);
                continue;
            }
            if (part.reversibility > slot.maximumReversibility) {
                reject(rejectionCounts, RejectionReason::reversibilityTooWeak);
                continue;
            }
            if (slot.requireIndependentControl && !part.independentlyControlled) {
                reject(rejectionCounts, RejectionReason::missingIndependentControl);
                continue;
            }
            if (problem.constraints.requireContextInsulation && !part.contextInsulated) {
                reject(rejectionCounts, RejectionReason::signalIncompatible);
                continue;
            }
            if (problem.constraints.requireResourceBuffering && !part.resourceBuffered) {
                reject(rejectionCounts, RejectionReason::signalIncompatible);
                continue;
            }
            if (!std::all_of(slot.requiredTags.begin(), slot.requiredTags.end(), [&](const std::string& tag) {
                    return hasRequiredTag(part, tag);
                })) {
                reject(rejectionCounts, RejectionReason::signalIncompatible);
                continue;
            }
            candidates.parts.push_back(&part);
        }
        std::sort(candidates.parts.begin(), candidates.parts.end(), [](const Part* left, const Part* right) {
            return left->id < right->id;
        });
        if (candidates.parts.empty() && !slot.optional) {
            result.diagnostics.error(
                "NVM015",
                "No mechanism part satisfies required slot '" + slot.id + "'.",
                "$.slots[" + std::to_string(slotIndex) + "]"
            );
        }
        choices.push_back(std::move(candidates));
    }
    if (result.diagnostics.hasErrors()) return result;

    std::stable_sort(choices.begin(), choices.end(), [](const SlotCandidates& left, const SlotCandidates& right) {
        const std::size_t leftCount = left.parts.size() + (left.slot->optional ? 1 : 0);
        const std::size_t rightCount = right.parts.size() + (right.slot->optional ? 1 : 0);
        if (leftCount != rightCount) return leftCount < rightCount;
        return left.slot->id < right.slot->id;
    });

    std::vector<const Part*> selected;
    selected.reserve(choices.size());
    std::map<std::string, std::string, std::less<>> assignments;

    std::function<void(std::size_t, const CandidateMetrics&)> search;
    search = [&](std::size_t depth, const CandidateMetrics& metrics) {
        if (result.searchBudgetExhausted) return;
        if (++result.visitedNodes > problem.constraints.maximumVisitedNodes) {
            result.searchBudgetExhausted = true;
            reject(rejectionCounts, RejectionReason::searchBudgetExceeded);
            return;
        }

        if (result.candidates.size() >= problem.constraints.maximumSolutions &&
            metrics.objective >= result.candidates.back().metrics.objective) {
            return;
        }

        if (depth == choices.size()) {
            ++result.completedAssignments;
            if (!allDependenciesSelected(selected)) {
                reject(rejectionCounts, RejectionReason::missingDependency);
                return;
            }
            if (problem.constraints.requireIndependentShutdown &&
                !hasIndependentShutdown(selected)) {
                reject(rejectionCounts, RejectionReason::requiredSafetyRoleMissing);
                return;
            }
            if (problem.constraints.requireMonitor && !hasRole(selected, Role::monitor)) {
                reject(rejectionCounts, RejectionReason::requiredSafetyRoleMissing);
                return;
            }
            if (const auto violation = violatesFinalBounds(metrics, problem.constraints)) {
                reject(rejectionCounts, *violation);
                return;
            }

            Candidate candidate;
            candidate.slotAssignments = assignments;
            candidate.metrics = metrics;
            candidate.fingerprint = candidateFingerprint(problem, assignments);
            for (const Part* part : selected) {
                if (part != nullptr) candidate.selectedPartIDs.push_back(part->id);
            }
            std::sort(candidate.selectedPartIDs.begin(), candidate.selectedPartIDs.end());
            insertCandidate(result.candidates, std::move(candidate), problem.constraints.maximumSolutions);
            return;
        }

        const SlotCandidates& current = choices[depth];
        if (current.slot->optional) {
            assignments[current.slot->id] = "";
            selected.push_back(nullptr);
            search(depth + 1, metrics);
            selected.pop_back();
            assignments.erase(current.slot->id);
        }

        for (const Part* part : current.parts) {
            RejectionReason pairReason = RejectionReason::incompatibleParts;
            if (incompatible(
                    *part,
                    selected,
                    problem.constraints.requireDistinctOrthogonalityGroups,
                    pairReason)) {
                reject(rejectionCounts, pairReason);
                continue;
            }
            const CandidateMetrics nextMetrics = addPart(
                metrics,
                *part,
                problem.constraints.weights
            );
            if (const auto violation = violatesUpperBounds(nextMetrics, problem.constraints)) {
                reject(rejectionCounts, *violation);
                continue;
            }
            selected.push_back(part);
            assignments[current.slot->id] = part->id;
            search(depth + 1, nextMetrics);
            assignments.erase(current.slot->id);
            selected.pop_back();
            if (result.searchBudgetExhausted) break;
        }
    };

    search(0, initialMetrics());
    if (result.searchBudgetExhausted) {
        result.diagnostics.warning(
            "NVM016",
            "Mechanism synthesis search budget was exhausted; returned candidates are valid but may not contain the global optimum.",
            "$.constraints.maximumVisitedNodes"
        );
    }
    if (result.candidates.empty()) {
        result.diagnostics.warning(
            "NVM017",
            "No complete mechanism assignment satisfies all declared constraints.",
            "$"
        );
    }

    for (std::size_t index = 0; index < rejectionCounts.size(); ++index) {
        if (rejectionCounts[index] == 0) continue;
        result.rejectionCounts.push_back({
            static_cast<RejectionReason>(index),
            rejectionCounts[index]
        });
    }
    return result;
}

std::string_view roleName(Role role) noexcept {
    switch (role) {
        case Role::sensor: return "sensor";
        case Role::transducer: return "transducer";
        case Role::logic: return "logic";
        case Role::temporal: return "temporal";
        case Role::memory: return "memory";
        case Role::effector: return "effector";
        case Role::communication: return "communication";
        case Role::containment: return "containment";
        case Role::shutdown: return "shutdown";
        case Role::monitor: return "monitor";
    }
    return "unknown";
}

std::string_view evidenceTierName(EvidenceTier tier) noexcept {
    switch (tier) {
        case EvidenceTier::observedInTargetContext: return "observed-target-context";
        case EvidenceTier::observedInRelatedContext: return "observed-related-context";
        case EvidenceTier::calibrated: return "calibrated";
        case EvidenceTier::inferred: return "inferred";
        case EvidenceTier::assumed: return "assumed";
        case EvidenceTier::hypothetical: return "hypothetical";
    }
    return "unknown";
}

std::string_view reversibilityName(Reversibility value) noexcept {
    switch (value) {
        case Reversibility::reversible: return "reversible";
        case Reversibility::resettable: return "resettable";
        case Reversibility::conditionallyIrreversible: return "conditionally-irreversible";
        case Reversibility::irreversible: return "irreversible";
    }
    return "unknown";
}

std::string_view rejectionReasonName(RejectionReason reason) noexcept {
    switch (reason) {
        case RejectionReason::hostIncompatible: return "host-incompatible";
        case RejectionReason::deliveryIncompatible: return "delivery-incompatible";
        case RejectionReason::signalIncompatible: return "signal-incompatible";
        case RejectionReason::evidenceTooWeak: return "evidence-too-weak";
        case RejectionReason::reversibilityTooWeak: return "reversibility-too-weak";
        case RejectionReason::missingIndependentControl: return "missing-independent-control";
        case RejectionReason::incompatibleParts: return "incompatible-parts";
        case RejectionReason::duplicateOrthogonalityGroup: return "duplicate-orthogonality-group";
        case RejectionReason::missingDependency: return "missing-dependency";
        case RejectionReason::payloadExceeded: return "payload-exceeded";
        case RejectionReason::burdenExceeded: return "burden-exceeded";
        case RejectionReason::latencyExceeded: return "latency-exceeded";
        case RejectionReason::leakageExceeded: return "leakage-exceeded";
        case RejectionReason::uncertaintyExceeded: return "uncertainty-exceeded";
        case RejectionReason::dynamicRangeInsufficient: return "dynamic-range-insufficient";
        case RejectionReason::specificityInsufficient: return "specificity-insufficient";
        case RejectionReason::robustnessInsufficient: return "robustness-insufficient";
        case RejectionReason::requiredSafetyRoleMissing: return "required-safety-role-missing";
        case RejectionReason::searchBudgetExceeded: return "search-budget-exceeded";
    }
    return "unknown";
}

} // namespace nvivo::mechanism
