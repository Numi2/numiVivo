#include "NumiVivoCore/NumiVivoMechanismSynthesis.h"
#include "NumiVivoCore/MechanismSynthesis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace {

using namespace nvivo;
using namespace nvivo::mechanism;

struct DecodeContext {
    Diagnostics diagnostics;
};

void clearBuffer(NVivoByteBuffer* buffer) {
    if (buffer == nullptr) return;
    buffer->data = nullptr;
    buffer->size = 0;
}

bool copyBuffer(std::string_view value, NVivoByteBuffer* output) {
    if (output == nullptr) return true;
    clearBuffer(output);
    if (value.empty()) return true;
    void* memory = std::malloc(value.size());
    if (memory == nullptr) return false;
    std::memcpy(memory, value.data(), value.size());
    output->data = static_cast<std::uint8_t*>(memory);
    output->size = value.size();
    return true;
}

void writeString(std::ostringstream& stream, std::string_view value) {
    stream << '"' << json::escape(value) << '"';
}

std::string pathJoin(std::string_view base, std::string_view component) {
    return std::string(base) + "." + std::string(component);
}

std::string pathIndex(std::string_view base, std::size_t index) {
    return std::string(base) + "[" + std::to_string(index) + "]";
}

bool requireObject(const json::Value& value,
                   std::string_view path,
                   DecodeContext& context) {
    if (value.isObject()) return true;
    context.diagnostics.error("NVMJ001", "Expected a JSON object.", std::string(path));
    return false;
}

const json::Value* optionalMember(const json::Value& object,
                                  std::string_view key) {
    return object.get(key);
}

const json::Value* requireMember(const json::Value& object,
                                 std::string_view key,
                                 std::string_view path,
                                 DecodeContext& context) {
    const auto* value = object.get(key);
    if (value == nullptr) {
        context.diagnostics.error("NVMJ002", "Required member is missing.", std::string(path));
    }
    return value;
}

std::string stringValue(const json::Value& object,
                        std::string_view key,
                        std::string fallback,
                        std::string_view path,
                        DecodeContext& context,
                        bool required = false) {
    const auto* value = required
        ? requireMember(object, key, path, context)
        : optionalMember(object, key);
    if (value == nullptr) return fallback;
    if (!value->isString()) {
        context.diagnostics.error("NVMJ003", "Expected a string.", std::string(path));
        return fallback;
    }
    return std::string(value->asString());
}

double numberValue(const json::Value& object,
                   std::string_view key,
                   double fallback,
                   std::string_view path,
                   DecodeContext& context) {
    const auto* value = optionalMember(object, key);
    if (value == nullptr) return fallback;
    if (!value->isNumber() || !std::isfinite(value->asNumber())) {
        context.diagnostics.error("NVMJ004", "Expected a finite number.", std::string(path));
        return fallback;
    }
    return value->asNumber();
}

bool boolValue(const json::Value& object,
               std::string_view key,
               bool fallback,
               std::string_view path,
               DecodeContext& context) {
    const auto* value = optionalMember(object, key);
    if (value == nullptr) return fallback;
    if (!value->isBool()) {
        context.diagnostics.error("NVMJ005", "Expected a boolean.", std::string(path));
        return fallback;
    }
    return value->asBool();
}

std::optional<std::uint64_t> unsignedIntegerValue(
    const json::Value& object,
    std::string_view key,
    std::string_view path,
    DecodeContext& context) {
    const auto* value = optionalMember(object, key);
    if (value == nullptr) return std::nullopt;
    const double number = value->asNumber(std::numeric_limits<double>::quiet_NaN());
    if (!value->isNumber() || !std::isfinite(number) || number < 0.0 ||
        std::floor(number) != number ||
        number > static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
        context.diagnostics.error("NVMJ006", "Expected an unsigned integer.", std::string(path));
        return std::nullopt;
    }
    return static_cast<std::uint64_t>(number);
}

std::set<std::string, std::less<>> stringSet(
    const json::Value& object,
    std::string_view key,
    std::string_view path,
    DecodeContext& context) {
    std::set<std::string, std::less<>> result;
    const auto* value = optionalMember(object, key);
    if (value == nullptr) return result;
    if (!value->isArray()) {
        context.diagnostics.error("NVMJ007", "Expected an array of strings.", std::string(path));
        return result;
    }
    if (value->asArray().size() > 65'536) {
        context.diagnostics.error("NVMJ008", "String set exceeds 65536 entries.", std::string(path));
        return result;
    }
    for (std::size_t index = 0; index < value->asArray().size(); ++index) {
        const auto& item = value->asArray()[index];
        if (!item.isString() || item.asString().empty() || item.asString().size() > 1'024) {
            context.diagnostics.error(
                "NVMJ009",
                "Set entries must be non-empty bounded strings.",
                pathIndex(path, index)
            );
            continue;
        }
        if (!result.emplace(item.asString()).second) {
            context.diagnostics.error(
                "NVMJ010",
                "String set contains a duplicate entry.",
                pathIndex(path, index)
            );
        }
    }
    return result;
}

Role parseRole(std::string_view value,
               std::string_view path,
               DecodeContext& context) {
    if (value == "sensor") return Role::sensor;
    if (value == "transducer") return Role::transducer;
    if (value == "logic") return Role::logic;
    if (value == "temporal") return Role::temporal;
    if (value == "memory") return Role::memory;
    if (value == "effector") return Role::effector;
    if (value == "communication") return Role::communication;
    if (value == "containment") return Role::containment;
    if (value == "shutdown") return Role::shutdown;
    if (value == "monitor") return Role::monitor;
    context.diagnostics.error(
        "NVMJ011",
        "Unknown mechanism role '" + std::string(value) + "'.",
        std::string(path)
    );
    return Role::sensor;
}

EvidenceTier parseEvidence(std::string_view value,
                           std::string_view path,
                           DecodeContext& context) {
    if (value == "observed-target-context") return EvidenceTier::observedInTargetContext;
    if (value == "observed-related-context") return EvidenceTier::observedInRelatedContext;
    if (value == "calibrated") return EvidenceTier::calibrated;
    if (value == "inferred") return EvidenceTier::inferred;
    if (value == "assumed") return EvidenceTier::assumed;
    if (value == "hypothetical") return EvidenceTier::hypothetical;
    context.diagnostics.error(
        "NVMJ012",
        "Unknown evidence tier '" + std::string(value) + "'.",
        std::string(path)
    );
    return EvidenceTier::assumed;
}

Reversibility parseReversibility(std::string_view value,
                                 std::string_view path,
                                 DecodeContext& context) {
    if (value == "reversible") return Reversibility::reversible;
    if (value == "resettable") return Reversibility::resettable;
    if (value == "conditionally-irreversible") return Reversibility::conditionallyIrreversible;
    if (value == "irreversible") return Reversibility::irreversible;
    context.diagnostics.error(
        "NVMJ013",
        "Unknown reversibility '" + std::string(value) + "'.",
        std::string(path)
    );
    return Reversibility::irreversible;
}

PerformanceEnvelope decodePerformance(const json::Value& value,
                                      std::string_view path,
                                      DecodeContext& context) {
    PerformanceEnvelope result;
    if (!requireObject(value, path, context)) return result;
    result.payloadBytes = numberValue(value, "payloadBytes", 0.0, pathJoin(path, "payloadBytes"), context);
    result.cellularBurden = numberValue(value, "cellularBurden", 0.0, pathJoin(path, "cellularBurden"), context);
    result.latencySeconds = numberValue(value, "latencySeconds", 0.0, pathJoin(path, "latencySeconds"), context);
    result.leakProbability = numberValue(value, "leakProbability", 0.0, pathJoin(path, "leakProbability"), context);
    result.dynamicRange = numberValue(value, "dynamicRange", 1.0, pathJoin(path, "dynamicRange"), context);
    result.specificity = numberValue(value, "specificity", 1.0, pathJoin(path, "specificity"), context);
    result.robustness = numberValue(value, "robustness", 1.0, pathJoin(path, "robustness"), context);
    result.relativeUncertainty = numberValue(value, "relativeUncertainty", 0.0, pathJoin(path, "relativeUncertainty"), context);
    return result;
}

Part decodePart(const json::Value& value,
                std::string_view path,
                DecodeContext& context) {
    Part part;
    if (!requireObject(value, path, context)) return part;
    part.id = stringValue(value, "id", {}, pathJoin(path, "id"), context, true);
    part.role = parseRole(
        stringValue(value, "role", {}, pathJoin(path, "role"), context, true),
        pathJoin(path, "role"),
        context
    );
    part.acceptedInputs = stringSet(value, "acceptedInputs", pathJoin(path, "acceptedInputs"), context);
    part.producedOutputs = stringSet(value, "producedOutputs", pathJoin(path, "producedOutputs"), context);
    part.supportedHosts = stringSet(value, "supportedHosts", pathJoin(path, "supportedHosts"), context);
    part.supportedDeliveryModes = stringSet(value, "supportedDeliveryModes", pathJoin(path, "supportedDeliveryModes"), context);
    part.dependencies = stringSet(value, "dependencies", pathJoin(path, "dependencies"), context);
    part.incompatibilities = stringSet(value, "incompatibilities", pathJoin(path, "incompatibilities"), context);
    part.orthogonalityGroup = stringValue(
        value,
        "orthogonalityGroup",
        {},
        pathJoin(path, "orthogonalityGroup"),
        context
    );
    part.evidence = parseEvidence(
        stringValue(value, "evidence", "assumed", pathJoin(path, "evidence"), context),
        pathJoin(path, "evidence"),
        context
    );
    part.reversibility = parseReversibility(
        stringValue(value, "reversibility", "reversible", pathJoin(path, "reversibility"), context),
        pathJoin(path, "reversibility"),
        context
    );
    part.independentlyControlled = boolValue(
        value,
        "independentlyControlled",
        false,
        pathJoin(path, "independentlyControlled"),
        context
    );
    part.contextInsulated = boolValue(
        value,
        "contextInsulated",
        false,
        pathJoin(path, "contextInsulated"),
        context
    );
    part.resourceBuffered = boolValue(
        value,
        "resourceBuffered",
        false,
        pathJoin(path, "resourceBuffered"),
        context
    );
    const auto* performance = requireMember(value, "performance", pathJoin(path, "performance"), context);
    if (performance != nullptr) {
        part.performance = decodePerformance(*performance, pathJoin(path, "performance"), context);
    }
    return part;
}

std::vector<Part> decodeLibrary(const json::Value& root,
                                DecodeContext& context) {
    std::vector<Part> library;
    if (!requireObject(root, "$library", context)) return library;
    const auto version = unsignedIntegerValue(root, "schemaVersion", "$library.schemaVersion", context);
    if (!version.has_value() || *version != 1) {
        context.diagnostics.error("NVMJ014", "Mechanism library schemaVersion must be 1.", "$library.schemaVersion");
    }
    const auto* parts = requireMember(root, "parts", "$library.parts", context);
    if (parts == nullptr) return library;
    if (!parts->isArray()) {
        context.diagnostics.error("NVMJ015", "Mechanism library parts must be an array.", "$library.parts");
        return library;
    }
    if (parts->asArray().empty() || parts->asArray().size() > 1'000'000) {
        context.diagnostics.error("NVMJ016", "Mechanism library part count is outside 1...1000000.", "$library.parts");
        return library;
    }
    library.reserve(parts->asArray().size());
    for (std::size_t index = 0; index < parts->asArray().size(); ++index) {
        library.push_back(decodePart(parts->asArray()[index], pathIndex("$library.parts", index), context));
    }
    return library;
}

Slot decodeSlot(const json::Value& value,
                std::string_view path,
                DecodeContext& context) {
    Slot slot;
    if (!requireObject(value, path, context)) return slot;
    slot.id = stringValue(value, "id", {}, pathJoin(path, "id"), context, true);
    slot.role = parseRole(
        stringValue(value, "role", {}, pathJoin(path, "role"), context, true),
        pathJoin(path, "role"),
        context
    );
    slot.requiredInput = stringValue(value, "requiredInput", {}, pathJoin(path, "requiredInput"), context);
    slot.requiredOutput = stringValue(value, "requiredOutput", {}, pathJoin(path, "requiredOutput"), context);
    slot.optional = boolValue(value, "optional", false, pathJoin(path, "optional"), context);
    slot.requireIndependentControl = boolValue(
        value,
        "requireIndependentControl",
        false,
        pathJoin(path, "requireIndependentControl"),
        context
    );
    slot.maximumReversibility = parseReversibility(
        stringValue(
            value,
            "maximumReversibility",
            "irreversible",
            pathJoin(path, "maximumReversibility"),
            context
        ),
        pathJoin(path, "maximumReversibility"),
        context
    );
    slot.requiredTags = stringSet(value, "requiredTags", pathJoin(path, "requiredTags"), context);
    return slot;
}

ObjectiveWeights decodeWeights(const json::Value* value,
                               std::string_view path,
                               DecodeContext& context) {
    ObjectiveWeights weights;
    if (value == nullptr) return weights;
    if (!requireObject(*value, path, context)) return weights;
    weights.payloadBytes = numberValue(*value, "payloadBytes", weights.payloadBytes, pathJoin(path, "payloadBytes"), context);
    weights.cellularBurden = numberValue(*value, "cellularBurden", weights.cellularBurden, pathJoin(path, "cellularBurden"), context);
    weights.latencySeconds = numberValue(*value, "latencySeconds", weights.latencySeconds, pathJoin(path, "latencySeconds"), context);
    weights.leakProbability = numberValue(*value, "leakProbability", weights.leakProbability, pathJoin(path, "leakProbability"), context);
    weights.uncertainty = numberValue(*value, "uncertainty", weights.uncertainty, pathJoin(path, "uncertainty"), context);
    weights.specificityPenalty = numberValue(*value, "specificityPenalty", weights.specificityPenalty, pathJoin(path, "specificityPenalty"), context);
    weights.robustnessPenalty = numberValue(*value, "robustnessPenalty", weights.robustnessPenalty, pathJoin(path, "robustnessPenalty"), context);
    weights.weakEvidencePenalty = numberValue(*value, "weakEvidencePenalty", weights.weakEvidencePenalty, pathJoin(path, "weakEvidencePenalty"), context);
    weights.irreversiblePenalty = numberValue(*value, "irreversiblePenalty", weights.irreversiblePenalty, pathJoin(path, "irreversiblePenalty"), context);
    return weights;
}

Constraints decodeConstraints(const json::Value* value,
                              DecodeContext& context) {
    Constraints constraints;
    if (value == nullptr) {
        context.diagnostics.error("NVMJ017", "Problem constraints are required.", "$problem.constraints");
        return constraints;
    }
    if (!requireObject(*value, "$problem.constraints", context)) return constraints;
    constraints.host = stringValue(*value, "host", {}, "$problem.constraints.host", context, true);
    constraints.deliveryMode = stringValue(*value, "deliveryMode", {}, "$problem.constraints.deliveryMode", context, true);
    constraints.maximumPayloadBytes = numberValue(*value, "maximumPayloadBytes", constraints.maximumPayloadBytes, "$problem.constraints.maximumPayloadBytes", context);
    constraints.maximumCellularBurden = numberValue(*value, "maximumCellularBurden", constraints.maximumCellularBurden, "$problem.constraints.maximumCellularBurden", context);
    constraints.maximumLatencySeconds = numberValue(*value, "maximumLatencySeconds", constraints.maximumLatencySeconds, "$problem.constraints.maximumLatencySeconds", context);
    constraints.maximumLeakProbability = numberValue(*value, "maximumLeakProbability", constraints.maximumLeakProbability, "$problem.constraints.maximumLeakProbability", context);
    constraints.minimumDynamicRange = numberValue(*value, "minimumDynamicRange", constraints.minimumDynamicRange, "$problem.constraints.minimumDynamicRange", context);
    constraints.minimumSpecificity = numberValue(*value, "minimumSpecificity", constraints.minimumSpecificity, "$problem.constraints.minimumSpecificity", context);
    constraints.minimumRobustness = numberValue(*value, "minimumRobustness", constraints.minimumRobustness, "$problem.constraints.minimumRobustness", context);
    constraints.maximumRelativeUncertainty = numberValue(*value, "maximumRelativeUncertainty", constraints.maximumRelativeUncertainty, "$problem.constraints.maximumRelativeUncertainty", context);
    constraints.weakestPermittedEvidence = parseEvidence(
        stringValue(*value, "weakestPermittedEvidence", "hypothetical", "$problem.constraints.weakestPermittedEvidence", context),
        "$problem.constraints.weakestPermittedEvidence",
        context
    );
    constraints.requireIndependentShutdown = boolValue(*value, "requireIndependentShutdown", true, "$problem.constraints.requireIndependentShutdown", context);
    constraints.requireMonitor = boolValue(*value, "requireMonitor", true, "$problem.constraints.requireMonitor", context);
    constraints.requireContextInsulation = boolValue(*value, "requireContextInsulation", false, "$problem.constraints.requireContextInsulation", context);
    constraints.requireResourceBuffering = boolValue(*value, "requireResourceBuffering", false, "$problem.constraints.requireResourceBuffering", context);
    constraints.requireDistinctOrthogonalityGroups = boolValue(*value, "requireDistinctOrthogonalityGroups", false, "$problem.constraints.requireDistinctOrthogonalityGroups", context);
    if (const auto maximumSolutions = unsignedIntegerValue(*value, "maximumSolutions", "$problem.constraints.maximumSolutions", context)) {
        if (*maximumSolutions == 0 || *maximumSolutions > std::numeric_limits<std::uint32_t>::max()) {
            context.diagnostics.error("NVMJ018", "maximumSolutions is outside the UInt32 range.", "$problem.constraints.maximumSolutions");
        } else {
            constraints.maximumSolutions = static_cast<std::uint32_t>(*maximumSolutions);
        }
    }
    if (const auto maximumNodes = unsignedIntegerValue(*value, "maximumVisitedNodes", "$problem.constraints.maximumVisitedNodes", context)) {
        constraints.maximumVisitedNodes = *maximumNodes;
    }
    constraints.weights = decodeWeights(
        optionalMember(*value, "weights"),
        "$problem.constraints.weights",
        context
    );
    return constraints;
}

Problem decodeProblem(const json::Value& root,
                      DecodeContext& context) {
    Problem problem;
    if (!requireObject(root, "$problem", context)) return problem;
    const auto version = unsignedIntegerValue(root, "schemaVersion", "$problem.schemaVersion", context);
    if (!version.has_value() || *version != 1) {
        context.diagnostics.error("NVMJ019", "Mechanism problem schemaVersion must be 1.", "$problem.schemaVersion");
    }
    problem.id = stringValue(root, "id", {}, "$problem.id", context, true);
    problem.constraints = decodeConstraints(optionalMember(root, "constraints"), context);
    const auto* slots = requireMember(root, "slots", "$problem.slots", context);
    if (slots == nullptr) return problem;
    if (!slots->isArray()) {
        context.diagnostics.error("NVMJ020", "Problem slots must be an array.", "$problem.slots");
        return problem;
    }
    if (slots->asArray().empty() || slots->asArray().size() > 65'536) {
        context.diagnostics.error("NVMJ021", "Problem slot count is outside 1...65536.", "$problem.slots");
        return problem;
    }
    problem.slots.reserve(slots->asArray().size());
    for (std::size_t index = 0; index < slots->asArray().size(); ++index) {
        problem.slots.push_back(decodeSlot(slots->asArray()[index], pathIndex("$problem.slots", index), context));
    }
    return problem;
}

bool decodeOptions(const NVivoMechanismSynthesisOptions* source,
                   Constraints& constraints,
                   Diagnostics& diagnostics) {
    if (source == nullptr) return true;
    if (source->struct_size < sizeof(NVivoMechanismSynthesisOptions) ||
        source->abi_version != NVIVO_MECHANISM_SYNTHESIS_ABI_VERSION ||
        source->maximum_solutions == 0 ||
        source->maximum_visited_nodes == 0) {
        diagnostics.error("NVMJ022", "Mechanism synthesis ABI options are invalid.", "$options");
        return false;
    }
    constraints.maximumSolutions = source->maximum_solutions;
    constraints.maximumVisitedNodes = source->maximum_visited_nodes;
    constraints.requireIndependentShutdown =
        (source->flags & NVIVO_SYNTHESIS_REQUIRE_INDEPENDENT_SHUTDOWN) != 0;
    constraints.requireMonitor =
        (source->flags & NVIVO_SYNTHESIS_REQUIRE_MONITOR) != 0;
    constraints.requireContextInsulation =
        (source->flags & NVIVO_SYNTHESIS_REQUIRE_CONTEXT_INSULATION) != 0;
    constraints.requireResourceBuffering =
        (source->flags & NVIVO_SYNTHESIS_REQUIRE_RESOURCE_BUFFERING) != 0;
    constraints.requireDistinctOrthogonalityGroups =
        (source->flags & NVIVO_SYNTHESIS_REQUIRE_DISTINCT_ORTHOGONALITY) != 0;
    return true;
}

std::string encodeResult(const Result& result) {
    std::ostringstream stream;
    stream << "{\"schemaVersion\":1"
           << ",\"visitedNodes\":" << result.visitedNodes
           << ",\"completedAssignments\":" << result.completedAssignments
           << ",\"searchBudgetExhausted\":" << (result.searchBudgetExhausted ? "true" : "false")
           << ",\"candidates\":[";
    for (std::size_t candidateIndex = 0; candidateIndex < result.candidates.size(); ++candidateIndex) {
        if (candidateIndex != 0) stream << ',';
        const Candidate& candidate = result.candidates[candidateIndex];
        stream << "{\"fingerprint\":";
        writeString(stream, hexFingerprint(candidate.fingerprint));
        stream << ",\"objective\":" << candidate.metrics.objective
               << ",\"weakestEvidence\":";
        writeString(stream, evidenceTierName(candidate.metrics.weakestEvidence));
        stream << ",\"worstReversibility\":";
        writeString(stream, reversibilityName(candidate.metrics.worstReversibility));
        const auto& metrics = candidate.metrics.aggregate;
        stream << ",\"metrics\":{"
               << "\"payloadBytes\":" << metrics.payloadBytes
               << ",\"cellularBurden\":" << metrics.cellularBurden
               << ",\"latencySeconds\":" << metrics.latencySeconds
               << ",\"leakProbability\":" << metrics.leakProbability
               << ",\"dynamicRange\":" << metrics.dynamicRange
               << ",\"specificity\":" << metrics.specificity
               << ",\"robustness\":" << metrics.robustness
               << ",\"relativeUncertainty\":" << metrics.relativeUncertainty
               << "},\"selectedPartIDs\":[";
        for (std::size_t index = 0; index < candidate.selectedPartIDs.size(); ++index) {
            if (index != 0) stream << ',';
            writeString(stream, candidate.selectedPartIDs[index]);
        }
        stream << "],\"slotAssignments\":{\"";
        bool firstAssignment = true;
        for (const auto& [slot, part] : candidate.slotAssignments) {
            if (!firstAssignment) stream << ',';
            firstAssignment = false;
            writeString(stream, slot);
            stream << ':';
            if (part.empty()) stream << "null";
            else writeString(stream, part);
        }
        stream << "}}";
    }
    stream << "],\"rejectionCounts\":[";
    for (std::size_t index = 0; index < result.rejectionCounts.size(); ++index) {
        if (index != 0) stream << ',';
        stream << "{\"reason\":";
        writeString(stream, rejectionReasonName(result.rejectionCounts[index].reason));
        stream << ",\"count\":" << result.rejectionCounts[index].count << '}';
    }
    stream << "]}";
    return stream.str();
}

NVivoStatus statusFor(const Diagnostics& diagnostics) {
    return diagnostics.hasErrors() ? NVIVO_STATUS_VALIDATION_ERROR : NVIVO_STATUS_OK;
}

} // namespace

extern "C" {

void nvivo_default_mechanism_synthesis_options(
    NVivoMechanismSynthesisOptions* options) {
    if (options == nullptr) return;
    std::memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->abi_version = NVIVO_MECHANISM_SYNTHESIS_ABI_VERSION;
    options->flags = NVIVO_SYNTHESIS_REQUIRE_INDEPENDENT_SHUTDOWN |
                     NVIVO_SYNTHESIS_REQUIRE_MONITOR;
    options->maximum_solutions = 64;
    options->maximum_visited_nodes = 10'000'000;
}

NVivoStatus nvivo_synthesize_mechanisms_json(
    const std::uint8_t* problemJson,
    std::size_t problemSize,
    const std::uint8_t* libraryJson,
    std::size_t librarySize,
    const NVivoMechanismSynthesisOptions* options,
    NVivoByteBuffer* resultJson,
    NVivoByteBuffer* diagnosticsJson) {
    clearBuffer(resultJson);
    clearBuffer(diagnosticsJson);
    if (resultJson == nullptr ||
        (problemJson == nullptr && problemSize != 0) ||
        (libraryJson == nullptr && librarySize != 0)) {
        return NVIVO_STATUS_INVALID_ARGUMENT;
    }

    try {
        const std::string_view problemText(
            problemJson == nullptr ? "" : reinterpret_cast<const char*>(problemJson),
            problemSize
        );
        const std::string_view libraryText(
            libraryJson == nullptr ? "" : reinterpret_cast<const char*>(libraryJson),
            librarySize
        );
        auto parsedProblem = json::parse(problemText);
        auto parsedLibrary = json::parse(libraryText);
        Diagnostics diagnostics;
        diagnostics.append(parsedProblem.diagnostics);
        diagnostics.append(parsedLibrary.diagnostics);
        if (!parsedProblem.root.has_value() || !parsedLibrary.root.has_value()) {
            const std::string report = diagnostics.toJson();
            if (!copyBuffer(report, diagnosticsJson)) return NVIVO_STATUS_OUT_OF_MEMORY;
            return NVIVO_STATUS_PARSE_ERROR;
        }

        DecodeContext decoding;
        Problem problem = decodeProblem(*parsedProblem.root, decoding);
        std::vector<Part> library = decodeLibrary(*parsedLibrary.root, decoding);
        diagnostics.append(decoding.diagnostics);
        if (!decodeOptions(options, problem.constraints, diagnostics) || diagnostics.hasErrors()) {
            const std::string report = diagnostics.toJson();
            if (!copyBuffer(report, diagnosticsJson)) return NVIVO_STATUS_OUT_OF_MEMORY;
            return NVIVO_STATUS_VALIDATION_ERROR;
        }

        Result result = Synthesizer().synthesize(problem, library);
        diagnostics.append(result.diagnostics);
        const std::string resultText = encodeResult(result);
        const std::string diagnosticText = diagnostics.toJson();
        if (!copyBuffer(resultText, resultJson)) return NVIVO_STATUS_OUT_OF_MEMORY;
        if (!copyBuffer(diagnosticText, diagnosticsJson)) {
            nvivo_buffer_release(resultJson);
            return NVIVO_STATUS_OUT_OF_MEMORY;
        }
        return statusFor(diagnostics);
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& error) {
        const std::string report = std::string("{\"status\":\"internal-error\",\"message\":\"") +
                                   json::escape(error.what()) + "\"}";
        if (!copyBuffer(report, diagnosticsJson)) return NVIVO_STATUS_OUT_OF_MEMORY;
        return NVIVO_STATUS_INTERNAL_ERROR;
    } catch (...) {
        const std::string report = "{\"status\":\"internal-error\",\"message\":\"unknown exception\"}";
        if (!copyBuffer(report, diagnosticsJson)) return NVIVO_STATUS_OUT_OF_MEMORY;
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}

} // extern "C"
