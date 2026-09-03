#include "NumiVivoCore/NumiVivoPackAccess.h"
#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <span>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using namespace nvivo;

struct SectionView {
    const PackSectionDescriptor* descriptor = nullptr;
    std::span<const std::byte> bytes;
};

SectionView section(std::span<const std::byte> pack,
                    const PackInspection& inspection,
                    PackSectionType type) {
    const auto rawType = static_cast<std::uint32_t>(type);
    const auto iterator = std::find_if(
        inspection.sections.begin(),
        inspection.sections.end(),
        [rawType](const PackSectionDescriptor& candidate) {
            return candidate.type == rawType;
        }
    );
    if (iterator == inspection.sections.end()) return {};
    return {
        &*iterator,
        pack.subspan(
            static_cast<std::size_t>(iterator->offset),
            static_cast<std::size_t>(iterator->size)
        )
    };
}

template <typename T>
bool loadRecord(const SectionView& view, std::uint32_t index, T& value) {
    if (view.descriptor == nullptr || view.descriptor->stride != sizeof(T) ||
        index >= view.descriptor->count) {
        return false;
    }
    const std::size_t offset = static_cast<std::size_t>(index) * sizeof(T);
    if (offset > view.bytes.size() || sizeof(T) > view.bytes.size() - offset) {
        return false;
    }
    std::memcpy(&value, view.bytes.data() + offset, sizeof(T));
    return true;
}

bool checkedRange(std::uint32_t offset,
                  std::uint32_t count,
                  std::uint32_t capacity) {
    return offset <= capacity && count <= capacity - offset;
}

bool validStringOffset(const SectionView& strings, std::uint32_t offset) {
    if (strings.descriptor == nullptr || offset >= strings.bytes.size()) return false;
    const auto* characters = reinterpret_cast<const char*>(strings.bytes.data());
    for (std::size_t index = offset; index < strings.bytes.size(); ++index) {
        if (characters[index] == '\0') return true;
    }
    return false;
}

void requireString(const SectionView& strings,
                   std::uint32_t offset,
                   std::string_view path,
                   Diagnostics& diagnostics) {
    if (!validStringOffset(strings, offset)) {
        diagnostics.error(
            "NVP100",
            "ProgramPack string offset is invalid or not NUL terminated.",
            std::string(path)
        );
    }
}

bool finiteOrdered(double minimum, double value, double maximum) {
    return std::isfinite(minimum) && std::isfinite(value) && std::isfinite(maximum) &&
           minimum <= value && value <= maximum;
}

struct SemanticContext {
    std::span<const std::byte> pack;
    PackInspection inspection;
    SectionView strings;
    SectionView species;
    SectionView parameters;
    SectionView reactionParameterIndices;
    SectionView stoichiometry;
    SectionView reactions;
    SectionView expressions;
    SectionView actions;
    SectionView rules;
    SectionView monitors;
    SectionView cohorts;
    SectionView incidenceOffsets;
    SectionView incidence;
    SectionView runtimeContract;
    RuntimeContractRecord contract{};
};

std::uint64_t pairKey(std::uint32_t species, std::uint32_t reaction) {
    return (static_cast<std::uint64_t>(species) << 32U) | reaction;
}

bool validateExpression(const SemanticContext& context,
                        std::uint32_t offset,
                        std::uint32_t count,
                        bool scanToEnd,
                        std::string_view path,
                        Diagnostics& diagnostics) {
    if (context.expressions.descriptor == nullptr) return false;
    const std::uint32_t capacity = context.expressions.descriptor->count;
    if (offset >= capacity) {
        diagnostics.error("NVP101", "Expression offset is outside the bytecode section.", std::string(path));
        return false;
    }

    std::uint32_t limit = 0;
    if (scanToEnd) {
        const std::uint32_t remaining = capacity - offset;
        limit = std::min<std::uint32_t>(remaining, 4'096U);
    } else {
        if (count == 0 || !checkedRange(offset, count, capacity)) {
            diagnostics.error("NVP102", "Expression instruction range is invalid.", std::string(path));
            return false;
        }
        limit = count;
    }

    std::int32_t depth = 0;
    std::uint32_t maximumDepth = 0;
    bool ended = false;
    for (std::uint32_t local = 0; local < limit; ++local) {
        ExpressionInstruction instruction{};
        if (!loadRecord(context.expressions, offset + local, instruction)) {
            diagnostics.error("NVP103", "Expression instruction cannot be decoded.", std::string(path));
            return false;
        }
        const std::uint16_t opcode = instruction.opcode;
        if (opcode == 255U) {
            if (depth != 1) {
                diagnostics.error("NVP104", "Expression terminates with an invalid stack depth.", std::string(path));
            }
            if (!scanToEnd && local + 1U != limit) {
                diagnostics.error("NVP105", "Expression contains instructions after its end opcode.", std::string(path));
            }
            ended = true;
            break;
        }

        if (opcode == 0U) {
            if (!std::isfinite(instruction.immediate)) {
                diagnostics.error("NVP106", "Expression literal is non-finite.", std::string(path));
            }
            ++depth;
        } else if (opcode == 1U) {
            if (instruction.operand >= context.contract.speciesCount) {
                diagnostics.error("NVP107", "Expression species operand is out of range.", std::string(path));
            }
            ++depth;
        } else if (opcode == 2U) {
            const bool timeReference = (instruction.flags & 1U) != 0U;
            if (!timeReference && instruction.operand >= context.contract.parameterCount) {
                diagnostics.error("NVP108", "Expression parameter operand is out of range.", std::string(path));
            }
            ++depth;
        } else if (opcode == 3U) {
            if (depth < 1) diagnostics.error("NVP109", "Unary expression stack underflow.", std::string(path));
        } else if (opcode >= 4U && opcode <= 17U) {
            if (depth < 2) {
                diagnostics.error("NVP110", "Binary expression stack underflow.", std::string(path));
            } else {
                --depth;
            }
        } else if (opcode == 18U) {
            if (depth < 3) {
                diagnostics.error("NVP111", "Clamp expression stack underflow.", std::string(path));
            } else {
                depth -= 2;
            }
        } else if (opcode >= 19U && opcode <= 22U) {
            if (depth < 1) diagnostics.error("NVP112", "Temporal expression stack underflow.", std::string(path));
            if (instruction.auxiliary >= context.contract.temporalStateCount) {
                diagnostics.error("NVP113", "Temporal-state operand is out of range.", std::string(path));
            }
            if (!std::isfinite(instruction.immediate) || instruction.immediate < 0.0F) {
                diagnostics.error("NVP114", "Temporal duration is invalid.", std::string(path));
            }
        } else {
            diagnostics.error("NVP115", "Expression contains an unknown opcode.", std::string(path));
            return false;
        }

        if (depth < 0) depth = 0;
        maximumDepth = std::max(maximumDepth, static_cast<std::uint32_t>(depth));
        if (maximumDepth > context.contract.maximumExpressionStack || maximumDepth > 64U) {
            diagnostics.error("NVP116", "Expression exceeds the declared or executable stack depth.", std::string(path));
            return false;
        }
    }

    if (!ended) {
        diagnostics.error("NVP117", "Expression does not terminate within its declared execution budget.", std::string(path));
        return false;
    }
    return true;
}

std::uint32_t requiredParameterCount(std::uint32_t rateLaw) {
    switch (rateLaw) {
        case 0U: return 1U;
        case 1U: return 1U;
        case 2U: return 3U;
        case 3U: return 3U;
        case 4U: return 2U;
        case 5U: return 2U;
        case 6U: return 1U;
        case 7U: return 2U;
        case 8U: return 1U;
        case 255U: return 0U;
        default: return std::numeric_limits<std::uint32_t>::max();
    }
}

bool copyReport(std::string_view report, NVivoByteBuffer* output) {
    if (output == nullptr) return true;
    output->data = nullptr;
    output->size = 0;
    if (report.empty()) return true;
    void* allocation = std::malloc(report.size());
    if (allocation == nullptr) return false;
    std::memcpy(allocation, report.data(), report.size());
    output->data = static_cast<std::uint8_t*>(allocation);
    output->size = report.size();
    return true;
}

NVivoStatus validateSemantic(std::span<const std::byte> pack,
                             Diagnostics& diagnostics) {
    SemanticContext context;
    context.pack = pack;
    context.inspection = inspectProgramPack(pack, true);
    diagnostics.append(context.inspection.diagnostics);
    if (!context.inspection.valid) return NVIVO_STATUS_INVALID_PACK;

    context.strings = section(pack, context.inspection, PackSectionType::strings);
    context.species = section(pack, context.inspection, PackSectionType::species);
    context.parameters = section(pack, context.inspection, PackSectionType::parameters);
    context.reactionParameterIndices = section(pack, context.inspection, PackSectionType::reactionParameterIndices);
    context.stoichiometry = section(pack, context.inspection, PackSectionType::stoichiometry);
    context.reactions = section(pack, context.inspection, PackSectionType::reactions);
    context.expressions = section(pack, context.inspection, PackSectionType::expressions);
    context.actions = section(pack, context.inspection, PackSectionType::actions);
    context.rules = section(pack, context.inspection, PackSectionType::rules);
    context.monitors = section(pack, context.inspection, PackSectionType::monitors);
    context.cohorts = section(pack, context.inspection, PackSectionType::cohorts);
    context.incidenceOffsets = section(pack, context.inspection, PackSectionType::speciesIncidenceOffsets);
    context.incidence = section(pack, context.inspection, PackSectionType::speciesIncidence);
    context.runtimeContract = section(pack, context.inspection, PackSectionType::runtimeContract);

    if (!loadRecord(context.runtimeContract, 0U, context.contract)) {
        diagnostics.error("NVP118", "Runtime contract cannot be decoded.");
        return NVIVO_STATUS_INVALID_PACK;
    }

    const auto countMatches = [&](const SectionView& view,
                                  std::uint32_t expected,
                                  std::string_view name) {
        if (view.descriptor == nullptr || view.descriptor->count != expected) {
            diagnostics.error(
                "NVP119",
                std::string(name) + " count does not match the runtime contract."
            );
        }
    };
    countMatches(context.species, context.contract.speciesCount, "Species");
    countMatches(context.parameters, context.contract.parameterCount, "Parameter");
    countMatches(context.reactions, context.contract.reactionCount, "Reaction");
    countMatches(context.rules, context.contract.ruleCount, "Rule");
    countMatches(context.monitors, context.contract.monitorCount, "Monitor");
    countMatches(context.cohorts, context.contract.cohortCount, "Cohort");
    if (context.incidenceOffsets.descriptor == nullptr ||
        context.incidenceOffsets.descriptor->count != context.contract.speciesCount + 1U) {
        diagnostics.error("NVP120", "Species-incidence offset count must equal speciesCount + 1.");
    }
    if (context.contract.maximumExpressionStack > 64U) {
        diagnostics.error("NVP121", "Runtime contract exceeds the Metal expression-stack limit.");
    }
    if (context.contract.authoritativeScalarBytes != sizeof(float)) {
        diagnostics.error("NVP122", "Runtime contract scalar width is not executable by Metal ABI v1.");
    }
    if (context.contract.randomStreamVersion != 1U) {
        diagnostics.error("NVP123", "ProgramPack requires an unsupported random-stream version.");
    }
    if (context.strings.bytes.empty() || context.strings.bytes.front() != std::byte{0}) {
        diagnostics.error("NVP124", "String table must begin with the canonical empty string.");
    }
    if (context.contract.reserved[0] != 0U) {
        if (context.contract.reserved[0] > std::numeric_limits<std::uint32_t>::max()) {
            diagnostics.error("NVP125", "Manifest string offset exceeds the ProgramPack string ABI.");
        } else {
            requireString(
                context.strings,
                static_cast<std::uint32_t>(context.contract.reserved[0]),
                "runtimeContract.manifest",
                diagnostics
            );
        }
    }

    for (std::uint32_t index = 0; index < context.contract.speciesCount; ++index) {
        SpeciesRecord record{};
        if (!loadRecord(context.species, index, record)) {
            diagnostics.error("NVP126", "Species record cannot be decoded.", "species[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "species[" + std::to_string(index) + "]";
        requireString(context.strings, record.nameOffset, path + ".name", diagnostics);
        requireString(context.strings, record.compartmentOffset, path + ".compartment", diagnostics);
        requireString(context.strings, record.unitOffset, path + ".unit", diagnostics);
        if (!finiteOrdered(record.minimum, record.initialValue, record.maximum)) {
            diagnostics.error("NVP127", "Species initial value and bounds are invalid.", path);
        }
    }

    for (std::uint32_t index = 0; index < context.contract.parameterCount; ++index) {
        ParameterRecord record{};
        if (!loadRecord(context.parameters, index, record)) {
            diagnostics.error("NVP128", "Parameter record cannot be decoded.", "parameters[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "parameters[" + std::to_string(index) + "]";
        requireString(context.strings, record.nameOffset, path + ".name", diagnostics);
        requireString(context.strings, record.unitOffset, path + ".unit", diagnostics);
        requireString(context.strings, record.evidenceSourceOffset, path + ".evidenceSource", diagnostics);
        requireString(context.strings, record.evidenceDetailOffset, path + ".evidenceDetail", diagnostics);
        if (!finiteOrdered(record.minimum, record.value, record.maximum)) {
            diagnostics.error("NVP129", "Parameter value and bounds are invalid.", path);
        }
    }

    std::unordered_map<std::uint64_t, std::int64_t> expectedIncidence;
    expectedIncidence.reserve(
        context.stoichiometry.descriptor == nullptr ? 0U : context.stoichiometry.descriptor->count
    );
    std::vector<bool> reactionCovered(context.contract.reactionCount, false);

    for (std::uint32_t index = 0; index < context.contract.reactionCount; ++index) {
        ReactionRecord reaction{};
        if (!loadRecord(context.reactions, index, reaction)) {
            diagnostics.error("NVP130", "Reaction record cannot be decoded.", "reactions[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "reactions[" + std::to_string(index) + "]";
        requireString(context.strings, reaction.nameOffset, path + ".name", diagnostics);
        requireString(context.strings, reaction.compartmentOffset, path + ".compartment", diagnostics);

        const auto stoichiometryCount = context.stoichiometry.descriptor == nullptr
            ? 0U
            : context.stoichiometry.descriptor->count;
        if (!checkedRange(reaction.reactantOffset, reaction.reactantCount, stoichiometryCount) ||
            !checkedRange(reaction.productOffset, reaction.productCount, stoichiometryCount)) {
            diagnostics.error("NVP131", "Reaction stoichiometry range is invalid.", path);
        }
        const auto parameterIndexCount = context.reactionParameterIndices.descriptor == nullptr
            ? 0U
            : context.reactionParameterIndices.descriptor->count;
        if (!checkedRange(reaction.parameterOffset, reaction.parameterCount, parameterIndexCount)) {
            diagnostics.error("NVP132", "Reaction parameter-index range is invalid.", path);
        } else {
            for (std::uint32_t local = 0; local < reaction.parameterCount; ++local) {
                std::uint32_t parameterIndex = 0;
                if (!loadRecord(context.reactionParameterIndices, reaction.parameterOffset + local, parameterIndex) ||
                    parameterIndex >= context.contract.parameterCount) {
                    diagnostics.error("NVP133", "Reaction parameter index is invalid.", path);
                }
            }
        }

        const auto minimumParameters = requiredParameterCount(reaction.rateLaw);
        if (minimumParameters == std::numeric_limits<std::uint32_t>::max()) {
            diagnostics.error("NVP134", "Reaction uses an unknown rate law.", path);
        } else if (reaction.rateLaw != 255U && reaction.parameterCount < minimumParameters) {
            diagnostics.error("NVP135", "Reaction rate law has insufficient parameters.", path);
        }
        if (reaction.rateLaw == 255U) {
            validateExpression(
                context,
                reaction.expressionOffset,
                reaction.expressionCount,
                false,
                path + ".rateExpression",
                diagnostics
            );
        } else if (reaction.expressionCount != 0U) {
            diagnostics.error("NVP136", "Non-custom reaction unexpectedly contains rate bytecode.", path);
        }
        if ((reaction.flags & (1U << 1U)) != 0U) {
            validateExpression(
                context,
                reaction.reserved,
                0U,
                true,
                path + ".gate",
                diagnostics
            );
        }
        if (!std::isfinite(reaction.delaySeconds) || reaction.delaySeconds < 0.0F ||
            !std::isfinite(reaction.characteristicRate) || reaction.characteristicRate < 0.0F) {
            diagnostics.error("NVP137", "Reaction delay or characteristic rate is invalid.", path);
        }
        if (reaction.cohortIndex >= context.contract.cohortCount && context.contract.reactionCount != 0U) {
            diagnostics.error("NVP138", "Reaction cohort index is out of range.", path);
        }

        const auto accumulateTerms = [&](std::uint32_t offset,
                                         std::uint32_t count,
                                         std::int32_t sign,
                                         std::string_view role) {
            for (std::uint32_t local = 0; local < count; ++local) {
                StoichiometryRecord term{};
                if (!loadRecord(context.stoichiometry, offset + local, term)) {
                    diagnostics.error("NVP139", "Stoichiometry record cannot be decoded.", path);
                    continue;
                }
                if (term.speciesIndex >= context.contract.speciesCount || term.coefficient <= 0) {
                    diagnostics.error("NVP140", "Stoichiometry species or coefficient is invalid.", path + "." + std::string(role));
                    continue;
                }
                expectedIncidence[pairKey(term.speciesIndex, index)] +=
                    static_cast<std::int64_t>(sign) * term.coefficient;
            }
        };
        accumulateTerms(reaction.reactantOffset, reaction.reactantCount, -1, "reactants");
        accumulateTerms(reaction.productOffset, reaction.productCount, 1, "products");
    }

    for (std::uint32_t index = 0; index < context.contract.cohortCount; ++index) {
        CohortRecord cohort{};
        if (!loadRecord(context.cohorts, index, cohort)) {
            diagnostics.error("NVP141", "Cohort record cannot be decoded.", "cohorts[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "cohorts[" + std::to_string(index) + "]";
        if (!checkedRange(cohort.reactionOffset, cohort.reactionCount, context.contract.reactionCount)) {
            diagnostics.error("NVP142", "Cohort reaction range is invalid.", path);
            continue;
        }
        if (!std::isfinite(cohort.maximumStableStep) || cohort.maximumStableStep <= 0.0F ||
            !std::isfinite(cohort.stiffnessEstimate) || cohort.stiffnessEstimate < 0.0F ||
            cohort.preferredThreads == 0U || cohort.preferredThreads > 1'024U) {
            diagnostics.error("NVP143", "Cohort execution metadata is invalid.", path);
        }
        for (std::uint32_t local = 0; local < cohort.reactionCount; ++local) {
            const std::uint32_t reactionIndex = cohort.reactionOffset + local;
            if (reactionCovered[reactionIndex]) {
                diagnostics.error("NVP144", "Reaction appears in more than one cohort.", path);
            }
            reactionCovered[reactionIndex] = true;
            ReactionRecord reaction{};
            if (loadRecord(context.reactions, reactionIndex, reaction) && reaction.cohortIndex != index) {
                diagnostics.error("NVP145", "Reaction and cohort indices disagree.", path);
            }
        }
    }
    for (std::uint32_t index = 0; index < reactionCovered.size(); ++index) {
        if (!reactionCovered[index]) {
            diagnostics.error("NVP146", "Reaction is not covered by any execution cohort.", "reactions[" + std::to_string(index) + "]");
        }
    }

    for (std::uint32_t index = 0; index < context.contract.ruleCount; ++index) {
        RuleRecord rule{};
        if (!loadRecord(context.rules, index, rule)) {
            diagnostics.error("NVP147", "Rule record cannot be decoded.", "rules[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "rules[" + std::to_string(index) + "]";
        requireString(context.strings, rule.nameOffset, path + ".name", diagnostics);
        validateExpression(context, rule.conditionOffset, rule.conditionCount, false, path + ".condition", diagnostics);
        const auto actionCount = context.actions.descriptor == nullptr ? 0U : context.actions.descriptor->count;
        if (!checkedRange(rule.actionOffset, rule.actionCount, actionCount)) {
            diagnostics.error("NVP148", "Rule action range is invalid.", path);
        }
        if (!std::isfinite(rule.refractorySeconds) || rule.refractorySeconds < 0.0F) {
            diagnostics.error("NVP149", "Rule refractory duration is invalid.", path);
        }
    }

    const std::uint32_t actionCount = context.actions.descriptor == nullptr ? 0U : context.actions.descriptor->count;
    for (std::uint32_t index = 0; index < actionCount; ++index) {
        ActionRecord action{};
        if (!loadRecord(context.actions, index, action)) {
            diagnostics.error("NVP150", "Action record cannot be decoded.", "actions[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "actions[" + std::to_string(index) + "]";
        if (action.kind > 11U) diagnostics.error("NVP151", "Action kind is unknown.", path);
        const bool stringTarget = (action.flags & 1U) != 0U;
        if (stringTarget) {
            requireString(context.strings, action.targetIndex, path + ".target", diagnostics);
        } else if (action.targetIndex >= context.contract.speciesCount && action.kind < 10U) {
            diagnostics.error("NVP152", "Action species target is out of range.", path);
        }
        if (action.expressionCount > 0U) {
            validateExpression(context, action.expressionOffset, action.expressionCount, false, path + ".value", diagnostics);
        }
        requireString(context.strings, action.unitOffset, path + ".unit", diagnostics);
        if (!std::isfinite(action.constantValue) || !std::isfinite(action.maximumRate) || action.maximumRate < 0.0F) {
            diagnostics.error("NVP153", "Action value or maximum rate is invalid.", path);
        }
    }

    for (std::uint32_t index = 0; index < context.contract.monitorCount; ++index) {
        MonitorRecord monitor{};
        if (!loadRecord(context.monitors, index, monitor)) {
            diagnostics.error("NVP154", "Monitor record cannot be decoded.", "monitors[" + std::to_string(index) + "]");
            continue;
        }
        const std::string path = "monitors[" + std::to_string(index) + "]";
        requireString(context.strings, monitor.nameOffset, path + ".name", diagnostics);
        requireString(context.strings, monitor.messageOffset, path + ".message", diagnostics);
        validateExpression(context, monitor.expressionOffset, monitor.expressionCount, false, path + ".expression", diagnostics);
        if (monitor.response > 5U) diagnostics.error("NVP155", "Monitor response is unknown.", path);
    }

    std::unordered_map<std::uint64_t, std::int64_t> actualIncidence;
    actualIncidence.reserve(context.incidence.descriptor == nullptr ? 0U : context.incidence.descriptor->count);
    std::uint32_t previousOffset = 0U;
    for (std::uint32_t speciesIndex = 0; speciesIndex <= context.contract.speciesCount; ++speciesIndex) {
        std::uint32_t offset = 0U;
        if (!loadRecord(context.incidenceOffsets, speciesIndex, offset)) {
            diagnostics.error("NVP156", "Species-incidence offset cannot be decoded.");
            break;
        }
        const auto incidenceCount = context.incidence.descriptor == nullptr ? 0U : context.incidence.descriptor->count;
        if (offset < previousOffset || offset > incidenceCount) {
            diagnostics.error("NVP157", "Species-incidence offsets are not monotonic or in bounds.");
        }
        if (speciesIndex > 0U) {
            for (std::uint32_t item = previousOffset; item < offset; ++item) {
                IncidenceRecord incidence{};
                if (!loadRecord(context.incidence, item, incidence)) {
                    diagnostics.error("NVP158", "Incidence record cannot be decoded.");
                    continue;
                }
                if (incidence.reactionIndex >= context.contract.reactionCount || incidence.netCoefficient == 0) {
                    diagnostics.error("NVP159", "Incidence record contains an invalid reaction or coefficient.");
                    continue;
                }
                actualIncidence[pairKey(speciesIndex - 1U, incidence.reactionIndex)] += incidence.netCoefficient;
            }
        }
        previousOffset = offset;
    }
    if (context.incidence.descriptor != nullptr && previousOffset != context.incidence.descriptor->count) {
        diagnostics.error("NVP160", "Final species-incidence offset does not equal incidence count.");
    }

    for (auto iterator = expectedIncidence.begin(); iterator != expectedIncidence.end();) {
        if (iterator->second == 0) iterator = expectedIncidence.erase(iterator);
        else ++iterator;
    }
    for (auto iterator = actualIncidence.begin(); iterator != actualIncidence.end();) {
        if (iterator->second == 0) iterator = actualIncidence.erase(iterator);
        else ++iterator;
    }
    if (expectedIncidence != actualIncidence) {
        diagnostics.error("NVP161", "Sparse species-incidence graph does not match reaction stoichiometry.");
    }

    return diagnostics.hasErrors() ? NVIVO_STATUS_INVALID_PACK : NVIVO_STATUS_OK;
}

} // namespace

extern "C" NVivoStatus nvivo_validate_program_pack_semantics(
    const std::uint8_t* pack,
    std::size_t packSize,
    NVivoByteBuffer* diagnosticsJson) {
    if (diagnosticsJson != nullptr) {
        diagnosticsJson->data = nullptr;
        diagnosticsJson->size = 0;
    }
    if (pack == nullptr && packSize != 0U) return NVIVO_STATUS_INVALID_ARGUMENT;

    try {
        Diagnostics diagnostics;
        const auto bytes = std::span<const std::byte>(
            reinterpret_cast<const std::byte*>(pack),
            packSize
        );
        const NVivoStatus status = validateSemantic(bytes, diagnostics);
        std::ostringstream report;
        report << "{\"valid\":" << (status == NVIVO_STATUS_OK ? "true" : "false")
               << ",\"diagnostics\":" << diagnostics.toJson() << '}';
        const std::string text = report.str();
        if (!copyReport(text, diagnosticsJson)) return NVIVO_STATUS_OUT_OF_MEMORY;
        return status;
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (...) {
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}
