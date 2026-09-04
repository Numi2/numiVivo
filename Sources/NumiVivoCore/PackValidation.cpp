#include "PackValidation.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <unordered_map>

namespace nvivo {
namespace {
constexpr std::uint32_t invalid = UINT32_MAX;
constexpr std::uint32_t maximumRecords = 1'048'576;

struct InvalidPack final : std::runtime_error {
    std::string path;
    InvalidPack(std::string message, std::string subject)
        : std::runtime_error(std::move(message)), path(std::move(subject)) {}
};
void require(bool condition, std::string_view message, std::string_view path) {
    if (!condition) throw InvalidPack(std::string(message), std::string(path));
}

// Reject overlong sequences, surrogate scalars and codepoints above U+10FFFF.
bool validUTF8(std::string_view text) {
    std::size_t i = 0;
    while (i < text.size()) {
        const auto lead = static_cast<unsigned char>(text[i++]);
        if (lead < 0x80) continue;
        unsigned remaining = 0;
        std::uint32_t scalar = 0, minimum = 0;
        if (lead >= 0xc2 && lead <= 0xdf) { remaining = 1; scalar = lead & 0x1f; minimum = 0x80; }
        else if (lead >= 0xe0 && lead <= 0xef) { remaining = 2; scalar = lead & 0x0f; minimum = 0x800; }
        else if (lead >= 0xf0 && lead <= 0xf4) { remaining = 3; scalar = lead & 7; minimum = 0x10000; }
        else return false;
        if (remaining > text.size() - i) return false;
        while (remaining-- > 0) {
            const auto next = static_cast<unsigned char>(text[i++]);
            if ((next & 0xc0) != 0x80) return false;
            scalar = (scalar << 6) | (next & 0x3f);
        }
        if (scalar < minimum || scalar > 0x10ffff || (scalar >= 0xd800 && scalar <= 0xdfff)) return false;
    }
    return true;
}

class Tables {
    std::span<const std::byte> bytes_;
    std::map<PackSectionType, PackSectionDescriptor> sections_;
    std::string_view strings_;
    std::vector<std::uint32_t> stringStarts_;
public:
    Tables(std::span<const std::byte> bytes, const std::vector<PackSectionDescriptor>& sections) : bytes_(bytes) {
        for (const auto& section : sections) {
            if (section.type >= 1 && section.type <= 14)
                sections_.emplace(static_cast<PackSectionType>(section.type), section);
        }
        for (std::uint32_t kind = 1; kind <= 14; ++kind)
            require(sections_.contains(static_cast<PackSectionType>(kind)), "Required v1 table is absent.", "sections");
        const auto& strings = descriptor(PackSectionType::strings);
        require(strings.size > 0 && strings.size <= 16U * 1024U * 1024U, "String table exceeds its bound.", "strings");
        strings_ = {reinterpret_cast<const char*>(bytes.data() + strings.offset), static_cast<std::size_t>(strings.size)};
        require(strings_[0] == '\0' && strings_.back() == '\0' && validUTF8(strings_), "Invalid UTF-8 or unterminated string table.", "strings");
        stringStarts_.push_back(0);
        for (std::size_t i = 1; i < strings_.size(); ++i) {
            if (strings_[i - 1] == '\0') stringStarts_.push_back(static_cast<std::uint32_t>(i));
        }
    }
    const PackSectionDescriptor& descriptor(PackSectionType type) const { return sections_.at(type); }
    std::uint32_t count(PackSectionType type) const { return descriptor(type).count; }
    template<class T> T read(PackSectionType type, std::uint32_t index) const {
        static_assert(std::is_trivially_copyable_v<T>);
        const auto& section = descriptor(type);
        require(section.stride == sizeof(T) && index < section.count, "Table read is out of bounds.", "tables");
        T value{};
        const auto offset = section.offset + std::uint64_t(index) * section.stride;
        require(offset <= bytes_.size() && sizeof(T) <= bytes_.size() - offset, "Payload read is out of bounds.", "tables");
        std::memcpy(&value, bytes_.data() + offset, sizeof(T));
        return value;
    }
    void range(PackSectionType type, std::uint32_t offset, std::uint32_t length, std::string_view path) const {
        const auto size = count(type);
        require(offset <= size && length <= size - offset, "Referenced table range is out of bounds.", path);
    }
    std::string_view string(std::uint32_t offset, std::string_view path) const {
        const auto found = std::lower_bound(stringStarts_.begin(), stringStarts_.end(), offset);
        require(found != stringStarts_.end() && *found == offset, "String reference is not a string boundary.", path);
        const auto next = std::next(found);
        const std::size_t end = next == stringStarts_.end() ? strings_.size() : *next;
        return strings_.substr(offset, end - offset - 1);
    }
};

struct ExpressionSlice { std::uint32_t count; std::uint32_t maximumStack; };
using Expressions = std::map<std::uint32_t, ExpressionSlice>;

Expressions validateExpressions(const Tables& t, const RuntimeContractRecord& contract) {
    Expressions result;
    std::set<std::uint32_t> temporalSlots;
    std::uint32_t start = 0, depth = 0, maximum = 0, globalMaximum = 0;
    const auto size = t.count(PackSectionType::expressions);
    for (std::uint32_t i = 0; i < size; ++i) {
        const auto instruction = t.read<ExpressionInstruction>(PackSectionType::expressions, i);
        const auto opcode = instruction.opcode;
        require(std::isfinite(instruction.immediate), "Expression immediate is nonfinite.", "expressions");
        require(instruction.flags == 0 || (opcode == 2 && instruction.flags == 1), "Unknown expression flags.", "expressions");
        if (opcode <= 2) {
            require(depth < 256, "Expression stack overflows the GPU VM.", "expressions");
            if (opcode == 1) require(instruction.operand < contract.speciesCount, "Invalid species operand.", "expressions");
            if (opcode == 2 && instruction.flags == 0)
                require(instruction.operand < contract.parameterCount, "Invalid parameter operand.", "expressions");
            if (opcode == 2 && instruction.flags == 1)
                require(instruction.operand == invalid, "Time reference must use the time sentinel.", "expressions");
            maximum = std::max(maximum, ++depth);
        } else if (opcode == 255) {
            require(depth == 1, "Expression must terminate with exactly one value.", "expressions");
            result.emplace(start, ExpressionSlice{i - start + 1, maximum});
            globalMaximum = std::max(globalMaximum, maximum);
            start = i + 1; depth = 0; maximum = 0;
        } else if (opcode == 3 || (opcode >= 19 && opcode <= 22)) {
            require(depth >= 1, "Unary expression stack underflow.", "expressions");
            if (opcode >= 19) {
                require(instruction.auxiliary < contract.temporalStateCount, "Temporal slot is out of bounds.", "expressions");
                require(temporalSlots.insert(instruction.auxiliary).second, "Temporal instructions alias mutable state.", "expressions");
                if (opcode <= 20) require(instruction.immediate > 0, "Temporal duration must be positive.", "expressions");
            }
        } else if (opcode >= 4 && opcode <= 17) {
            require(depth >= 2, "Binary expression stack underflow.", "expressions");
            --depth;
        } else if (opcode == 18) {
            require(depth >= 3, "Clamp expression stack underflow.", "expressions");
            depth -= 2;
        } else require(false, "Unknown expression opcode.", "expressions");
    }
    require(start == size, "Expression table has an unterminated tail.", "expressions");
    require(globalMaximum == contract.maximumExpressionStack, "Declared VM stack differs from bytecode.", "runtimeContract");
    require(temporalSlots.size() == contract.temporalStateCount, "Temporal-state count differs from bytecode ownership.", "runtimeContract");
    return result;
}
void expressionReference(const Expressions& expressions, std::uint32_t start, std::uint32_t count,
                         std::string_view path, bool gate = false) {
    const auto found = expressions.find(start);
    require(found != expressions.end(), "Expression reference is not a program boundary.", path);
    require(gate ? (found->second.count <= 4096) : (count == found->second.count),
            "Expression length disagrees with its program boundary or gate budget.", path);
}

void validateTables(const Tables& t, const PackHeader& header) {
    using S = PackSectionType;
    require(t.count(S::runtimeContract) == 1, "Runtime contract must have one record.", "runtimeContract");
    const auto c = t.read<RuntimeContractRecord>(S::runtimeContract, 0);
    const ResourceLimits limits;
    require(c.speciesCount <= limits.maximumSpecies && c.parameterCount <= limits.maximumParameters &&
            c.reactionCount <= limits.maximumReactions && c.ruleCount <= limits.maximumRules &&
            c.monitorCount <= limits.maximumConstraints && c.temporalStateCount <= limits.maximumTemporalStates,
            "Pack exceeds supported compiled resource capacities.", "runtimeContract");
    for (std::uint32_t type = 2; type <= 14; ++type)
        require(t.count(static_cast<S>(type)) <= maximumRecords, "Table exceeds its bounded record count.", "sections");
    require(c.speciesCount == t.count(S::species) && c.parameterCount == t.count(S::parameters) &&
            c.reactionCount == t.count(S::reactions) && c.ruleCount == t.count(S::rules) &&
            c.monitorCount == t.count(S::monitors) && c.cohortCount == t.count(S::cohorts),
            "Runtime counts disagree with descriptor counts.", "runtimeContract");
    require(c.authoritativeScalarBytes == 4 && c.randomStreamVersion == 1 &&
            c.maximumExpressionStack <= 256 && c.reserved0 == 0 &&
            c.featureFlags == header.flags && (header.flags & ~0x7fu) == 0,
            "Unsupported runtime ABI, precision, random stream or feature flags.", "runtimeContract");
    require(c.reserved[1] == 0 && c.reserved[2] == 0 && c.reserved[3] == 0,
            "Unknown required runtime contract extension.", "runtimeContract");

    std::set<std::string_view> symbols;
    for (std::uint32_t i = 0; i < c.speciesCount; ++i) {
        const auto value = t.read<SpeciesRecord>(S::species, i);
        auto name = t.string(value.nameOffset, "species.name");
        require(isValidIdentifier(name) && symbols.insert(name).second, "Invalid or duplicate species identifier.", "species");
        t.string(value.compartmentOffset, "species.compartment");
        t.string(value.unitOffset, "species.unit");
        require((value.flags & ~0x3fu) == 0 && value.reserved == 0, "Unknown species flags or extension.", "species");
        require(std::isfinite(value.initialValue) && std::isfinite(value.minimum) && std::isfinite(value.maximum) &&
                value.minimum <= value.initialValue && value.initialValue <= value.maximum,
                "Invalid species value or bounds.", "species");
        if ((value.flags & speciesCountValued) != 0)
            require(value.initialValue >= 0 && value.initialValue <= 16'777'216.0f && std::trunc(value.initialValue) == value.initialValue,
                    "FP32 count state is not an exactly representable nonnegative integer.", "species");
    }
    for (std::uint32_t i = 0; i < c.parameterCount; ++i) {
        const auto value = t.read<ParameterRecord>(S::parameters, i);
        auto name = t.string(value.nameOffset, "parameters.name");
        require(isValidIdentifier(name) && symbols.insert(name).second, "Invalid or duplicate parameter identifier.", "parameters");
        t.string(value.unitOffset, "parameters.unit"); t.string(value.evidenceSourceOffset, "parameters.evidence");
        require(value.evidenceClass <= 5 && (value.flags & ~3u) == 0 && value.reserved == 0,
                "Unknown parameter evidence, flags or extension.", "parameters");
        require(std::isfinite(value.value) && std::isfinite(value.minimum) && std::isfinite(value.maximum) &&
                value.minimum <= value.value && value.value <= value.maximum &&
                std::abs(value.value) <= std::numeric_limits<float>::max(), "Invalid parameter value or bounds.", "parameters");
        const auto rounded = static_cast<float>(value.value);
        require(value.value == 0 || rounded != 0, "Parameter underflows the FP32 runtime.", "parameters");
    }
    for (std::uint32_t i = 0; i < t.count(S::reactionParameterIndices); ++i)
        require(t.read<std::uint32_t>(S::reactionParameterIndices, i) < c.parameterCount, "Invalid reaction parameter index.", "reactionParameterIndices");
    for (std::uint32_t i = 0; i < t.count(S::stoichiometry); ++i) {
        const auto term = t.read<StoichiometryRecord>(S::stoichiometry, i);
        require(term.speciesIndex < c.speciesCount && term.coefficient > 0 && term.role <= 1,
                "Invalid stoichiometry species, coefficient or role.", "stoichiometry");
    }
    const auto expressions = validateExpressions(t, c);
    std::map<std::uint64_t, std::int64_t> expectedIncidence;
    std::vector<bool> reactionCovered(c.reactionCount, false);
    std::set<std::string_view> reactionNames;
    for (std::uint32_t i = 0; i < c.reactionCount; ++i) {
        const auto r = t.read<ReactionRecord>(S::reactions, i);
        auto name = t.string(r.nameOffset, "reactions.name");
        require(isValidIdentifier(name) && reactionNames.insert(name).second, "Invalid or duplicate reaction identifier.", "reactions");
        t.string(r.compartmentOffset, "reactions.compartment");
        require(r.rateLaw <= 8 || r.rateLaw == 255, "Unknown reaction rate law.", "reactions");
        require((r.flags & ~0x1fu) == 0 && std::isfinite(r.delaySeconds) && r.delaySeconds >= 0 &&
                std::isfinite(r.characteristicRate) && r.characteristicRate >= 0 && r.cohortIndex < c.cohortCount,
                "Invalid reaction flags, rate, delay or cohort.", "reactions");
        require(((r.flags & reactionDelayed) != 0) == (r.delaySeconds > 0), "Delayed flag differs from delay duration.", "reactions");
        t.range(S::reactionParameterIndices, r.parameterOffset, r.parameterCount, "reactions.parameters");
        static constexpr std::uint32_t parameterCounts[] = {1, 1, 3, 3, 2, 2, 1, 2, 1};
        if (r.rateLaw <= 8) require(r.parameterCount == parameterCounts[r.rateLaw], "Wrong kinetic parameter arity.", "reactions");
        if (r.rateLaw == 255) require(r.expressionCount > 0, "Custom kinetics have no expression.", "reactions");
        if (r.expressionCount) expressionReference(expressions, r.expressionOffset, r.expressionCount, "reactions.expression");
        if ((r.flags & reactionHasGate) != 0) expressionReference(expressions, r.reserved, 0, "reactions.gate", true);
        else require(r.reserved == invalid, "Ungated reaction has a gate offset.", "reactions");
        require(r.reactantCount > 0 || r.productCount > 0, "Reaction has no stoichiometry.", "reactions");
        auto accumulate = [&](std::uint32_t offset, std::uint32_t count, std::uint16_t role, int sign) {
            t.range(S::stoichiometry, offset, count, "reactions.stoichiometry");
            for (std::uint32_t j = 0; j < count; ++j) {
                const auto term = t.read<StoichiometryRecord>(S::stoichiometry, offset + j);
                require(term.role == role, "Reactant/product role mismatch.", "reactions.stoichiometry");
                expectedIncidence[(std::uint64_t(term.speciesIndex) << 32) | i] += sign * std::int64_t(term.coefficient);
            }
        };
        accumulate(r.reactantOffset, r.reactantCount, 0, -1);
        accumulate(r.productOffset, r.productCount, 1, 1);
    }
    for (std::uint32_t i = 0; i < c.cohortCount; ++i) {
        const auto cohort = t.read<CohortRecord>(S::cohorts, i);
        t.range(S::reactions, cohort.reactionOffset, cohort.reactionCount, "cohorts");
        require(cohort.reactionCount > 0 && std::isfinite(cohort.maximumStableStep) && cohort.maximumStableStep > 0 &&
                std::isfinite(cohort.stiffnessEstimate) && cohort.stiffnessEstimate >= 0 && cohort.reserved == 0,
                "Invalid cohort descriptor.", "cohorts");
        for (std::uint32_t j = 0; j < cohort.reactionCount; ++j) {
            const auto index = cohort.reactionOffset + j;
            const auto reaction = t.read<ReactionRecord>(S::reactions, index);
            require(!reactionCovered[index] && reaction.cohortIndex == i && reaction.rateLaw == cohort.rateLaw,
                    "Cohort coverage or kinetic class is inconsistent.", "cohorts");
            reactionCovered[index] = true;
        }
    }
    require(std::all_of(reactionCovered.begin(), reactionCovered.end(), [](bool v) { return v; }), "Uncovered reaction.", "cohorts");

    require(std::uint64_t(t.count(S::speciesIncidenceOffsets)) == std::uint64_t(c.speciesCount) + 1 &&
            t.read<std::uint32_t>(S::speciesIncidenceOffsets, 0) == 0, "Invalid CSR offset shape.", "incidence");
    std::uint32_t end = 0;
    for (std::uint32_t s = 0; s < c.speciesCount; ++s) {
        const auto begin = end;
        end = t.read<std::uint32_t>(S::speciesIncidenceOffsets, s + 1);
        require(end >= begin && end <= t.count(S::speciesIncidence), "CSR offsets are not monotone or bounded.", "incidence");
        std::uint32_t previous = 0; bool first = true;
        for (auto i = begin; i < end; ++i) {
            const auto edge = t.read<IncidenceRecord>(S::speciesIncidence, i);
            require(edge.reactionIndex < c.reactionCount && edge.netCoefficient != 0 &&
                    (first || edge.reactionIndex > previous) && edge.reserved16 == 0 && edge.reserved0 == 0 && edge.reserved1 == 0,
                    "Invalid, duplicate or unsorted CSR incidence.", "incidence");
            const auto key = (std::uint64_t(s) << 32) | edge.reactionIndex;
            const auto expected = expectedIncidence.find(key);
            require(expected != expectedIncidence.end() && expected->second == edge.netCoefficient,
                    "CSR incidence does not equal reaction stoichiometry.", "incidence");
            expectedIncidence.erase(expected);
            first = false; previous = edge.reactionIndex;
        }
    }
    require(end == t.count(S::speciesIncidence), "Unreachable CSR tail.", "incidence");
    for (const auto& [key, coefficient] : expectedIncidence) {
        (void)key;
        require(coefficient == 0, "Reaction stoichiometry is missing from CSR.", "incidence");
    }
    for (std::uint32_t i = 0; i < t.count(S::actions); ++i) {
        const auto a = t.read<ActionRecord>(S::actions, i);
        require(a.kind <= 11 && (a.flags & ~1u) == 0 && std::isfinite(a.constantValue) &&
                std::isfinite(a.maximumRate) && a.maximumRate >= 0, "Invalid action opcode, flags or numeric value.", "actions");
        t.string(a.unitOffset, "actions.unit");
        if (a.expressionCount) expressionReference(expressions, a.expressionOffset, a.expressionCount, "actions.expression");
        if (a.kind >= 10) require(a.flags == 0 && a.targetIndex == invalid, "Shutdown action has a writable target.", "actions");
        else if (a.flags == 1) {
            require(a.kind >= 7 && a.kind <= 9, "String action has an invalid opcode.", "actions");
            t.string(a.targetIndex, "actions.target");
        } else {
            require(a.targetIndex < c.speciesCount, "Action target is out of bounds.", "actions");
            require((t.read<SpeciesRecord>(S::species, a.targetIndex).flags & speciesExternallyOwned) == 0,
                    "Action writes externally owned state.", "actions");
        }
    }
    for (std::uint32_t i = 0; i < c.ruleCount; ++i) {
        const auto r = t.read<RuleRecord>(S::rules, i);
        t.string(r.nameOffset, "rules.name");
        require(std::isfinite(r.refractorySeconds) && r.refractorySeconds >= 0 && r.temporalStateOffset <= c.temporalStateCount,
                "Invalid rule timing or temporal offset.", "rules");
        t.range(S::actions, r.actionOffset, r.actionCount, "rules.actions");
        expressionReference(expressions, r.conditionOffset, r.conditionCount, "rules.condition");
    }
    for (std::uint32_t i = 0; i < c.monitorCount; ++i) {
        const auto m = t.read<MonitorRecord>(S::monitors, i);
        t.string(m.nameOffset, "monitors.name"); t.string(m.messageOffset, "monitors.message");
        require(m.severity <= 3 && m.response <= 5 && (m.flags & ~1u) == 0 && m.temporalStateOffset <= c.temporalStateCount,
                "Invalid monitor response, severity, flags or temporal offset.", "monitors");
        if (m.flags == 1) require(m.response == 4 || m.response == 5, "Termination monitor does not terminate.", "monitors");
        expressionReference(expressions, m.expressionOffset, m.expressionCount, "monitors.expression");
    }
    // v1 hashes exclude the header. Bind its scientific identity to the hashed
    // manifest and contract; integrity hashing is still not a digital signature.
    require(c.reserved[0] <= UINT32_MAX, "Manifest offset exceeds the string ABI.", "manifest");
    const auto manifest = t.string(static_cast<std::uint32_t>(c.reserved[0]), "manifest");
    auto parsed = json::parse(manifest);
    require(parsed.root.has_value() && parsed.root->isObject(), "Missing or invalid compilation manifest.", "manifest");
    const auto* program = parsed.root->get("program");
    const auto* fidelity = parsed.root->get("fidelity");
    const auto* fingerprint = program ? program->get("sourceFingerprint") : nullptr;
    require(fingerprint && fingerprint->isString() && fingerprint->asString() == hexFingerprint(header.sourceFingerprint) &&
            fidelity && fidelity->isString() && fidelity->asString() == fidelityName(static_cast<FidelityLevel>(header.fidelity)),
            "Header identity disagrees with the hashed manifest.", "manifest");
}
} // namespace

void validateProgramPackTables(std::span<const std::byte> bytes, const PackHeader& header,
                               const std::vector<PackSectionDescriptor>& sections, Diagnostics& diagnostics) {
    try {
        const Tables tables(bytes, sections);
        validateTables(tables, header);
    } catch (const InvalidPack& error) {
        diagnostics.error("NVK100", error.what(), error.path);
    }
}
} // namespace nvivo
