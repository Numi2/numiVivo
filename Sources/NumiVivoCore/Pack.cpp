#include "NumiVivoCore/Core.hpp"
#include "PackValidation.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>
#include <type_traits>

namespace nvivo {
namespace {
static_assert(std::endian::native == std::endian::little, "ProgramPack v1 requires a little-endian host");
static_assert(sizeof(PackHeader) == 128 && sizeof(PackSectionDescriptor) == 72);
constexpr std::array<char, 8> magic = {'N','V','I','V','O','P','K','\0'};
constexpr std::uint64_t maximumPackBytes = 256ULL * 1024ULL * 1024ULL;
struct SectionABI { std::uint32_t stride; std::uint32_t alignment; const char* name; };
constexpr SectionABI abi[] = {
    {0,1,"unknown"}, {1,1,"strings"}, {sizeof(SpeciesRecord),alignof(SpeciesRecord),"species"},
    {sizeof(ParameterRecord),alignof(ParameterRecord),"parameters"}, {4,4,"reactionParameterIndices"},
    {sizeof(StoichiometryRecord),alignof(StoichiometryRecord),"stoichiometry"},
    {sizeof(ReactionRecord),alignof(ReactionRecord),"reactions"},
    {sizeof(ExpressionInstruction),alignof(ExpressionInstruction),"expressions"},
    {sizeof(ActionRecord),alignof(ActionRecord),"actions"},
    {sizeof(RuleRecord),alignof(RuleRecord),"rules"},
    {sizeof(MonitorRecord),alignof(MonitorRecord),"monitors"},
    {sizeof(CohortRecord),alignof(CohortRecord),"cohorts"}, {4,4,"speciesIncidenceOffsets"},
    {sizeof(IncidenceRecord),alignof(IncidenceRecord),"speciesIncidence"},
    {sizeof(RuntimeContractRecord),alignof(RuntimeContractRecord),"runtimeContract"}
};
bool powerOfTwo(std::uint64_t value) { return value != 0 && (value & (value - 1)) == 0; }
const char* name(std::uint32_t type) { return type < std::size(abi) ? abi[type].name : "unknown"; }
template<class T> std::span<const std::byte> asBytes(const std::vector<T>& values) {
    static_assert(std::is_trivially_copyable_v<T>);
    return {reinterpret_cast<const std::byte*>(values.data()), values.size() * sizeof(T)};
}
struct Input { PackSectionType type; std::span<const std::byte> bytes; std::uint32_t count; };
}

PackBuildResult serializeProgramPack(const ProgramIR& ir) {
    PackBuildResult result;
    auto count = [&](std::size_t value) -> std::uint32_t {
        if (value > UINT32_MAX) { result.diagnostics.error("NVK001", "Table count exceeds UInt32."); return 0; }
        return static_cast<std::uint32_t>(value);
    };
    RuntimeContractRecord contract{};
    contract.speciesCount = count(ir.species.size());
    contract.parameterCount = count(ir.parameters.size());
    contract.reactionCount = count(ir.reactions.size());
    contract.ruleCount = count(ir.rules.size());
    contract.monitorCount = count(ir.monitors.size());
    contract.cohortCount = count(ir.cohorts.size());
    contract.temporalStateCount = ir.temporalStateCount;
    contract.maximumExpressionStack = ir.maximumExpressionStack;
    contract.featureFlags = ir.featureFlags;
    contract.authoritativeScalarBytes = 4;
    contract.randomStreamVersion = 1;
    if (ir.strings.size() > 1) {
        std::size_t start = ir.strings.size() - 2;
        while (start > 0 && ir.strings[start - 1] != '\0') --start;
        contract.reserved[0] = count(start);
    }
    const std::array<Input,14> inputs = {{
        {PackSectionType::strings, asBytes(ir.strings), count(ir.strings.size())},
        {PackSectionType::species, asBytes(ir.species), contract.speciesCount},
        {PackSectionType::parameters, asBytes(ir.parameters), contract.parameterCount},
        {PackSectionType::reactionParameterIndices, asBytes(ir.reactionParameterIndices), count(ir.reactionParameterIndices.size())},
        {PackSectionType::stoichiometry, asBytes(ir.stoichiometry), count(ir.stoichiometry.size())},
        {PackSectionType::reactions, asBytes(ir.reactions), contract.reactionCount},
        {PackSectionType::expressions, asBytes(ir.expressions), count(ir.expressions.size())},
        {PackSectionType::actions, asBytes(ir.actions), count(ir.actions.size())},
        {PackSectionType::rules, asBytes(ir.rules), contract.ruleCount},
        {PackSectionType::monitors, asBytes(ir.monitors), contract.monitorCount},
        {PackSectionType::cohorts, asBytes(ir.cohorts), contract.cohortCount},
        {PackSectionType::speciesIncidenceOffsets, asBytes(ir.speciesIncidenceOffsets), count(ir.speciesIncidenceOffsets.size())},
        {PackSectionType::speciesIncidence, asBytes(ir.speciesIncidence), count(ir.speciesIncidence.size())},
        {PackSectionType::runtimeContract, {reinterpret_cast<const std::byte*>(&contract), sizeof(contract)}, 1}
    }};
    if (result.diagnostics.hasErrors()) return result;
    const std::size_t headerBytes = sizeof(PackHeader) + inputs.size() * sizeof(PackSectionDescriptor);
    std::vector<PackSectionDescriptor> sections;
    sections.reserve(inputs.size());
    std::uint64_t total = headerBytes;
    for (const auto& input : inputs) {
        const auto type = static_cast<std::uint32_t>(input.type);
        const std::uint32_t alignment = type == 1 || type == 14 ? 16 : 256;
        const std::uint64_t aligned = (total + alignment - 1) & ~std::uint64_t(alignment - 1);
        if (aligned > maximumPackBytes || input.bytes.size() > maximumPackBytes - aligned) {
            result.diagnostics.error("NVK002", "ProgramPack exceeds the 256 MiB format limit.");
            return result;
        }
        PackSectionDescriptor section{};
        section.type = type;
        section.flags = type == 1 || type == 2 || type == 3 || type == 6 || type == 10 || type == 11 || type == 12 || type == 14 ? 1 : 0;
        section.offset = aligned;
        section.size = input.bytes.size();
        section.stride = abi[type].stride;
        section.count = input.count;
        section.alignment = alignment;
        section.fingerprint = sha256(input.bytes);
        sections.push_back(section);
        total = aligned + input.bytes.size();
    }
    result.bytes.resize(static_cast<std::size_t>(total), std::byte{0});
    for (std::size_t i = 0; i < inputs.size(); ++i) {
        if (!inputs[i].bytes.empty()) std::memcpy(result.bytes.data() + sections[i].offset, inputs[i].bytes.data(), inputs[i].bytes.size());
    }
    PackHeader header{};
    header.magic = magic;
    header.major = kProgramPackMajor;
    header.minor = kProgramPackMinor;
    header.headerBytes = static_cast<std::uint32_t>(headerBytes);
    header.compilerABI = kCompilerABIVersion;
    header.flags = ir.featureFlags;
    header.fidelity = static_cast<std::uint32_t>(ir.fidelity);
    header.sectionCount = static_cast<std::uint32_t>(sections.size());
    header.totalBytes = total;
    header.sourceFingerprint = ir.sourceFingerprint;
    std::memcpy(result.bytes.data(), &header, sizeof(header));
    std::memcpy(result.bytes.data() + sizeof(header), sections.data(), sections.size() * sizeof(PackSectionDescriptor));
    // Preserve v1 wire identity. Inspection binds the unhashed identity fields
    // to the hashed runtime contract and compilation manifest.
    header.contentFingerprint = sha256(std::span<const std::byte>(result.bytes).subspan(sizeof(header)));
    std::memcpy(result.bytes.data(), &header, sizeof(header));
    const auto inspection = inspectProgramPack(result.bytes, true);
    result.diagnostics.append(inspection.diagnostics);
    if (!inspection.valid) { result.bytes.clear(); return result; }
    result.fingerprint = header.contentFingerprint;
    return result;
}

PackInspection inspectProgramPack(std::span<const std::byte> bytes, bool verifySectionHashes) {
    PackInspection result;
    auto error = [&](std::string_view message, std::string_view path = "pack") {
        result.diagnostics.error("NVK010", std::string(message), std::string(path));
    };
    if (bytes.size() < sizeof(PackHeader) || bytes.size() > maximumPackBytes) {
        error("Input size is outside the bounded ProgramPack format."); return result;
    }
    std::memcpy(&result.header, bytes.data(), sizeof(PackHeader));
    const auto& h = result.header;
    if (h.magic != magic || h.major != kProgramPackMajor || h.minor != kProgramPackMinor || h.compilerABI != kCompilerABIVersion) {
        error("Unknown magic, format version or compiler ABI."); return result;
    }
    if (h.fidelity > 4 || (h.flags & ~0x7fu) != 0 ||
        !std::all_of(h.reserved.begin(), h.reserved.end(), [](auto v) { return v == 0; })) {
        error("Unsupported header fidelity, feature flags or extension."); return result;
    }
    if (h.sectionCount > 1024 || h.sectionCount < 14 ||
        std::uint64_t(h.headerBytes) != sizeof(PackHeader) + std::uint64_t(h.sectionCount) * sizeof(PackSectionDescriptor) ||
        h.headerBytes > bytes.size() || h.totalBytes != bytes.size()) {
        error("Invalid header/descriptor/total size."); return result;
    }
    result.sections.resize(h.sectionCount);
    std::memcpy(result.sections.data(), bytes.data() + sizeof(PackHeader), result.sections.size() * sizeof(PackSectionDescriptor));
    std::set<std::uint32_t> seen;
    std::vector<std::pair<std::uint64_t,std::uint64_t>> ranges;
    for (const auto& s : result.sections) {
        if (!seen.insert(s.type).second || (s.flags & ~1u) != 0 || s.reserved != 0) {
            error("Duplicate section type or unknown descriptor extension.", name(s.type)); continue;
        }
        if (s.type == 0 || s.type >= std::size(abi)) {
            if ((s.flags & 1) != 0) error("Unknown required section.");
            // Unknown optional sections are still structurally and hash checked.
        } else if (s.stride != abi[s.type].stride || s.alignment < abi[s.type].alignment) {
            error("Section stride or minimum alignment differs from the record ABI.", name(s.type));
        }
        if (!powerOfTwo(s.alignment) || s.alignment > 65536 || s.offset % s.alignment != 0 ||
            s.offset < h.headerBytes || s.offset > bytes.size() || s.size > bytes.size() - s.offset ||
            std::uint64_t(s.count) * s.stride != s.size) {
            error("Section has an invalid range, alignment, stride or count.", name(s.type)); continue;
        }
        if (s.size) ranges.emplace_back(s.offset, s.offset + s.size);
        if (verifySectionHashes && sha256(bytes.subspan(static_cast<std::size_t>(s.offset), static_cast<std::size_t>(s.size))) != s.fingerprint)
            error("Section fingerprint mismatch.", name(s.type));
    }
    std::sort(ranges.begin(), ranges.end());
    for (std::size_t i = 1; i < ranges.size(); ++i)
        if (ranges[i].first < ranges[i - 1].second) error("Payload sections overlap.");
    for (std::uint32_t type = 1; type <= 14; ++type)
        if (!seen.contains(type)) error("Required v1 table is missing.", name(type));
    if (sha256(bytes.subspan(sizeof(PackHeader))) != h.contentFingerprint) error("Content fingerprint mismatch.");
    if (result.diagnostics.hasErrors()) return result;
    validateProgramPackTables(bytes, h, result.sections, result.diagnostics);
    result.valid = !result.diagnostics.hasErrors();
    return result;
}

std::string PackInspection::toJson() const {
    std::ostringstream out;
    out << "{\"valid\":" << (valid ? "true" : "false")
        << ",\"header\":{\"major\":" << header.major << ",\"minor\":" << header.minor
        << ",\"compilerABI\":" << header.compilerABI << ",\"flags\":" << header.flags
        << ",\"fidelity\":" << header.fidelity << ",\"sectionCount\":" << header.sectionCount
        << ",\"totalBytes\":" << header.totalBytes << ",\"sourceFingerprint\":\"" << hexFingerprint(header.sourceFingerprint)
        << "\",\"contentFingerprint\":\"" << hexFingerprint(header.contentFingerprint) << "\"},\"sections\":[";
    for (std::size_t i = 0; i < sections.size(); ++i) {
        if (i) out << ',';
        const auto& s = sections[i];
        out << "{\"type\":" << s.type << ",\"name\":\"" << name(s.type) << "\",\"flags\":" << s.flags
            << ",\"offset\":" << s.offset << ",\"size\":" << s.size << ",\"stride\":" << s.stride
            << ",\"count\":" << s.count << ",\"alignment\":" << s.alignment
            << ",\"fingerprint\":\"" << hexFingerprint(s.fingerprint) << "\"}";
    }
    out << "],\"diagnostics\":" << diagnostics.toJson() << '}';
    return out.str();
}
} // namespace nvivo
