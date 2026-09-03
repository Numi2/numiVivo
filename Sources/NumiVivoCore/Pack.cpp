#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <set>
#include <sstream>
#include <type_traits>

namespace nvivo {

namespace {

constexpr std::array<char, 8> kPackMagic = {'N', 'V', 'I', 'V', 'O', 'P', 'K', '\0'};
constexpr std::uint32_t kSectionRequired = 1U << 0U;
constexpr std::uint64_t kGPUSectionAlignment = 256;
constexpr std::uint64_t kSmallSectionAlignment = 16;
constexpr std::uint32_t kMaximumSectionCount = 1'024;

bool checkedAdd(std::uint64_t left, std::uint64_t right, std::uint64_t& result) {
    if (left > std::numeric_limits<std::uint64_t>::max() - right) return false;
    result = left + right;
    return true;
}

bool checkedMultiply(std::uint64_t left, std::uint64_t right, std::uint64_t& result) {
    if (left != 0 && right > std::numeric_limits<std::uint64_t>::max() / left) return false;
    result = left * right;
    return true;
}

bool isPowerOfTwo(std::uint64_t value) {
    return value != 0 && (value & (value - 1U)) == 0;
}

std::uint64_t alignUp(std::uint64_t value, std::uint64_t alignment, bool& valid) {
    if (!isPowerOfTwo(alignment)) {
        valid = false;
        return value;
    }
    const std::uint64_t mask = alignment - 1U;
    if (value > std::numeric_limits<std::uint64_t>::max() - mask) {
        valid = false;
        return value;
    }
    return (value + mask) & ~mask;
}

template <typename T>
std::span<const std::byte> asBytes(const std::vector<T>& values) {
    static_assert(std::is_trivially_copyable_v<T>);
    if (values.empty()) return {};
    return {
        reinterpret_cast<const std::byte*>(values.data()),
        values.size() * sizeof(T)
    };
}

template <typename T>
std::span<const std::byte> asBytes(const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    return {reinterpret_cast<const std::byte*>(&value), sizeof(T)};
}

std::uint32_t lastStringOffset(const std::vector<char>& strings) {
    if (strings.size() <= 1) return 0;
    std::size_t start = strings.size() - 2;
    while (start > 0 && strings[start - 1] != '\0') --start;
    return start <= std::numeric_limits<std::uint32_t>::max()
        ? static_cast<std::uint32_t>(start)
        : 0;
}

struct SectionInput {
    PackSectionType type;
    std::uint32_t flags;
    std::uint32_t stride;
    std::uint32_t count;
    std::uint32_t alignment;
    std::span<const std::byte> bytes;
};

std::string sectionName(std::uint32_t type) {
    switch (static_cast<PackSectionType>(type)) {
        case PackSectionType::strings: return "strings";
        case PackSectionType::species: return "species";
        case PackSectionType::parameters: return "parameters";
        case PackSectionType::reactionParameterIndices: return "reactionParameterIndices";
        case PackSectionType::stoichiometry: return "stoichiometry";
        case PackSectionType::reactions: return "reactions";
        case PackSectionType::expressions: return "expressions";
        case PackSectionType::actions: return "actions";
        case PackSectionType::rules: return "rules";
        case PackSectionType::monitors: return "monitors";
        case PackSectionType::cohorts: return "cohorts";
        case PackSectionType::speciesIncidenceOffsets: return "speciesIncidenceOffsets";
        case PackSectionType::speciesIncidence: return "speciesIncidence";
        case PackSectionType::runtimeContract: return "runtimeContract";
    }
    return "unknown";
}

std::optional<std::uint32_t> expectedStride(std::uint32_t type) {
    switch (static_cast<PackSectionType>(type)) {
        case PackSectionType::strings: return 1;
        case PackSectionType::species: return sizeof(SpeciesRecord);
        case PackSectionType::parameters: return sizeof(ParameterRecord);
        case PackSectionType::reactionParameterIndices: return sizeof(std::uint32_t);
        case PackSectionType::stoichiometry: return sizeof(StoichiometryRecord);
        case PackSectionType::reactions: return sizeof(ReactionRecord);
        case PackSectionType::expressions: return sizeof(ExpressionInstruction);
        case PackSectionType::actions: return sizeof(ActionRecord);
        case PackSectionType::rules: return sizeof(RuleRecord);
        case PackSectionType::monitors: return sizeof(MonitorRecord);
        case PackSectionType::cohorts: return sizeof(CohortRecord);
        case PackSectionType::speciesIncidenceOffsets: return sizeof(std::uint32_t);
        case PackSectionType::speciesIncidence: return sizeof(IncidenceRecord);
        case PackSectionType::runtimeContract: return sizeof(RuntimeContractRecord);
    }
    return std::nullopt;
}

bool appendAligned(Bytes& output,
                   std::span<const std::byte> payload,
                   std::uint32_t alignment,
                   PackSectionDescriptor& descriptor,
                   Diagnostics& diagnostics) {
    bool valid = true;
    const std::uint64_t aligned = alignUp(output.size(), alignment, valid);
    if (!valid || aligned > std::numeric_limits<std::size_t>::max()) {
        diagnostics.fatal("NVK001", "ProgramPack alignment overflows the host address space.");
        return false;
    }
    output.resize(static_cast<std::size_t>(aligned), std::byte{0});
    descriptor.offset = aligned;
    descriptor.size = payload.size();
    descriptor.alignment = alignment;
    descriptor.fingerprint = sha256(payload);

    if (payload.size() > std::numeric_limits<std::size_t>::max() - output.size()) {
        diagnostics.fatal("NVK002", "ProgramPack payload exceeds the host address space.");
        return false;
    }
    output.insert(output.end(), payload.begin(), payload.end());
    return true;
}

void writeJsonString(std::ostringstream& stream, std::string_view value) {
    stream << '"' << json::escape(value) << '"';
}

} // namespace

PackBuildResult serializeProgramPack(const ProgramIR& ir) {
    PackBuildResult result;

    if (ir.speciesIncidenceOffsets.size() != ir.species.size() + 1U) {
        result.diagnostics.error(
            "NVK003",
            "Species incidence offset table must contain speciesCount + 1 entries."
        );
        return result;
    }

    RuntimeContractRecord contract;
    contract.speciesCount = static_cast<std::uint32_t>(ir.species.size());
    contract.parameterCount = static_cast<std::uint32_t>(ir.parameters.size());
    contract.reactionCount = static_cast<std::uint32_t>(ir.reactions.size());
    contract.ruleCount = static_cast<std::uint32_t>(ir.rules.size());
    contract.monitorCount = static_cast<std::uint32_t>(ir.monitors.size());
    contract.cohortCount = static_cast<std::uint32_t>(ir.cohorts.size());
    contract.temporalStateCount = ir.temporalStateCount;
    contract.maximumExpressionStack = ir.maximumExpressionStack;
    contract.featureFlags = ir.featureFlags;
    contract.authoritativeScalarBytes = sizeof(float);
    contract.randomStreamVersion = 1;
    contract.reserved[0] = lastStringOffset(ir.strings); // ProgramPack v1 manifest string offset.

    const auto stringBytes = ir.strings.empty()
        ? std::span<const std::byte>{}
        : std::span<const std::byte>{
              reinterpret_cast<const std::byte*>(ir.strings.data()),
              ir.strings.size()
          };

    const std::array<SectionInput, 14> inputs = {{
        {PackSectionType::strings, kSectionRequired, 1, static_cast<std::uint32_t>(ir.strings.size()), static_cast<std::uint32_t>(kSmallSectionAlignment), stringBytes},
        {PackSectionType::species, kSectionRequired, sizeof(SpeciesRecord), static_cast<std::uint32_t>(ir.species.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.species)},
        {PackSectionType::parameters, kSectionRequired, sizeof(ParameterRecord), static_cast<std::uint32_t>(ir.parameters.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.parameters)},
        {PackSectionType::reactionParameterIndices, 0, sizeof(std::uint32_t), static_cast<std::uint32_t>(ir.reactionParameterIndices.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.reactionParameterIndices)},
        {PackSectionType::stoichiometry, 0, sizeof(StoichiometryRecord), static_cast<std::uint32_t>(ir.stoichiometry.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.stoichiometry)},
        {PackSectionType::reactions, kSectionRequired, sizeof(ReactionRecord), static_cast<std::uint32_t>(ir.reactions.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.reactions)},
        {PackSectionType::expressions, 0, sizeof(ExpressionInstruction), static_cast<std::uint32_t>(ir.expressions.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.expressions)},
        {PackSectionType::actions, 0, sizeof(ActionRecord), static_cast<std::uint32_t>(ir.actions.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.actions)},
        {PackSectionType::rules, 0, sizeof(RuleRecord), static_cast<std::uint32_t>(ir.rules.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.rules)},
        {PackSectionType::monitors, kSectionRequired, sizeof(MonitorRecord), static_cast<std::uint32_t>(ir.monitors.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.monitors)},
        {PackSectionType::cohorts, kSectionRequired, sizeof(CohortRecord), static_cast<std::uint32_t>(ir.cohorts.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.cohorts)},
        {PackSectionType::speciesIncidenceOffsets, kSectionRequired, sizeof(std::uint32_t), static_cast<std::uint32_t>(ir.speciesIncidenceOffsets.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.speciesIncidenceOffsets)},
        {PackSectionType::speciesIncidence, 0, sizeof(IncidenceRecord), static_cast<std::uint32_t>(ir.speciesIncidence.size()), static_cast<std::uint32_t>(kGPUSectionAlignment), asBytes(ir.speciesIncidence)},
        {PackSectionType::runtimeContract, kSectionRequired, sizeof(RuntimeContractRecord), 1, static_cast<std::uint32_t>(kSmallSectionAlignment), asBytes(contract)}
    }};

    const std::uint64_t descriptorBytes = inputs.size() * sizeof(PackSectionDescriptor);
    std::uint64_t headerBytes = 0;
    if (!checkedAdd(sizeof(PackHeader), descriptorBytes, headerBytes) ||
        headerBytes > std::numeric_limits<std::uint32_t>::max()) {
        result.diagnostics.fatal("NVK004", "ProgramPack header size overflows its ABI field.");
        return result;
    }

    result.bytes.resize(static_cast<std::size_t>(headerBytes), std::byte{0});
    std::vector<PackSectionDescriptor> descriptors;
    descriptors.reserve(inputs.size());

    for (const auto& input : inputs) {
        PackSectionDescriptor descriptor;
        descriptor.type = static_cast<std::uint32_t>(input.type);
        descriptor.flags = input.flags;
        descriptor.stride = input.stride;
        descriptor.count = input.count;
        if (!appendAligned(result.bytes, input.bytes, input.alignment, descriptor, result.diagnostics)) {
            result.bytes.clear();
            return result;
        }
        descriptors.push_back(descriptor);
    }

    PackHeader header;
    header.magic = kPackMagic;
    header.major = kProgramPackMajor;
    header.minor = kProgramPackMinor;
    header.headerBytes = static_cast<std::uint32_t>(headerBytes);
    header.compilerABI = kCompilerABIVersion;
    header.flags = ir.featureFlags;
    header.fidelity = static_cast<std::uint32_t>(ir.fidelity);
    header.sectionCount = static_cast<std::uint32_t>(descriptors.size());
    header.totalBytes = result.bytes.size();
    header.sourceFingerprint = ir.sourceFingerprint;

    std::memcpy(result.bytes.data(), &header, sizeof(header));
    std::memcpy(
        result.bytes.data() + sizeof(PackHeader),
        descriptors.data(),
        descriptors.size() * sizeof(PackSectionDescriptor)
    );

    // The content fingerprint covers every descriptor, alignment byte, and
    // payload byte. Header fields are validated independently and the source
    // document is bound through sourceFingerprint.
    const auto content = std::span<const std::byte>(result.bytes).subspan(sizeof(PackHeader));
    header.contentFingerprint = sha256(content);
    result.fingerprint = header.contentFingerprint;
    std::memcpy(result.bytes.data(), &header, sizeof(header));
    return result;
}

PackInspection inspectProgramPack(std::span<const std::byte> bytes,
                                  bool verifySectionHashes) {
    PackInspection inspection;
    if (bytes.size() < sizeof(PackHeader)) {
        inspection.diagnostics.error("NVK005", "Input is smaller than a ProgramPack header.");
        return inspection;
    }

    std::memcpy(&inspection.header, bytes.data(), sizeof(PackHeader));
    const auto& header = inspection.header;
    if (header.magic != kPackMagic) {
        inspection.diagnostics.error("NVK006", "ProgramPack magic is invalid.");
    }
    if (header.major != kProgramPackMajor) {
        inspection.diagnostics.error(
            "NVK007",
            "Unsupported ProgramPack major version " + std::to_string(header.major) + "."
        );
    }
    if (header.compilerABI > kCompilerABIVersion) {
        inspection.diagnostics.error(
            "NVK008",
            "ProgramPack requires a newer compiler/runtime ABI."
        );
    }
    if (header.sectionCount > kMaximumSectionCount) {
        inspection.diagnostics.error("NVK009", "ProgramPack section count exceeds the hard reader limit.");
        return inspection;
    }

    std::uint64_t descriptorsBytes = 0;
    std::uint64_t minimumHeaderBytes = 0;
    if (!checkedMultiply(header.sectionCount, sizeof(PackSectionDescriptor), descriptorsBytes) ||
        !checkedAdd(sizeof(PackHeader), descriptorsBytes, minimumHeaderBytes)) {
        inspection.diagnostics.error("NVK010", "ProgramPack descriptor table size overflows.");
        return inspection;
    }
    if (header.headerBytes != minimumHeaderBytes || header.headerBytes > bytes.size()) {
        inspection.diagnostics.error("NVK011", "ProgramPack headerBytes does not match its descriptor table.");
        return inspection;
    }
    if (header.totalBytes != bytes.size()) {
        inspection.diagnostics.error(
            "NVK012",
            "ProgramPack totalBytes does not match the supplied byte count."
        );
    }

    inspection.sections.resize(header.sectionCount);
    if (!inspection.sections.empty()) {
        std::memcpy(
            inspection.sections.data(),
            bytes.data() + sizeof(PackHeader),
            inspection.sections.size() * sizeof(PackSectionDescriptor)
        );
    }

    std::set<std::uint32_t> sectionTypes;
    struct Range { std::uint64_t begin; std::uint64_t end; std::uint32_t type; };
    std::vector<Range> ranges;

    for (std::size_t index = 0; index < inspection.sections.size(); ++index) {
        const auto& section = inspection.sections[index];
        const std::string path = "sections[" + std::to_string(index) + "]";

        if (!sectionTypes.insert(section.type).second) {
            inspection.diagnostics.error("NVK013", "Duplicate ProgramPack section type.", path);
        }
        const auto stride = expectedStride(section.type);
        if (!stride.has_value()) {
            if ((section.flags & kSectionRequired) != 0) {
                inspection.diagnostics.error("NVK014", "Unknown required ProgramPack section.", path);
            } else {
                inspection.diagnostics.warning("NVK015", "Unknown optional ProgramPack section ignored.", path);
            }
        } else if (section.stride != *stride) {
            inspection.diagnostics.error(
                "NVK016",
                "ProgramPack section stride does not match the ABI for " + sectionName(section.type) + ".",
                path
            );
        }

        if (!isPowerOfTwo(section.alignment) || section.alignment > 65'536U) {
            inspection.diagnostics.error("NVK017", "ProgramPack section alignment is invalid.", path);
        } else if (section.offset % section.alignment != 0) {
            inspection.diagnostics.error("NVK018", "ProgramPack section offset violates its alignment.", path);
        }
        if (section.offset < header.headerBytes) {
            inspection.diagnostics.error("NVK019", "ProgramPack section overlaps the header.", path);
        }

        std::uint64_t end = 0;
        if (!checkedAdd(section.offset, section.size, end) || end > bytes.size()) {
            inspection.diagnostics.error("NVK020", "ProgramPack section range is out of bounds.", path);
            continue;
        }
        std::uint64_t expectedSize = 0;
        if (!checkedMultiply(section.count, section.stride, expectedSize) || expectedSize != section.size) {
            inspection.diagnostics.error("NVK021", "ProgramPack section count and stride do not match its size.", path);
        }
        if (section.size > 0) ranges.push_back({section.offset, end, section.type});

        if (verifySectionHashes && section.offset <= bytes.size() && end <= bytes.size()) {
            const auto payload = bytes.subspan(
                static_cast<std::size_t>(section.offset),
                static_cast<std::size_t>(section.size)
            );
            if (sha256(payload) != section.fingerprint) {
                inspection.diagnostics.error(
                    "NVK022",
                    "ProgramPack section fingerprint mismatch for " + sectionName(section.type) + ".",
                    path
                );
            }
        }
    }

    std::sort(ranges.begin(), ranges.end(), [](const Range& left, const Range& right) {
        return std::tie(left.begin, left.end, left.type) < std::tie(right.begin, right.end, right.type);
    });
    for (std::size_t index = 1; index < ranges.size(); ++index) {
        if (ranges[index].begin < ranges[index - 1].end) {
            inspection.diagnostics.error("NVK023", "ProgramPack payload sections overlap.");
            break;
        }
    }

    if (bytes.size() >= sizeof(PackHeader)) {
        const auto content = bytes.subspan(sizeof(PackHeader));
        if (sha256(content) != header.contentFingerprint) {
            inspection.diagnostics.error("NVK024", "ProgramPack content fingerprint mismatch.");
        }
    }

    const std::set<std::uint32_t> requiredSections = {
        static_cast<std::uint32_t>(PackSectionType::strings),
        static_cast<std::uint32_t>(PackSectionType::species),
        static_cast<std::uint32_t>(PackSectionType::parameters),
        static_cast<std::uint32_t>(PackSectionType::reactions),
        static_cast<std::uint32_t>(PackSectionType::monitors),
        static_cast<std::uint32_t>(PackSectionType::cohorts),
        static_cast<std::uint32_t>(PackSectionType::speciesIncidenceOffsets),
        static_cast<std::uint32_t>(PackSectionType::runtimeContract)
    };
    for (const auto required : requiredSections) {
        if (!sectionTypes.contains(required)) {
            inspection.diagnostics.error(
                "NVK025",
                "ProgramPack is missing required section " + sectionName(required) + "."
            );
        }
    }

    inspection.valid = !inspection.diagnostics.hasErrors();
    return inspection;
}

std::string PackInspection::toJson() const {
    std::ostringstream stream;
    stream << "{\"valid\":" << (valid ? "true" : "false")
           << ",\"header\":{\"major\":" << header.major
           << ",\"minor\":" << header.minor
           << ",\"compilerABI\":" << header.compilerABI
           << ",\"flags\":" << header.flags
           << ",\"fidelity\":" << header.fidelity
           << ",\"sectionCount\":" << header.sectionCount
           << ",\"totalBytes\":" << header.totalBytes
           << ",\"sourceFingerprint\":";
    writeJsonString(stream, hexFingerprint(header.sourceFingerprint));
    stream << ",\"contentFingerprint\":";
    writeJsonString(stream, hexFingerprint(header.contentFingerprint));
    stream << "},\"sections\":[";

    for (std::size_t index = 0; index < sections.size(); ++index) {
        if (index != 0) stream << ',';
        const auto& section = sections[index];
        stream << "{\"type\":" << section.type << ",\"name\":";
        writeJsonString(stream, sectionName(section.type));
        stream << ",\"flags\":" << section.flags
               << ",\"offset\":" << section.offset
               << ",\"size\":" << section.size
               << ",\"stride\":" << section.stride
               << ",\"count\":" << section.count
               << ",\"alignment\":" << section.alignment
               << ",\"fingerprint\":";
        writeJsonString(stream, hexFingerprint(section.fingerprint));
        stream << '}';
    }
    stream << "],\"diagnostics\":" << diagnostics.toJson() << '}';
    return stream.str();
}

} // namespace nvivo
