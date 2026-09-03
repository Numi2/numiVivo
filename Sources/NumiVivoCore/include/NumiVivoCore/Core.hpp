#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace nvivo {

using Bytes = std::vector<std::byte>;

constexpr std::uint32_t kCompilerABIVersion = 1;
constexpr std::uint16_t kProgramPackMajor = 1;
constexpr std::uint16_t kProgramPackMinor = 0;
constexpr std::size_t kFingerprintBytes = 32;

// MARK: - Diagnostics

enum class Severity : std::uint8_t {
    note = 0,
    warning = 1,
    error = 2,
    fatal = 3
};

struct SourceRange {
    std::size_t offset = 0;
    std::size_t length = 0;
    std::size_t line = 1;
    std::size_t column = 1;
};

struct Diagnostic {
    Severity severity = Severity::error;
    std::string code;
    std::string message;
    std::string path;
    std::string hint;
    SourceRange source;
};

class Diagnostics {
public:
    void add(Diagnostic diagnostic);
    void note(std::string code, std::string message, std::string path = {}, std::string hint = {});
    void warning(std::string code, std::string message, std::string path = {}, std::string hint = {});
    void error(std::string code, std::string message, std::string path = {}, std::string hint = {});
    void fatal(std::string code, std::string message, std::string path = {}, std::string hint = {});

    [[nodiscard]] bool hasErrors() const noexcept;
    [[nodiscard]] bool hasFatal() const noexcept;
    [[nodiscard]] const std::vector<Diagnostic>& entries() const noexcept;
    [[nodiscard]] std::string toJson() const;

    void append(const Diagnostics& other);

private:
    std::vector<Diagnostic> entries_;
};

// MARK: - Bounded JSON representation

namespace json {

class Value {
public:
    using Array = std::vector<Value>;
    using Object = std::map<std::string, Value, std::less<>>;
    using Storage = std::variant<std::nullptr_t, bool, double, std::string, Array, Object>;

    Value() noexcept;
    explicit Value(std::nullptr_t) noexcept;
    explicit Value(bool value) noexcept;
    explicit Value(double value) noexcept;
    explicit Value(std::string value);
    explicit Value(Array value);
    explicit Value(Object value);

    [[nodiscard]] bool isNull() const noexcept;
    [[nodiscard]] bool isBool() const noexcept;
    [[nodiscard]] bool isNumber() const noexcept;
    [[nodiscard]] bool isString() const noexcept;
    [[nodiscard]] bool isArray() const noexcept;
    [[nodiscard]] bool isObject() const noexcept;

    [[nodiscard]] bool asBool(bool fallback = false) const noexcept;
    [[nodiscard]] double asNumber(double fallback = 0.0) const noexcept;
    [[nodiscard]] std::string_view asString() const noexcept;
    [[nodiscard]] const Array& asArray() const;
    [[nodiscard]] const Object& asObject() const;

    [[nodiscard]] const Value* get(std::string_view key) const noexcept;
    [[nodiscard]] const Value* at(std::size_t index) const noexcept;
    [[nodiscard]] const Storage& storage() const noexcept;

private:
    Storage storage_;
};

struct ParseLimits {
    std::size_t maxBytes = 16U * 1024U * 1024U;
    std::size_t maxDepth = 128;
    std::size_t maxNodes = 1'000'000;
    std::size_t maxStringBytes = 4U * 1024U * 1024U;
    std::size_t maxArrayElements = 1'000'000;
    std::size_t maxObjectMembers = 250'000;
};

struct ParseResult {
    std::optional<Value> root;
    Diagnostics diagnostics;
};

[[nodiscard]] ParseResult parse(std::string_view input, const ParseLimits& limits = {});
[[nodiscard]] std::string escape(std::string_view value);

} // namespace json

// MARK: - Units and evidence

struct Dimension {
    std::int8_t length = 0;
    std::int8_t mass = 0;
    std::int8_t time = 0;
    std::int8_t amount = 0;
    std::int8_t temperature = 0;
    std::int8_t count = 0;

    [[nodiscard]] bool operator==(const Dimension&) const noexcept = default;
    [[nodiscard]] bool isDimensionless() const noexcept;
};

struct UnitDefinition {
    std::string symbol;
    Dimension dimension;
    double scaleToSI = 1.0;
    double offsetToSI = 0.0;
    bool contextDependent = false;
};

class UnitRegistry {
public:
    UnitRegistry();

    [[nodiscard]] const UnitDefinition* find(std::string_view symbol) const noexcept;
    [[nodiscard]] bool compatible(std::string_view lhs, std::string_view rhs) const noexcept;
    [[nodiscard]] std::optional<double> convert(double value,
                                                std::string_view from,
                                                std::string_view to) const noexcept;
    [[nodiscard]] std::optional<double> parseDurationSeconds(std::string_view text) const noexcept;

private:
    std::map<std::string, UnitDefinition, std::less<>> units_;
};

enum class EvidenceClass : std::uint8_t {
    observed = 0,
    derived = 1,
    calibrated = 2,
    inferred = 3,
    assumed = 4,
    hypothetical = 5
};

struct Evidence {
    EvidenceClass classification = EvidenceClass::assumed;
    std::string source;
    std::string dataset;
    std::string context;
    std::string note;
};

struct Bounds {
    std::optional<double> minimum;
    std::optional<double> maximum;
};

// MARK: - Program source model

enum class FidelityLevel : std::uint8_t {
    f0Logic = 0,
    f1Deterministic = 1,
    f2Stochastic = 2,
    f3Spatial = 3,
    f4Tissue = 4
};

enum class SignalSource : std::uint8_t {
    intracellular = 0,
    extracellular = 1,
    membrane = 2,
    host = 3,
    numanX = 4,
    numiTissue = 5,
    numiBrain = 6,
    experiment = 7
};

enum class SpeciesKind : std::uint8_t {
    molecularCount = 0,
    concentration = 1,
    activity = 2,
    occupancy = 3,
    output = 4,
    externalField = 5,
    latent = 6
};

enum class StateKind : std::uint8_t {
    scalar = 0,
    leakyIntegrator = 1,
    latch = 2,
    counter = 3,
    timer = 4,
    finiteState = 5,
    permanentMemory = 6
};

enum class RateLawKind : std::uint16_t {
    zeroOrder = 0,
    massAction = 1,
    hillActivation = 2,
    hillRepression = 3,
    michaelisMenten = 4,
    reversibleMassAction = 5,
    passiveTransport = 6,
    saturableTransport = 7,
    degradation = 8,
    customBytecode = 255
};

enum class ExpressionKind : std::uint16_t {
    literal = 0,
    reference = 1,
    logicalNot = 2,
    logicalAll = 3,
    logicalAny = 4,
    greater = 5,
    greaterEqual = 6,
    less = 7,
    lessEqual = 8,
    equal = 9,
    notEqual = 10,
    add = 11,
    subtract = 12,
    multiply = 13,
    divide = 14,
    minimum = 15,
    maximum = 16,
    clamp = 17,
    sustained = 18,
    within = 19,
    risingEdge = 20,
    fallingEdge = 21
};

enum class ReferenceKind : std::uint8_t {
    input = 0,
    species = 1,
    state = 2,
    parameter = 3,
    hostChannel = 4,
    time = 5
};

struct Expression {
    ExpressionKind kind = ExpressionKind::literal;
    ReferenceKind referenceKind = ReferenceKind::input;
    std::string reference;
    std::string unit;
    double value = 0.0;
    double durationSeconds = 0.0;
    std::vector<Expression> children;
};

struct Metadata {
    std::string name;
    std::string version = "0.1.0";
    std::string description;
    std::string namespaceURI;
    std::map<std::string, std::string, std::less<>> labels;
};

struct TargetDefinition {
    std::string cellType;
    std::string tissue;
    std::string species;
    std::string developmentalStage;
    std::string diseaseState;
    std::string deliveryMode;
};

struct InputDefinition {
    std::string id;
    SignalSource source = SignalSource::extracellular;
    std::string unit = "1";
    Bounds bounds;
    double defaultValue = 0.0;
    Evidence evidence;
};

struct SpeciesDefinition {
    std::string id;
    SpeciesKind kind = SpeciesKind::concentration;
    std::string compartment = "cell";
    std::string unit = "nM";
    double initialValue = 0.0;
    Bounds bounds;
    bool conserved = false;
    bool externallyOwned = false;
    Evidence evidence;
};

struct StateDefinition {
    std::string id;
    StateKind kind = StateKind::scalar;
    std::string unit = "1";
    double initialValue = 0.0;
    double halfLifeSeconds = 0.0;
    Bounds bounds;
    std::vector<std::string> states;
};

struct ParameterDefinition {
    std::string id;
    std::string unit = "1";
    double value = 0.0;
    Bounds bounds;
    Evidence evidence;
};

struct StoichiometryTerm {
    std::string species;
    std::int16_t coefficient = 1;
};

struct RateDefinition {
    RateLawKind law = RateLawKind::massAction;
    std::vector<std::string> parameters;
    std::optional<Expression> expression;
};

struct ReactionDefinition {
    std::string id;
    std::string compartment = "cell";
    std::vector<StoichiometryTerm> reactants;
    std::vector<StoichiometryTerm> products;
    RateDefinition rate;
    std::optional<Expression> gate;
    double delaySeconds = 0.0;
    bool critical = false;
};

enum class ActionKind : std::uint16_t {
    setOutput = 0,
    addOutput = 1,
    express = 2,
    suppress = 3,
    degrade = 4,
    setState = 5,
    incrementState = 6,
    emitEvent = 7,
    requestDifferentiation = 8,
    requestMigration = 9,
    reversibleShutdown = 10,
    permanentShutdown = 11
};

struct ActionDefinition {
    ActionKind kind = ActionKind::setOutput;
    std::string target;
    std::optional<Expression> value;
    double constantValue = 0.0;
    std::string unit = "1";
    double maximumRate = 0.0;
};

struct RuleDefinition {
    std::string id;
    Expression condition;
    std::vector<ActionDefinition> actions;
    std::int32_t priority = 0;
    double refractorySeconds = 0.0;
};

enum class ConstraintResponse : std::uint16_t {
    record = 0,
    clamp = 1,
    rejectStep = 2,
    substep = 3,
    reversibleShutdown = 4,
    permanentShutdown = 5
};

struct ConstraintDefinition {
    std::string id;
    Expression condition;
    Severity severity = Severity::error;
    ConstraintResponse response = ConstraintResponse::rejectStep;
    std::string message;
};

struct TerminationDefinition {
    std::string id;
    Expression condition;
    ActionKind action = ActionKind::permanentShutdown;
    std::string reason;
};

struct Program {
    std::string apiVersion;
    std::string kind;
    Metadata metadata;
    TargetDefinition target;
    FidelityLevel minimumFidelity = FidelityLevel::f0Logic;
    std::vector<InputDefinition> inputs;
    std::vector<SpeciesDefinition> species;
    std::vector<StateDefinition> state;
    std::vector<ParameterDefinition> parameters;
    std::vector<ReactionDefinition> reactions;
    std::vector<RuleDefinition> rules;
    std::vector<ConstraintDefinition> constraints;
    std::vector<TerminationDefinition> termination;
};

struct DecodeResult {
    std::optional<Program> program;
    Diagnostics diagnostics;
};

[[nodiscard]] DecodeResult decodeProgram(const json::Value& root, const UnitRegistry& units);

// MARK: - Compiled IR and pack ABI

enum SpeciesFlags : std::uint32_t {
    speciesConserved = 1U << 0U,
    speciesExternallyOwned = 1U << 1U,
    speciesInput = 1U << 2U,
    speciesState = 1U << 3U,
    speciesOutput = 1U << 4U,
    speciesCountValued = 1U << 5U
};

enum ReactionFlags : std::uint32_t {
    reactionCritical = 1U << 0U,
    reactionHasGate = 1U << 1U,
    reactionDelayed = 1U << 2U,
    reactionStochasticEligible = 1U << 3U,
    reactionSpatial = 1U << 4U
};

enum class ExpressionOpcode : std::uint16_t {
    pushConstant = 0,
    loadSpecies = 1,
    loadParameter = 2,
    logicalNot = 3,
    logicalAnd = 4,
    logicalOr = 5,
    greater = 6,
    greaterEqual = 7,
    less = 8,
    lessEqual = 9,
    equal = 10,
    notEqual = 11,
    add = 12,
    subtract = 13,
    multiply = 14,
    divide = 15,
    minimum = 16,
    maximum = 17,
    clamp = 18,
    sustained = 19,
    within = 20,
    risingEdge = 21,
    fallingEdge = 22,
    end = 255
};

struct alignas(16) SpeciesRecord {
    std::uint32_t nameOffset = 0;
    std::uint32_t compartmentOffset = 0;
    std::uint32_t unitOffset = 0;
    std::uint32_t flags = 0;
    float initialValue = 0.0F;
    float minimum = 0.0F;
    float maximum = 0.0F;
    float reserved = 0.0F;
};

struct alignas(16) ParameterRecord {
    std::uint32_t nameOffset = 0;
    std::uint32_t unitOffset = 0;
    std::uint32_t evidenceSourceOffset = 0;
    std::uint32_t flags = 0;
    double value = 0.0;
    double minimum = 0.0;
    double maximum = 0.0;
    std::uint32_t evidenceClass = 0;
    std::uint32_t reserved = 0;
};

struct alignas(8) StoichiometryRecord {
    std::uint32_t speciesIndex = 0;
    std::int16_t coefficient = 0;
    std::uint16_t role = 0; // 0 reactant, 1 product, 2 net incidence
};

struct alignas(16) ReactionRecord {
    std::uint32_t nameOffset = 0;
    std::uint32_t compartmentOffset = 0;
    std::uint32_t reactantOffset = 0;
    std::uint32_t productOffset = 0;
    std::uint32_t reactantCount = 0;
    std::uint32_t productCount = 0;
    std::uint32_t parameterOffset = 0;
    std::uint32_t parameterCount = 0;
    std::uint32_t expressionOffset = 0;
    std::uint32_t expressionCount = 0;
    std::uint32_t rateLaw = 0;
    std::uint32_t flags = 0;
    float delaySeconds = 0.0F;
    float characteristicRate = 0.0F;
    std::uint32_t cohortIndex = 0;
    std::uint32_t reserved = 0;
};

struct alignas(16) ExpressionInstruction {
    std::uint16_t opcode = 0;
    std::uint16_t flags = 0;
    std::uint32_t operand = 0;
    float immediate = 0.0F;
    std::uint32_t auxiliary = 0;
};

struct alignas(16) ActionRecord {
    std::uint32_t targetIndex = 0;
    std::uint32_t expressionOffset = 0;
    std::uint32_t expressionCount = 0;
    std::uint32_t kind = 0;
    float constantValue = 0.0F;
    float maximumRate = 0.0F;
    std::uint32_t unitOffset = 0;
    std::uint32_t flags = 0;
};

struct alignas(16) RuleRecord {
    std::uint32_t nameOffset = 0;
    std::uint32_t conditionOffset = 0;
    std::uint32_t conditionCount = 0;
    std::uint32_t actionOffset = 0;
    std::uint32_t actionCount = 0;
    std::int32_t priority = 0;
    float refractorySeconds = 0.0F;
    std::uint32_t temporalStateOffset = 0;
};

struct alignas(16) MonitorRecord {
    std::uint32_t nameOffset = 0;
    std::uint32_t expressionOffset = 0;
    std::uint32_t expressionCount = 0;
    std::uint32_t messageOffset = 0;
    std::uint32_t severity = 0;
    std::uint32_t response = 0;
    std::uint32_t temporalStateOffset = 0;
    std::uint32_t flags = 0;
};

struct alignas(16) CohortRecord {
    std::uint32_t reactionOffset = 0;
    std::uint32_t reactionCount = 0;
    std::uint32_t rateLaw = 0;
    std::uint32_t flags = 0;
    float maximumStableStep = 0.0F;
    float stiffnessEstimate = 0.0F;
    std::uint32_t preferredThreads = 0;
    std::uint32_t reserved = 0;
};

struct alignas(16) IncidenceRecord {
    std::uint32_t reactionIndex = 0;
    std::int16_t netCoefficient = 0;
    std::uint16_t reserved16 = 0;
    std::uint32_t reserved0 = 0;
    std::uint32_t reserved1 = 0;
};

struct ProgramIR {
    FidelityLevel fidelity = FidelityLevel::f0Logic;
    std::array<std::uint8_t, kFingerprintBytes> sourceFingerprint{};
    std::vector<char> strings;
    std::vector<SpeciesRecord> species;
    std::vector<ParameterRecord> parameters;
    std::vector<std::uint32_t> reactionParameterIndices;
    std::vector<StoichiometryRecord> stoichiometry;
    std::vector<ReactionRecord> reactions;
    std::vector<ExpressionInstruction> expressions;
    std::vector<ActionRecord> actions;
    std::vector<RuleRecord> rules;
    std::vector<MonitorRecord> monitors;
    std::vector<CohortRecord> cohorts;
    std::vector<std::uint32_t> speciesIncidenceOffsets;
    std::vector<IncidenceRecord> speciesIncidence;
    std::uint32_t temporalStateCount = 0;
    std::uint32_t maximumExpressionStack = 0;
    std::uint32_t featureFlags = 0;
};

struct ResourceLimits {
    std::uint32_t maximumSpecies = 16'384;
    std::uint32_t maximumParameters = 65'536;
    std::uint32_t maximumReactions = 65'536;
    std::uint32_t maximumRules = 16'384;
    std::uint32_t maximumConstraints = 16'384;
    std::uint32_t maximumExpressionInstructions = 1'048'576;
    std::uint32_t maximumStoichiometryTerms = 1'048'576;
    std::uint32_t maximumTemporalStates = 262'144;
};

struct CompileOptions {
    FidelityLevel requestedFidelity = FidelityLevel::f2Stochastic;
    bool strictUnits = true;
    bool strictSafety = true;
    bool deterministicPack = true;
    bool permitHypotheticalParameters = true;
    bool requireTermination = true;
    ResourceLimits limits;
};

enum class FindingCategory : std::uint16_t {
    missingBound = 0,
    missingTermination = 1,
    singleSignalActivation = 2,
    irreversibleAction = 3,
    unsupportedContext = 4,
    hypotheticalParameter = 5,
    contradictoryConstraint = 6,
    unmonitoredOutput = 7,
    unboundedAccumulation = 8,
    unitAmbiguity = 9,
    unreachableRule = 10,
    resourceHazard = 11
};

struct SafetyFinding {
    std::string id;
    FindingCategory category = FindingCategory::missingBound;
    Severity severity = Severity::warning;
    std::string subject;
    std::string message;
    std::string requiredMitigation;
};

struct SafetyReport {
    std::vector<SafetyFinding> findings;
    bool hasBlockingFinding = false;
    std::uint32_t monitoredOutputCount = 0;
    std::uint32_t irreversibleActionCount = 0;

    [[nodiscard]] std::string toJson() const;
};

class SafetyAnalyzer {
public:
    [[nodiscard]] SafetyReport analyze(const Program& program,
                                       const CompileOptions& options,
                                       Diagnostics& diagnostics) const;
};

struct CompileResult {
    std::optional<ProgramIR> ir;
    SafetyReport safety;
    Diagnostics diagnostics;
};

class Compiler {
public:
    explicit Compiler(UnitRegistry units = {});

    [[nodiscard]] CompileResult compile(const Program& program,
                                        std::span<const std::byte> source,
                                        const CompileOptions& options) const;

private:
    UnitRegistry units_;
};

enum class PackSectionType : std::uint32_t {
    strings = 1,
    species = 2,
    parameters = 3,
    reactionParameterIndices = 4,
    stoichiometry = 5,
    reactions = 6,
    expressions = 7,
    actions = 8,
    rules = 9,
    monitors = 10,
    cohorts = 11,
    speciesIncidenceOffsets = 12,
    speciesIncidence = 13,
    runtimeContract = 14
};

#pragma pack(push, 1)
struct PackHeader {
    std::array<char, 8> magic{};
    std::uint16_t major = kProgramPackMajor;
    std::uint16_t minor = kProgramPackMinor;
    std::uint32_t headerBytes = 0;
    std::uint32_t compilerABI = kCompilerABIVersion;
    std::uint32_t flags = 0;
    std::uint32_t fidelity = 0;
    std::uint32_t sectionCount = 0;
    std::uint64_t totalBytes = 0;
    std::array<std::uint8_t, kFingerprintBytes> sourceFingerprint{};
    std::array<std::uint8_t, kFingerprintBytes> contentFingerprint{};
    std::array<std::uint8_t, 24> reserved{};
};

struct PackSectionDescriptor {
    std::uint32_t type = 0;
    std::uint32_t flags = 0;
    std::uint64_t offset = 0;
    std::uint64_t size = 0;
    std::uint32_t stride = 0;
    std::uint32_t count = 0;
    std::uint32_t alignment = 0;
    std::uint32_t reserved = 0;
    std::array<std::uint8_t, kFingerprintBytes> fingerprint{};
};

struct RuntimeContractRecord {
    std::uint32_t speciesCount = 0;
    std::uint32_t parameterCount = 0;
    std::uint32_t reactionCount = 0;
    std::uint32_t ruleCount = 0;
    std::uint32_t monitorCount = 0;
    std::uint32_t cohortCount = 0;
    std::uint32_t temporalStateCount = 0;
    std::uint32_t maximumExpressionStack = 0;
    std::uint32_t featureFlags = 0;
    std::uint32_t authoritativeScalarBytes = 4;
    std::uint32_t randomStreamVersion = 1;
    std::uint32_t reserved0 = 0;
    std::array<std::uint64_t, 4> reserved{};
};
#pragma pack(pop)

struct PackBuildResult {
    Bytes bytes;
    std::array<std::uint8_t, kFingerprintBytes> fingerprint{};
    Diagnostics diagnostics;
};

struct PackInspection {
    bool valid = false;
    PackHeader header;
    std::vector<PackSectionDescriptor> sections;
    Diagnostics diagnostics;

    [[nodiscard]] std::string toJson() const;
};

[[nodiscard]] PackBuildResult serializeProgramPack(const ProgramIR& ir);
[[nodiscard]] PackInspection inspectProgramPack(std::span<const std::byte> bytes,
                                                bool verifySectionHashes = true);

[[nodiscard]] std::array<std::uint8_t, kFingerprintBytes>
sha256(std::span<const std::byte> bytes) noexcept;

[[nodiscard]] std::string hexFingerprint(
    const std::array<std::uint8_t, kFingerprintBytes>& fingerprint);

[[nodiscard]] bool isValidIdentifier(std::string_view identifier) noexcept;
[[nodiscard]] std::string_view severityName(Severity severity) noexcept;
[[nodiscard]] std::string_view fidelityName(FidelityLevel fidelity) noexcept;

} // namespace nvivo
