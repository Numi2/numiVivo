#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>

namespace nvivo {

namespace {

constexpr std::uint32_t kInvalidIndex = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint32_t kActionTargetIsString = 1U << 0U;
constexpr std::uint32_t kMonitorIsTermination = 1U << 0U;
constexpr std::uint32_t kExpressionReferenceIsTime = 1U << 0U;
constexpr std::uint32_t kFeatureDeterministic = 1U << 0U;
constexpr std::uint32_t kFeatureStochastic = 1U << 1U;
constexpr std::uint32_t kFeatureSpatial = 1U << 2U;
constexpr std::uint32_t kFeatureTissue = 1U << 3U;
constexpr std::uint32_t kFeatureTemporal = 1U << 4U;
constexpr std::uint32_t kFeatureIrreversible = 1U << 5U;
constexpr std::uint32_t kFeatureExternalCoupling = 1U << 6U;

struct Symbol {
    enum class Kind : std::uint8_t { dynamic, parameter };
    Kind kind = Kind::dynamic;
    std::uint32_t index = 0;
    std::string unit;
    bool externallyOwned = false;
    SpeciesKind speciesKind = SpeciesKind::latent;
};

class StringTable {
public:
    explicit StringTable(std::vector<char>& storage) : storage_(storage) {
        storage_.clear();
        storage_.push_back('\0');
        offsets_.emplace("", 0U);
    }

    std::uint32_t intern(std::string_view value, Diagnostics& diagnostics) {
        const auto found = offsets_.find(value);
        if (found != offsets_.end()) {
            return found->second;
        }
        if (storage_.size() > std::numeric_limits<std::uint32_t>::max() - value.size() - 1U) {
            diagnostics.fatal("NVC001", "Compiled string table exceeds the 32-bit pack address space.");
            return 0;
        }
        const auto offset = static_cast<std::uint32_t>(storage_.size());
        storage_.insert(storage_.end(), value.begin(), value.end());
        storage_.push_back('\0');
        offsets_.emplace(std::string(value), offset);
        return offset;
    }

private:
    std::vector<char>& storage_;
    std::map<std::string, std::uint32_t, std::less<>> offsets_;
};

float finiteFloat(double value, double fallback = 0.0) {
    if (!std::isfinite(value)) {
        return static_cast<float>(fallback);
    }
    const double limit = static_cast<double>(std::numeric_limits<float>::max());
    return static_cast<float>(std::clamp(value, -limit, limit));
}

float lowerBound(const Bounds& bounds, bool nonnegative) {
    if (bounds.minimum.has_value()) {
        return finiteFloat(*bounds.minimum, nonnegative ? 0.0 : -std::numeric_limits<float>::max());
    }
    return nonnegative ? 0.0F : -std::numeric_limits<float>::max();
}

float upperBound(const Bounds& bounds) {
    if (bounds.maximum.has_value()) {
        return finiteFloat(*bounds.maximum, std::numeric_limits<float>::max());
    }
    return std::numeric_limits<float>::max();
}

double lowerBoundDouble(const Bounds& bounds) {
    return bounds.minimum.value_or(-std::numeric_limits<double>::max());
}

double upperBoundDouble(const Bounds& bounds) {
    return bounds.maximum.value_or(std::numeric_limits<double>::max());
}

std::optional<std::size_t> requiredParameterCount(RateLawKind law) {
    switch (law) {
        case RateLawKind::zeroOrder: return 1;
        case RateLawKind::massAction: return 1;
        case RateLawKind::hillActivation: return 3;
        case RateLawKind::hillRepression: return 3;
        case RateLawKind::michaelisMenten: return 2;
        case RateLawKind::reversibleMassAction: return 2;
        case RateLawKind::passiveTransport: return 1;
        case RateLawKind::saturableTransport: return 2;
        case RateLawKind::degradation: return 1;
        case RateLawKind::customBytecode: return std::nullopt;
    }
    return std::nullopt;
}

bool isSpatialRateLaw(RateLawKind law) {
    return law == RateLawKind::passiveTransport || law == RateLawKind::saturableTransport;
}

bool isStochasticEligible(RateLawKind law) {
    return law != RateLawKind::customBytecode;
}

int rateBucket(double value) {
    value = std::abs(value);
    if (value < 1.0e-4) return 0;
    if (value < 1.0e-2) return 1;
    if (value < 1.0) return 2;
    if (value < 100.0) return 3;
    if (value < 10'000.0) return 4;
    return 5;
}

std::uint32_t speciesFlags(const SpeciesDefinition& species) {
    std::uint32_t flags = 0;
    if (species.conserved) flags |= speciesConserved;
    if (species.externallyOwned) flags |= speciesExternallyOwned;
    if (species.kind == SpeciesKind::output) flags |= speciesOutput;
    if (species.kind == SpeciesKind::molecularCount) flags |= speciesCountValued;
    return flags;
}

class ExpressionCompiler {
public:
    struct Slice {
        std::uint32_t offset = 0;
        std::uint32_t count = 0;
        std::string unit;
        bool logical = false;
    };

    ExpressionCompiler(ProgramIR& ir,
                       const std::map<std::string, Symbol, std::less<>>& symbols,
                       const UnitRegistry& units,
                       bool strictUnits,
                       const ResourceLimits& limits,
                       Diagnostics& diagnostics)
        : ir_(ir),
          symbols_(symbols),
          units_(units),
          strictUnits_(strictUnits),
          limits_(limits),
          diagnostics_(diagnostics) {}

    Slice append(const Expression& expression, std::string_view path) {
        Slice slice;
        if (ir_.expressions.size() > std::numeric_limits<std::uint32_t>::max()) {
            diagnostics_.fatal("NVC002", "Expression program exceeds 32-bit offsets.", std::string(path));
            return slice;
        }

        slice.offset = static_cast<std::uint32_t>(ir_.expressions.size());
        currentStack_ = 0;
        currentMaximumStack_ = 0;
        const ValueInfo value = compile(expression, path);
        emit(ExpressionOpcode::end, 0, 0, 0.0F, 0, 0);
        slice.count = static_cast<std::uint32_t>(ir_.expressions.size() - slice.offset);
        slice.unit = value.unit;
        slice.logical = value.logical;
        ir_.maximumExpressionStack = std::max(ir_.maximumExpressionStack, currentMaximumStack_);

        if (ir_.expressions.size() > limits_.maximumExpressionInstructions) {
            diagnostics_.fatal(
                "NVC003",
                "Compiled expression instructions exceed the configured resource limit.",
                std::string(path)
            );
        }
        return slice;
    }

private:
    struct ValueInfo {
        std::string unit = "1";
        bool logical = false;
    };

    ValueInfo compile(const Expression& expression, std::string_view path) {
        switch (expression.kind) {
            case ExpressionKind::literal: {
                const std::string unit = expression.unit.empty() ? "1" : expression.unit;
                emit(ExpressionOpcode::pushConstant, 0, 0, finiteFloat(expression.value), 0, +1);
                return {unit, false};
            }
            case ExpressionKind::reference:
                return compileReference(expression, path);
            case ExpressionKind::logicalNot: {
                if (!requireChildren(expression, 1, path)) return {};
                const auto child = compile(expression.children[0], childPath(path, 0));
                validateLogical(child, path);
                emit(ExpressionOpcode::logicalNot, 0, 0, 0.0F, 0, 0);
                return {"1", true};
            }
            case ExpressionKind::logicalAll:
                return compileFold(expression, ExpressionOpcode::logicalAnd, true, path);
            case ExpressionKind::logicalAny:
                return compileFold(expression, ExpressionOpcode::logicalOr, true, path);
            case ExpressionKind::greater:
                return compileBinary(expression, ExpressionOpcode::greater, BinaryClass::comparison, path);
            case ExpressionKind::greaterEqual:
                return compileBinary(expression, ExpressionOpcode::greaterEqual, BinaryClass::comparison, path);
            case ExpressionKind::less:
                return compileBinary(expression, ExpressionOpcode::less, BinaryClass::comparison, path);
            case ExpressionKind::lessEqual:
                return compileBinary(expression, ExpressionOpcode::lessEqual, BinaryClass::comparison, path);
            case ExpressionKind::equal:
                return compileBinary(expression, ExpressionOpcode::equal, BinaryClass::comparison, path);
            case ExpressionKind::notEqual:
                return compileBinary(expression, ExpressionOpcode::notEqual, BinaryClass::comparison, path);
            case ExpressionKind::add:
                return compileBinary(expression, ExpressionOpcode::add, BinaryClass::sameDimension, path);
            case ExpressionKind::subtract:
                return compileBinary(expression, ExpressionOpcode::subtract, BinaryClass::sameDimension, path);
            case ExpressionKind::minimum:
                return compileBinary(expression, ExpressionOpcode::minimum, BinaryClass::sameDimension, path);
            case ExpressionKind::maximum:
                return compileBinary(expression, ExpressionOpcode::maximum, BinaryClass::sameDimension, path);
            case ExpressionKind::multiply:
                return compileBinary(expression, ExpressionOpcode::multiply, BinaryClass::product, path);
            case ExpressionKind::divide:
                return compileBinary(expression, ExpressionOpcode::divide, BinaryClass::quotient, path);
            case ExpressionKind::clamp:
                return compileClamp(expression, path);
            case ExpressionKind::sustained:
                return compileTemporal(expression, ExpressionOpcode::sustained, path);
            case ExpressionKind::within:
                return compileTemporal(expression, ExpressionOpcode::within, path);
            case ExpressionKind::risingEdge:
                return compileEdge(expression, ExpressionOpcode::risingEdge, path);
            case ExpressionKind::fallingEdge:
                return compileEdge(expression, ExpressionOpcode::fallingEdge, path);
        }
        diagnostics_.error("NVC004", "Unsupported expression kind.", std::string(path));
        emit(ExpressionOpcode::pushConstant, 0, 0, 0.0F, 0, +1);
        return {};
    }

    ValueInfo compileReference(const Expression& expression, std::string_view path) {
        if (expression.referenceKind == ReferenceKind::time) {
            emit(
                ExpressionOpcode::loadParameter,
                kExpressionReferenceIsTime,
                kInvalidIndex,
                0.0F,
                0,
                +1
            );
            return {"s", false};
        }

        const auto iterator = symbols_.find(expression.reference);
        if (iterator == symbols_.end()) {
            diagnostics_.error(
                "NVC005",
                "Unresolved expression reference '" + expression.reference + "'.",
                std::string(path),
                "Declare the referenced input, species, state, or parameter."
            );
            emit(ExpressionOpcode::pushConstant, 0, 0, 0.0F, 0, +1);
            return {};
        }

        const Symbol& symbol = iterator->second;
        const auto opcode = symbol.kind == Symbol::Kind::parameter
            ? ExpressionOpcode::loadParameter
            : ExpressionOpcode::loadSpecies;
        emit(opcode, 0, symbol.index, 0.0F, 0, +1);

        if (!expression.unit.empty() && expression.unit != symbol.unit) {
            const auto scale = units_.convert(1.0, symbol.unit, expression.unit);
            if (!scale.has_value()) {
                diagnostics_.error(
                    "NVC006",
                    "Reference unit '" + expression.unit + "' is not compatible with declared unit '" + symbol.unit + "'.",
                    std::string(path)
                );
            } else {
                emit(ExpressionOpcode::pushConstant, 0, 0, finiteFloat(*scale), 0, +1);
                emit(ExpressionOpcode::multiply, 0, 0, 0.0F, 0, -1);
                return {expression.unit, false};
            }
        }
        return {symbol.unit, false};
    }

    enum class BinaryClass { comparison, sameDimension, product, quotient };

    ValueInfo compileBinary(const Expression& expression,
                            ExpressionOpcode opcode,
                            BinaryClass operation,
                            std::string_view path) {
        if (!requireChildren(expression, 2, path)) return {};
        const ValueInfo left = compile(expression.children[0], childPath(path, 0));
        ValueInfo right = compile(expression.children[1], childPath(path, 1));

        if (operation == BinaryClass::comparison || operation == BinaryClass::sameDimension) {
            right = convertTop(right, left.unit, childPath(path, 1));
        }

        emit(opcode, 0, 0, 0.0F, 0, -1);
        if (operation == BinaryClass::comparison) {
            return {"1", true};
        }
        if (operation == BinaryClass::sameDimension) {
            return {left.unit, false};
        }
        if (operation == BinaryClass::product) {
            if (left.unit == "1") return right;
            if (right.unit == "1") return left;
            return {left.unit + "*" + right.unit, false};
        }
        if (right.unit == "1") return left;
        if (left.unit == right.unit) return {"1", false};
        return {left.unit + "/" + right.unit, false};
    }

    ValueInfo compileFold(const Expression& expression,
                          ExpressionOpcode opcode,
                          bool logical,
                          std::string_view path) {
        if (expression.children.empty()) {
            diagnostics_.error("NVC007", "Fold expression has no operands.", std::string(path));
            emit(ExpressionOpcode::pushConstant, 0, 0, 0.0F, 0, +1);
            return {"1", logical};
        }

        ValueInfo result = compile(expression.children.front(), childPath(path, 0));
        if (logical) validateLogical(result, childPath(path, 0));
        for (std::size_t index = 1; index < expression.children.size(); ++index) {
            const auto item = compile(expression.children[index], childPath(path, index));
            if (logical) validateLogical(item, childPath(path, index));
            emit(opcode, 0, 0, 0.0F, 0, -1);
        }
        return logical ? ValueInfo{"1", true} : result;
    }

    ValueInfo compileClamp(const Expression& expression, std::string_view path) {
        if (!requireChildren(expression, 3, path)) return {};
        const auto source = compile(expression.children[0], childPath(path, 0));
        auto minimum = compile(expression.children[1], childPath(path, 1));
        minimum = convertTop(minimum, source.unit, childPath(path, 1));
        auto maximum = compile(expression.children[2], childPath(path, 2));
        maximum = convertTop(maximum, source.unit, childPath(path, 2));
        emit(ExpressionOpcode::clamp, 0, 0, 0.0F, 0, -2);
        return {source.unit, false};
    }

    ValueInfo compileTemporal(const Expression& expression,
                              ExpressionOpcode opcode,
                              std::string_view path) {
        if (!requireChildren(expression, 1, path)) return {};
        const auto condition = compile(expression.children[0], childPath(path, 0));
        validateLogical(condition, childPath(path, 0));
        if (ir_.temporalStateCount >= limits_.maximumTemporalStates) {
            diagnostics_.fatal(
                "NVC008",
                "Temporal state count exceeds the configured resource limit.",
                std::string(path)
            );
            return {"1", true};
        }
        const std::uint32_t stateIndex = ir_.temporalStateCount++;
        emit(opcode, 0, 0, finiteFloat(expression.durationSeconds), stateIndex, 0);
        ir_.featureFlags |= kFeatureTemporal;
        return {"1", true};
    }

    ValueInfo compileEdge(const Expression& expression,
                          ExpressionOpcode opcode,
                          std::string_view path) {
        if (!requireChildren(expression, 1, path)) return {};
        const auto condition = compile(expression.children[0], childPath(path, 0));
        validateLogical(condition, childPath(path, 0));
        if (ir_.temporalStateCount >= limits_.maximumTemporalStates) {
            diagnostics_.fatal(
                "NVC009",
                "Temporal edge state count exceeds the configured resource limit.",
                std::string(path)
            );
            return {"1", true};
        }
        const std::uint32_t stateIndex = ir_.temporalStateCount++;
        emit(opcode, 0, 0, 0.0F, stateIndex, 0);
        ir_.featureFlags |= kFeatureTemporal;
        return {"1", true};
    }

    ValueInfo convertTop(ValueInfo value, std::string_view expected, std::string_view path) {
        if (value.unit == expected) return value;

        if (value.unit == "1" && expected != "1") {
            const std::string message = "Unitless operand is compared or combined with unit '" + std::string(expected) + "'.";
            if (strictUnits_) {
                diagnostics_.error("NVC010", message, std::string(path), "Attach an explicit unit to the literal or reference.");
            } else {
                diagnostics_.warning("NVC010", message, std::string(path), "The operand is interpreted in the other operand's unit.");
            }
            value.unit = std::string(expected);
            return value;
        }

        const auto scale = units_.convert(1.0, value.unit, expected);
        if (!scale.has_value()) {
            diagnostics_.error(
                "NVC011",
                "Expression units '" + value.unit + "' and '" + std::string(expected) + "' are incompatible.",
                std::string(path)
            );
            return value;
        }
        emit(ExpressionOpcode::pushConstant, 0, 0, finiteFloat(*scale), 0, +1);
        emit(ExpressionOpcode::multiply, 0, 0, 0.0F, 0, -1);
        value.unit = std::string(expected);
        return value;
    }

    void validateLogical(const ValueInfo& value, std::string_view path) {
        if (!value.logical && value.unit != "1") {
            diagnostics_.error(
                "NVC012",
                "Logical expression operand has non-dimensionless unit '" + value.unit + "'.",
                std::string(path)
            );
        }
    }

    bool requireChildren(const Expression& expression,
                         std::size_t required,
                         std::string_view path) {
        if (expression.children.size() == required) return true;
        diagnostics_.error(
            "NVC013",
            "Expression has " + std::to_string(expression.children.size()) +
                " operands but requires " + std::to_string(required) + ".",
            std::string(path)
        );
        emit(ExpressionOpcode::pushConstant, 0, 0, 0.0F, 0, +1);
        return false;
    }

    void emit(ExpressionOpcode opcode,
              std::uint16_t flags,
              std::uint32_t operand,
              float immediate,
              std::uint32_t auxiliary,
              int stackDelta) {
        ir_.expressions.push_back({
            static_cast<std::uint16_t>(opcode),
            flags,
            operand,
            immediate,
            auxiliary
        });
        if (stackDelta > 0) {
            currentStack_ += static_cast<std::uint32_t>(stackDelta);
            currentMaximumStack_ = std::max(currentMaximumStack_, currentStack_);
        } else if (stackDelta < 0) {
            const auto decrement = static_cast<std::uint32_t>(-stackDelta);
            currentStack_ = currentStack_ >= decrement ? currentStack_ - decrement : 0;
        }
    }

    static std::string childPath(std::string_view path, std::size_t index) {
        return std::string(path) + ".children[" + std::to_string(index) + "]";
    }

    ProgramIR& ir_;
    const std::map<std::string, Symbol, std::less<>>& symbols_;
    const UnitRegistry& units_;
    bool strictUnits_;
    const ResourceLimits& limits_;
    Diagnostics& diagnostics_;
    std::uint32_t currentStack_ = 0;
    std::uint32_t currentMaximumStack_ = 0;
};

void validateResourceCount(std::size_t actual,
                           std::uint32_t limit,
                           std::string_view name,
                           Diagnostics& diagnostics) {
    if (actual > limit) {
        diagnostics.fatal(
            "NVC014",
            std::string(name) + " count " + std::to_string(actual) +
                " exceeds configured limit " + std::to_string(limit) + "."
        );
    }
    if (actual > std::numeric_limits<std::uint32_t>::max()) {
        diagnostics.fatal("NVC015", std::string(name) + " count exceeds the pack ABI address space.");
    }
}

template <typename Collection, typename Name>
void validateUnique(const Collection& collection,
                    Name&& name,
                    std::string_view path,
                    Diagnostics& diagnostics) {
    std::set<std::string, std::less<>> identifiers;
    for (std::size_t index = 0; index < collection.size(); ++index) {
        const std::string identifier(name(collection[index]));
        if (!identifiers.insert(identifier).second) {
            diagnostics.error(
                "NVC016",
                "Duplicate identifier '" + identifier + "'.",
                std::string(path) + "[" + std::to_string(index) + "].id"
            );
        }
    }
}

std::string manifestJson(const Program& program,
                         const SafetyReport& safety,
                         const ProgramIR& ir) {
    std::ostringstream stream;
    stream << "{\"apiVersion\":\"numivivo.org/pack-manifest/v1\""
           << ",\"program\":{\"name\":\"" << json::escape(program.metadata.name)
           << "\",\"version\":\"" << json::escape(program.metadata.version)
           << "\",\"sourceFingerprint\":\"" << hexFingerprint(ir.sourceFingerprint) << "\"}"
           << ",\"target\":{\"cellType\":\"" << json::escape(program.target.cellType)
           << "\",\"tissue\":\"" << json::escape(program.target.tissue)
           << "\",\"species\":\"" << json::escape(program.target.species)
           << "\",\"developmentalStage\":\"" << json::escape(program.target.developmentalStage)
           << "\",\"diseaseState\":\"" << json::escape(program.target.diseaseState)
           << "\",\"deliveryMode\":\"" << json::escape(program.target.deliveryMode) << "\"}"
           << ",\"fidelity\":\"" << fidelityName(ir.fidelity) << "\""
           << ",\"counts\":{\"species\":" << ir.species.size()
           << ",\"parameters\":" << ir.parameters.size()
           << ",\"reactions\":" << ir.reactions.size()
           << ",\"rules\":" << ir.rules.size()
           << ",\"monitors\":" << ir.monitors.size()
           << ",\"cohorts\":" << ir.cohorts.size()
           << ",\"temporalStates\":" << ir.temporalStateCount << "}"
           << ",\"safety\":" << safety.toJson() << '}';
    return stream.str();
}

} // namespace

Compiler::Compiler(UnitRegistry units) : units_(std::move(units)) {}

CompileResult Compiler::compile(const Program& program,
                                std::span<const std::byte> source,
                                const CompileOptions& options) const {
    CompileResult result;

    validateResourceCount(program.inputs.size() + program.species.size() + program.state.size(),
                          options.limits.maximumSpecies,
                          "Dynamic species",
                          result.diagnostics);
    validateResourceCount(program.parameters.size(), options.limits.maximumParameters, "Parameter", result.diagnostics);
    validateResourceCount(program.reactions.size(), options.limits.maximumReactions, "Reaction", result.diagnostics);
    validateResourceCount(program.rules.size(), options.limits.maximumRules, "Rule", result.diagnostics);
    validateResourceCount(program.constraints.size() + program.termination.size(),
                          options.limits.maximumConstraints,
                          "Monitor",
                          result.diagnostics);

    validateUnique(program.inputs, [](const auto& item) { return item.id; }, "$.spec.inputs", result.diagnostics);
    validateUnique(program.species, [](const auto& item) { return item.id; }, "$.spec.species", result.diagnostics);
    validateUnique(program.state, [](const auto& item) { return item.id; }, "$.spec.state", result.diagnostics);
    validateUnique(program.parameters, [](const auto& item) { return item.id; }, "$.spec.parameters", result.diagnostics);
    validateUnique(program.reactions, [](const auto& item) { return item.id; }, "$.spec.reactions", result.diagnostics);
    validateUnique(program.rules, [](const auto& item) { return item.id; }, "$.spec.rules", result.diagnostics);
    validateUnique(program.constraints, [](const auto& item) { return item.id; }, "$.spec.constraints", result.diagnostics);
    validateUnique(program.termination, [](const auto& item) { return item.id; }, "$.spec.termination", result.diagnostics);

    if (static_cast<std::uint8_t>(options.requestedFidelity) < static_cast<std::uint8_t>(program.minimumFidelity)) {
        result.diagnostics.error(
            "NVC017",
            "Requested fidelity " + std::string(fidelityName(options.requestedFidelity)) +
                " is below the program minimum " + std::string(fidelityName(program.minimumFidelity)) + ".",
            "$.spec.minimumFidelity"
        );
    }

    result.safety = SafetyAnalyzer().analyze(program, options, result.diagnostics);
    if (result.diagnostics.hasErrors()) {
        return result;
    }

    ProgramIR ir;
    ir.fidelity = options.requestedFidelity;
    ir.sourceFingerprint = sha256(source);
    ir.featureFlags = kFeatureDeterministic;
    if (static_cast<std::uint8_t>(ir.fidelity) >= static_cast<std::uint8_t>(FidelityLevel::f2Stochastic)) {
        ir.featureFlags |= kFeatureStochastic;
    }
    if (static_cast<std::uint8_t>(ir.fidelity) >= static_cast<std::uint8_t>(FidelityLevel::f3Spatial)) {
        ir.featureFlags |= kFeatureSpatial;
    }
    if (static_cast<std::uint8_t>(ir.fidelity) >= static_cast<std::uint8_t>(FidelityLevel::f4Tissue)) {
        ir.featureFlags |= kFeatureTissue;
    }

    StringTable strings(ir.strings);
    std::map<std::string, Symbol, std::less<>> symbols;

    auto addDynamic = [&](std::string_view id,
                          std::string_view unit,
                          float initial,
                          float minimum,
                          float maximum,
                          std::uint32_t flags,
                          std::string_view compartment,
                          bool externallyOwned,
                          SpeciesKind kind,
                          std::string_view path) {
        if (symbols.contains(id)) {
            result.diagnostics.error(
                "NVC018",
                "Identifier '" + std::string(id) + "' collides with another dynamic or parameter symbol.",
                std::string(path)
            );
            return;
        }
        if (ir.species.size() >= std::numeric_limits<std::uint32_t>::max()) {
            result.diagnostics.fatal("NVC019", "Dynamic symbol index exceeds the pack ABI.");
            return;
        }
        const auto index = static_cast<std::uint32_t>(ir.species.size());
        ir.species.push_back({
            strings.intern(id, result.diagnostics),
            strings.intern(compartment, result.diagnostics),
            strings.intern(unit, result.diagnostics),
            flags,
            initial,
            minimum,
            maximum,
            0.0F
        });
        symbols.emplace(std::string(id), Symbol{
            Symbol::Kind::dynamic,
            index,
            std::string(unit),
            externallyOwned,
            kind
        });
    };

    for (std::size_t index = 0; index < program.inputs.size(); ++index) {
        const auto& input = program.inputs[index];
        std::uint32_t flags = speciesInput | speciesExternallyOwned;
        if (input.unit == "count" || input.unit == "molecule" || input.unit == "copy") {
            flags |= speciesCountValued;
        }
        addDynamic(
            input.id,
            input.unit,
            finiteFloat(input.defaultValue),
            lowerBound(input.bounds, (flags & speciesCountValued) != 0),
            upperBound(input.bounds),
            flags,
            "external",
            true,
            SpeciesKind::externalField,
            "$.spec.inputs[" + std::to_string(index) + "].id"
        );
        if (input.source == SignalSource::numanX || input.source == SignalSource::numiTissue ||
            input.source == SignalSource::numiBrain || input.source == SignalSource::host) {
            ir.featureFlags |= kFeatureExternalCoupling;
        }
    }

    for (std::size_t index = 0; index < program.species.size(); ++index) {
        const auto& species = program.species[index];
        addDynamic(
            species.id,
            species.unit,
            finiteFloat(species.initialValue),
            lowerBound(species.bounds, species.kind == SpeciesKind::molecularCount),
            upperBound(species.bounds),
            speciesFlags(species),
            species.compartment,
            species.externallyOwned,
            species.kind,
            "$.spec.species[" + std::to_string(index) + "].id"
        );
    }

    for (std::size_t index = 0; index < program.state.size(); ++index) {
        const auto& state = program.state[index];
        std::uint32_t flags = speciesState;
        if (state.kind == StateKind::counter) flags |= speciesCountValued;
        if (state.kind == StateKind::permanentMemory) ir.featureFlags |= kFeatureIrreversible;
        addDynamic(
            state.id,
            state.unit,
            finiteFloat(state.initialValue),
            lowerBound(state.bounds, state.kind == StateKind::counter),
            upperBound(state.bounds),
            flags,
            "program-state",
            false,
            SpeciesKind::latent,
            "$.spec.state[" + std::to_string(index) + "].id"
        );
    }

    for (std::size_t index = 0; index < program.parameters.size(); ++index) {
        const auto& parameter = program.parameters[index];
        if (symbols.contains(parameter.id)) {
            result.diagnostics.error(
                "NVC020",
                "Parameter identifier '" + parameter.id + "' collides with another symbol.",
                "$.spec.parameters[" + std::to_string(index) + "].id"
            );
            continue;
        }
        const auto parameterIndex = static_cast<std::uint32_t>(ir.parameters.size());
        std::uint32_t parameterFlags = 0;
        if (parameter.evidence.classification == EvidenceClass::hypothetical) parameterFlags |= 1U << 0U;
        if (parameter.evidence.classification == EvidenceClass::observed) parameterFlags |= 1U << 1U;
        ir.parameters.push_back({
            strings.intern(parameter.id, result.diagnostics),
            strings.intern(parameter.unit, result.diagnostics),
            strings.intern(parameter.evidence.source, result.diagnostics),
            parameterFlags,
            parameter.value,
            lowerBoundDouble(parameter.bounds),
            upperBoundDouble(parameter.bounds),
            static_cast<std::uint32_t>(parameter.evidence.classification),
            0
        });
        symbols.emplace(parameter.id, Symbol{
            Symbol::Kind::parameter,
            parameterIndex,
            parameter.unit,
            false,
            SpeciesKind::latent
        });
    }

    if (result.diagnostics.hasErrors()) return result;

    ExpressionCompiler expressionCompiler(
        ir,
        symbols,
        units_,
        options.strictUnits,
        options.limits,
        result.diagnostics
    );

    struct OrderedReaction {
        std::size_t sourceIndex = 0;
        std::tuple<std::uint32_t, int, bool, bool, bool> key;
        double characteristicRate = 0.0;
    };
    std::vector<OrderedReaction> orderedReactions;
    orderedReactions.reserve(program.reactions.size());

    for (std::size_t index = 0; index < program.reactions.size(); ++index) {
        const auto& reaction = program.reactions[index];
        const auto required = requiredParameterCount(reaction.rate.law);
        if (required.has_value() && reaction.rate.parameters.size() != *required) {
            result.diagnostics.error(
                "NVC021",
                "Rate law requires " + std::to_string(*required) + " parameter reference(s), but " +
                    std::to_string(reaction.rate.parameters.size()) + " were provided.",
                "$.spec.reactions[" + std::to_string(index) + "].rate.parameters"
            );
        }

        double characteristic = 0.0;
        for (std::size_t parameterPosition = 0; parameterPosition < reaction.rate.parameters.size(); ++parameterPosition) {
            const auto iterator = symbols.find(reaction.rate.parameters[parameterPosition]);
            if (iterator == symbols.end() || iterator->second.kind != Symbol::Kind::parameter) {
                result.diagnostics.error(
                    "NVC022",
                    "Reaction references unknown parameter '" + reaction.rate.parameters[parameterPosition] + "'.",
                    "$.spec.reactions[" + std::to_string(index) + "].rate.parameters[" +
                        std::to_string(parameterPosition) + "]"
                );
            } else if (parameterPosition == 0) {
                characteristic = std::abs(ir.parameters[iterator->second.index].value);
            }
        }

        for (const auto& term : reaction.reactants) {
            const auto iterator = symbols.find(term.species);
            if (iterator == symbols.end() || iterator->second.kind != Symbol::Kind::dynamic) {
                result.diagnostics.error("NVC023", "Unknown reactant species '" + term.species + "'.", "$.spec.reactions[" + std::to_string(index) + "]");
            }
        }
        for (const auto& term : reaction.products) {
            const auto iterator = symbols.find(term.species);
            if (iterator == symbols.end() || iterator->second.kind != Symbol::Kind::dynamic) {
                result.diagnostics.error("NVC024", "Unknown product species '" + term.species + "'.", "$.spec.reactions[" + std::to_string(index) + "]");
            } else if (iterator->second.externallyOwned) {
                result.diagnostics.error(
                    "NVC025",
                    "Reaction writes externally owned species '" + term.species + "'.",
                    "$.spec.reactions[" + std::to_string(index) + "]",
                    "Introduce an internal species and a declared coupling or transport boundary."
                );
            }
        }

        if (isSpatialRateLaw(reaction.rate.law) &&
            static_cast<std::uint8_t>(options.requestedFidelity) < static_cast<std::uint8_t>(FidelityLevel::f3Spatial)) {
            result.diagnostics.error(
                "NVC026",
                "Transport rate law requires F3 or F4 execution.",
                "$.spec.reactions[" + std::to_string(index) + "].rate.law"
            );
        }

        orderedReactions.push_back({
            index,
            {
                static_cast<std::uint32_t>(reaction.rate.law),
                rateBucket(characteristic),
                reaction.gate.has_value(),
                reaction.delaySeconds > 0.0,
                reaction.critical
            },
            characteristic
        });
    }

    if (result.diagnostics.hasErrors()) return result;

    std::stable_sort(orderedReactions.begin(), orderedReactions.end(), [](const auto& left, const auto& right) {
        return left.key < right.key;
    });

    std::vector<std::map<std::uint32_t, std::int32_t>> incidence(ir.species.size());
    std::optional<std::tuple<std::uint32_t, int, bool, bool, bool>> currentCohortKey;
    std::uint32_t currentCohortIndex = 0;

    for (std::size_t orderedIndex = 0; orderedIndex < orderedReactions.size(); ++orderedIndex) {
        const auto& item = orderedReactions[orderedIndex];
        const auto& reaction = program.reactions[item.sourceIndex];

        if (!currentCohortKey.has_value() || *currentCohortKey != item.key) {
            currentCohortKey = item.key;
            currentCohortIndex = static_cast<std::uint32_t>(ir.cohorts.size());
            const double rate = std::max(item.characteristicRate, 1.0e-9);
            CohortRecord cohort;
            cohort.reactionOffset = static_cast<std::uint32_t>(ir.reactions.size());
            cohort.rateLaw = static_cast<std::uint32_t>(reaction.rate.law);
            cohort.flags = 0;
            if (reaction.gate.has_value()) cohort.flags |= reactionHasGate;
            if (reaction.delaySeconds > 0.0) cohort.flags |= reactionDelayed;
            if (isSpatialRateLaw(reaction.rate.law)) cohort.flags |= reactionSpatial;
            cohort.maximumStableStep = finiteFloat(std::clamp(0.1 / rate, 1.0e-7, 60.0));
            cohort.stiffnessEstimate = finiteFloat(rate);
            cohort.preferredThreads = 256;
            ir.cohorts.push_back(cohort);
        }
        ++ir.cohorts[currentCohortIndex].reactionCount;

        ReactionRecord record;
        record.nameOffset = strings.intern(reaction.id, result.diagnostics);
        record.compartmentOffset = strings.intern(reaction.compartment, result.diagnostics);
        record.reactantOffset = static_cast<std::uint32_t>(ir.stoichiometry.size());
        record.reactantCount = static_cast<std::uint32_t>(reaction.reactants.size());
        for (const auto& term : reaction.reactants) {
            const auto symbol = symbols.find(term.species);
            ir.stoichiometry.push_back({symbol->second.index, term.coefficient, 0});
            incidence[symbol->second.index][static_cast<std::uint32_t>(ir.reactions.size())] -= term.coefficient;
        }
        record.productOffset = static_cast<std::uint32_t>(ir.stoichiometry.size());
        record.productCount = static_cast<std::uint32_t>(reaction.products.size());
        for (const auto& term : reaction.products) {
            const auto symbol = symbols.find(term.species);
            ir.stoichiometry.push_back({symbol->second.index, term.coefficient, 1});
            incidence[symbol->second.index][static_cast<std::uint32_t>(ir.reactions.size())] += term.coefficient;
        }

        record.parameterOffset = static_cast<std::uint32_t>(ir.reactionParameterIndices.size());
        record.parameterCount = static_cast<std::uint32_t>(reaction.rate.parameters.size());
        for (const auto& parameter : reaction.rate.parameters) {
            ir.reactionParameterIndices.push_back(symbols.find(parameter)->second.index);
        }

        record.expressionOffset = 0;
        record.expressionCount = 0;
        if (reaction.rate.expression.has_value()) {
            const auto slice = expressionCompiler.append(
                *reaction.rate.expression,
                "$.spec.reactions[" + std::to_string(item.sourceIndex) + "].rate.expression"
            );
            record.expressionOffset = slice.offset;
            record.expressionCount = slice.count;
        }

        // In ProgramPack v1, `reserved` is the gate bytecode offset. The gate
        // program is terminated by the `end` opcode; kInvalidIndex means absent.
        record.reserved = kInvalidIndex;
        if (reaction.gate.has_value()) {
            const auto gate = expressionCompiler.append(
                *reaction.gate,
                "$.spec.reactions[" + std::to_string(item.sourceIndex) + "].gate"
            );
            record.reserved = gate.offset;
        }

        record.rateLaw = static_cast<std::uint32_t>(reaction.rate.law);
        record.flags = 0;
        if (reaction.critical) record.flags |= reactionCritical;
        if (reaction.gate.has_value()) record.flags |= reactionHasGate;
        if (reaction.delaySeconds > 0.0) record.flags |= reactionDelayed;
        if (isStochasticEligible(reaction.rate.law)) record.flags |= reactionStochasticEligible;
        if (isSpatialRateLaw(reaction.rate.law)) record.flags |= reactionSpatial;
        record.delaySeconds = finiteFloat(reaction.delaySeconds);
        record.characteristicRate = finiteFloat(item.characteristicRate);
        record.cohortIndex = currentCohortIndex;
        ir.reactions.push_back(record);
    }

    validateResourceCount(ir.stoichiometry.size(), options.limits.maximumStoichiometryTerms, "Stoichiometry term", result.diagnostics);

    std::vector<std::size_t> orderedRules(program.rules.size());
    for (std::size_t index = 0; index < orderedRules.size(); ++index) orderedRules[index] = index;
    std::stable_sort(orderedRules.begin(), orderedRules.end(), [&program](std::size_t left, std::size_t right) {
        const auto& a = program.rules[left];
        const auto& b = program.rules[right];
        return std::tie(a.priority, b.id) > std::tie(b.priority, a.id);
    });

    for (const std::size_t sourceIndex : orderedRules) {
        const auto& rule = program.rules[sourceIndex];
        RuleRecord record;
        record.nameOffset = strings.intern(rule.id, result.diagnostics);
        record.temporalStateOffset = ir.temporalStateCount;
        const auto condition = expressionCompiler.append(
            rule.condition,
            "$.spec.rules[" + std::to_string(sourceIndex) + "].when"
        );
        record.conditionOffset = condition.offset;
        record.conditionCount = condition.count;
        record.actionOffset = static_cast<std::uint32_t>(ir.actions.size());
        record.actionCount = static_cast<std::uint32_t>(rule.actions.size());
        record.priority = rule.priority;
        record.refractorySeconds = finiteFloat(rule.refractorySeconds);

        for (std::size_t actionIndex = 0; actionIndex < rule.actions.size(); ++actionIndex) {
            const auto& action = rule.actions[actionIndex];
            ActionRecord actionRecord;
            actionRecord.kind = static_cast<std::uint32_t>(action.kind);
            actionRecord.constantValue = finiteFloat(action.constantValue);
            actionRecord.maximumRate = finiteFloat(action.maximumRate);
            actionRecord.unitOffset = strings.intern(action.unit, result.diagnostics);

            const bool stringTarget = action.kind == ActionKind::emitEvent ||
                                      action.kind == ActionKind::requestDifferentiation ||
                                      action.kind == ActionKind::requestMigration;
            const bool targetless = action.kind == ActionKind::reversibleShutdown ||
                                    action.kind == ActionKind::permanentShutdown;
            if (targetless) {
                actionRecord.targetIndex = kInvalidIndex;
            } else if (stringTarget) {
                actionRecord.targetIndex = strings.intern(action.target, result.diagnostics);
                actionRecord.flags |= kActionTargetIsString;
            } else {
                const auto target = symbols.find(action.target);
                if (target == symbols.end() || target->second.kind != Symbol::Kind::dynamic) {
                    result.diagnostics.error(
                        "NVC027",
                        "Action references unknown dynamic target '" + action.target + "'.",
                        "$.spec.rules[" + std::to_string(sourceIndex) + "].actions[" +
                            std::to_string(actionIndex) + "].target"
                    );
                    actionRecord.targetIndex = kInvalidIndex;
                } else if (target->second.externallyOwned) {
                    result.diagnostics.error(
                        "NVC028",
                        "Action writes externally owned target '" + action.target + "'.",
                        "$.spec.rules[" + std::to_string(sourceIndex) + "].actions[" +
                            std::to_string(actionIndex) + "].target"
                    );
                    actionRecord.targetIndex = target->second.index;
                } else {
                    actionRecord.targetIndex = target->second.index;
                }
            }

            if (action.value.has_value()) {
                const auto expression = expressionCompiler.append(
                    *action.value,
                    "$.spec.rules[" + std::to_string(sourceIndex) + "].actions[" +
                        std::to_string(actionIndex) + "].value"
                );
                actionRecord.expressionOffset = expression.offset;
                actionRecord.expressionCount = expression.count;
            }
            if (action.kind == ActionKind::permanentShutdown ||
                action.kind == ActionKind::requestDifferentiation) {
                ir.featureFlags |= kFeatureIrreversible;
            }
            ir.actions.push_back(actionRecord);
        }
        ir.rules.push_back(record);
    }

    for (std::size_t index = 0; index < program.constraints.size(); ++index) {
        const auto& constraint = program.constraints[index];
        MonitorRecord record;
        record.nameOffset = strings.intern(constraint.id, result.diagnostics);
        record.temporalStateOffset = ir.temporalStateCount;
        const auto expression = expressionCompiler.append(
            constraint.condition,
            "$.spec.constraints[" + std::to_string(index) + "].expression"
        );
        record.expressionOffset = expression.offset;
        record.expressionCount = expression.count;
        record.messageOffset = strings.intern(constraint.message, result.diagnostics);
        record.severity = static_cast<std::uint32_t>(constraint.severity);
        record.response = static_cast<std::uint32_t>(constraint.response);
        ir.monitors.push_back(record);
    }

    for (std::size_t index = 0; index < program.termination.size(); ++index) {
        const auto& termination = program.termination[index];
        MonitorRecord record;
        record.nameOffset = strings.intern(termination.id, result.diagnostics);
        record.temporalStateOffset = ir.temporalStateCount;
        const auto expression = expressionCompiler.append(
            termination.condition,
            "$.spec.termination[" + std::to_string(index) + "].when"
        );
        record.expressionOffset = expression.offset;
        record.expressionCount = expression.count;
        record.messageOffset = strings.intern(termination.reason, result.diagnostics);
        record.severity = static_cast<std::uint32_t>(Severity::fatal);
        record.response = termination.action == ActionKind::permanentShutdown
            ? static_cast<std::uint32_t>(ConstraintResponse::permanentShutdown)
            : static_cast<std::uint32_t>(ConstraintResponse::reversibleShutdown);
        record.flags = kMonitorIsTermination;
        ir.monitors.push_back(record);
    }

    ir.speciesIncidenceOffsets.reserve(ir.species.size() + 1U);
    ir.speciesIncidenceOffsets.push_back(0);
    for (std::size_t speciesIndex = 0; speciesIndex < incidence.size(); ++speciesIndex) {
        for (const auto& [reactionIndex, coefficient] : incidence[speciesIndex]) {
            if (coefficient == 0) continue;
            if (coefficient < std::numeric_limits<std::int16_t>::min() ||
                coefficient > std::numeric_limits<std::int16_t>::max()) {
                result.diagnostics.error(
                    "NVC029",
                    "Net stoichiometric coefficient exceeds the 16-bit GPU ABI.",
                    "species:" + std::to_string(speciesIndex)
                );
                continue;
            }
            ir.speciesIncidence.push_back({
                reactionIndex,
                static_cast<std::int16_t>(coefficient),
                0,
                0,
                0
            });
        }
        ir.speciesIncidenceOffsets.push_back(static_cast<std::uint32_t>(ir.speciesIncidence.size()));
    }

    validateResourceCount(ir.expressions.size(), options.limits.maximumExpressionInstructions, "Expression instruction", result.diagnostics);
    validateResourceCount(ir.temporalStateCount, options.limits.maximumTemporalStates, "Temporal state", result.diagnostics);

    // Preserve a compact manifest in the string table. Offset zero remains the
    // empty string; the manifest is the final string and is referenced by the
    // RuntimeContract extension in the pack serializer.
    const std::string manifest = manifestJson(program, result.safety, ir);
    strings.intern(manifest, result.diagnostics);

    if (!result.diagnostics.hasErrors()) {
        result.ir = std::move(ir);
    }
    return result;
}

} // namespace nvivo
