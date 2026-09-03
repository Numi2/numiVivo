#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace nvivo {

namespace {

class ProgramDecoder {
public:
    ProgramDecoder(const json::Value& root, const UnitRegistry& units)
        : root_(root), units_(units) {}

    DecodeResult run() {
        DecodeResult result;
        if (!root_.isObject()) {
            diagnostics_.error("NVP001", "VivoProgram root must be a JSON object.", "$", {});
            result.diagnostics = std::move(diagnostics_);
            return result;
        }

        Program program;
        program.apiVersion = requireString(root_, "apiVersion", "$.apiVersion");
        program.kind = requireString(root_, "kind", "$.kind");
        if (program.apiVersion != "numivivo.org/v1alpha1") {
            diagnostics_.error(
                "NVP002",
                "Unsupported apiVersion '" + program.apiVersion + "'.",
                "$.apiVersion",
                "Use numivivo.org/v1alpha1."
            );
        }
        if (program.kind != "VivoProgram") {
            diagnostics_.error(
                "NVP003",
                "Document kind must be VivoProgram.",
                "$.kind",
                {}
            );
        }

        if (const auto* metadata = root_.get("metadata")) {
            decodeMetadata(*metadata, program.metadata, "$.metadata");
        } else {
            diagnostics_.error("NVP004", "metadata is required.", "$.metadata", {});
        }

        const auto* specification = root_.get("spec");
        if (specification == nullptr || !specification->isObject()) {
            diagnostics_.error("NVP005", "spec must be an object.", "$.spec", {});
        } else {
            decodeSpecification(*specification, program);
        }

        if (program.metadata.name.empty()) {
            diagnostics_.error("NVP006", "metadata.name cannot be empty.", "$.metadata.name", {});
        } else if (!isValidIdentifier(program.metadata.name)) {
            diagnostics_.error(
                "NVP007",
                "metadata.name is not a valid NumiVivo identifier.",
                "$.metadata.name",
                "Use letters, digits, '_', '-', or '.', beginning with a letter or underscore."
            );
        }

        if (!diagnostics_.hasErrors()) {
            result.program = std::move(program);
        }
        result.diagnostics = std::move(diagnostics_);
        return result;
    }

private:
    void decodeMetadata(const json::Value& value, Metadata& metadata, std::string_view path) {
        if (!value.isObject()) {
            diagnostics_.error("NVP008", "metadata must be an object.", std::string(path), {});
            return;
        }
        metadata.name = requireString(value, "name", join(path, "name"));
        metadata.version = optionalString(value, "version", "0.1.0", join(path, "version"));
        metadata.description = optionalString(value, "description", {}, join(path, "description"));
        metadata.namespaceURI = optionalString(value, "namespace", {}, join(path, "namespace"));

        if (const auto* labels = value.get("labels")) {
            if (!labels->isObject()) {
                diagnostics_.error("NVP009", "metadata.labels must be an object.", join(path, "labels"), {});
            } else {
                for (const auto& [key, labelValue] : labels->asObject()) {
                    if (!labelValue.isString()) {
                        diagnostics_.error(
                            "NVP010",
                            "metadata label values must be strings.",
                            join(join(path, "labels"), key),
                            {}
                        );
                        continue;
                    }
                    metadata.labels.emplace(key, std::string(labelValue.asString()));
                }
            }
        }
    }

    void decodeSpecification(const json::Value& specification, Program& program) {
        if (const auto* target = specification.get("target")) {
            decodeTarget(*target, program.target, "$.spec.target");
        } else {
            diagnostics_.error("NVP011", "spec.target is required.", "$.spec.target", {});
        }

        if (const auto* fidelity = specification.get("minimumFidelity")) {
            program.minimumFidelity = decodeFidelity(*fidelity, "$.spec.minimumFidelity");
        }

        decodeArray(specification, "inputs", "$.spec.inputs", [this, &program](const json::Value& item, const std::string& path) {
            InputDefinition definition;
            decodeInput(item, definition, path);
            program.inputs.push_back(std::move(definition));
        });

        decodeArray(specification, "species", "$.spec.species", [this, &program](const json::Value& item, const std::string& path) {
            SpeciesDefinition definition;
            decodeSpecies(item, definition, path);
            program.species.push_back(std::move(definition));
        });

        decodeArray(specification, "state", "$.spec.state", [this, &program](const json::Value& item, const std::string& path) {
            StateDefinition definition;
            decodeState(item, definition, path);
            program.state.push_back(std::move(definition));
        });

        decodeArray(specification, "parameters", "$.spec.parameters", [this, &program](const json::Value& item, const std::string& path) {
            ParameterDefinition definition;
            decodeParameter(item, definition, path);
            program.parameters.push_back(std::move(definition));
        });

        decodeArray(specification, "reactions", "$.spec.reactions", [this, &program](const json::Value& item, const std::string& path) {
            ReactionDefinition definition;
            decodeReaction(item, definition, path);
            program.reactions.push_back(std::move(definition));
        });

        decodeArray(specification, "rules", "$.spec.rules", [this, &program](const json::Value& item, const std::string& path) {
            RuleDefinition definition;
            decodeRule(item, definition, path);
            program.rules.push_back(std::move(definition));
        });

        decodeArray(specification, "constraints", "$.spec.constraints", [this, &program](const json::Value& item, const std::string& path) {
            ConstraintDefinition definition;
            decodeConstraint(item, definition, path);
            program.constraints.push_back(std::move(definition));
        });

        decodeArray(specification, "termination", "$.spec.termination", [this, &program](const json::Value& item, const std::string& path) {
            TerminationDefinition definition;
            decodeTermination(item, definition, path);
            program.termination.push_back(std::move(definition));
        });

        if (program.reactions.empty() && program.rules.empty()) {
            diagnostics_.warning(
                "NVP012",
                "Program contains neither reactions nor behavioral rules.",
                "$.spec",
                "Add at least one reaction or rule to produce executable behavior."
            );
        }
    }

    void decodeTarget(const json::Value& value, TargetDefinition& target, std::string_view path) {
        if (!value.isObject()) {
            diagnostics_.error("NVP013", "target must be an object.", std::string(path), {});
            return;
        }
        target.cellType = requireString(value, "cellType", join(path, "cellType"));
        target.tissue = optionalString(value, "tissue", {}, join(path, "tissue"));
        target.species = optionalString(value, "species", {}, join(path, "species"));
        target.developmentalStage = optionalString(value, "developmentalStage", {}, join(path, "developmentalStage"));
        target.diseaseState = optionalString(value, "diseaseState", {}, join(path, "diseaseState"));
        target.deliveryMode = optionalString(value, "deliveryMode", {}, join(path, "deliveryMode"));
    }

    void decodeInput(const json::Value& value, InputDefinition& input, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        input.id = decodeIdentifier(value, "id", join(path, "id"));
        input.source = decodeSignalSource(optionalString(value, "source", "extracellular", join(path, "source")), join(path, "source"));
        input.unit = optionalString(value, "unit", "1", join(path, "unit"));
        validateUnit(input.unit, join(path, "unit"));
        input.defaultValue = optionalNumber(value, "default", 0.0, join(path, "default"));
        input.bounds = decodeBounds(value.get("bounds"), join(path, "bounds"));
        decodeEvidence(value.get("evidence"), input.evidence, join(path, "evidence"));
        validateBounds(input.defaultValue, input.bounds, path);
    }

    void decodeSpecies(const json::Value& value, SpeciesDefinition& species, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        species.id = decodeIdentifier(value, "id", join(path, "id"));
        species.kind = decodeSpeciesKind(optionalString(value, "kind", "concentration", join(path, "kind")), join(path, "kind"));
        species.compartment = optionalString(value, "compartment", "cell", join(path, "compartment"));
        species.unit = optionalString(value, "unit", species.kind == SpeciesKind::molecularCount ? "count" : "nM", join(path, "unit"));
        validateUnit(species.unit, join(path, "unit"));
        species.initialValue = decodeQuantity(value.get("initial"), 0.0, species.unit, join(path, "initial"));
        species.bounds = decodeBounds(value.get("bounds"), join(path, "bounds"));
        species.conserved = optionalBool(value, "conserved", false, join(path, "conserved"));
        species.externallyOwned = optionalBool(value, "externallyOwned", false, join(path, "externallyOwned"));
        decodeEvidence(value.get("evidence"), species.evidence, join(path, "evidence"));
        validateBounds(species.initialValue, species.bounds, path);
    }

    void decodeState(const json::Value& value, StateDefinition& state, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        state.id = decodeIdentifier(value, "id", join(path, "id"));
        state.kind = decodeStateKind(optionalString(value, "type", "scalar", join(path, "type")), join(path, "type"));
        state.unit = optionalString(value, "unit", "1", join(path, "unit"));
        validateUnit(state.unit, join(path, "unit"));
        state.initialValue = decodeQuantity(value.get("initial"), 0.0, state.unit, join(path, "initial"));
        state.bounds = decodeBounds(value.get("bounds"), join(path, "bounds"));
        if (const auto* halfLife = value.get("halfLife")) {
            state.halfLifeSeconds = decodeDuration(*halfLife, join(path, "halfLife"));
        }
        if (const auto* states = value.get("states")) {
            if (!states->isArray()) {
                diagnostics_.error("NVP014", "state.states must be an array of strings.", join(path, "states"), {});
            } else {
                for (std::size_t index = 0; index < states->asArray().size(); ++index) {
                    const auto& item = states->asArray()[index];
                    if (!item.isString()) {
                        diagnostics_.error("NVP015", "Finite-state labels must be strings.", indexPath(join(path, "states"), index), {});
                    } else {
                        state.states.emplace_back(item.asString());
                    }
                }
            }
        }
        if (state.kind == StateKind::leakyIntegrator && state.halfLifeSeconds <= 0.0) {
            diagnostics_.error(
                "NVP016",
                "A leaky-integrator state requires a positive halfLife.",
                join(path, "halfLife"),
                {}
            );
        }
        if (state.kind == StateKind::finiteState && state.states.empty()) {
            diagnostics_.error(
                "NVP017",
                "A finite-state variable requires at least one state label.",
                join(path, "states"),
                {}
            );
        }
        validateBounds(state.initialValue, state.bounds, path);
    }

    void decodeParameter(const json::Value& value, ParameterDefinition& parameter, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        parameter.id = decodeIdentifier(value, "id", join(path, "id"));
        parameter.unit = optionalString(value, "unit", "1", join(path, "unit"));
        validateUnit(parameter.unit, join(path, "unit"));
        parameter.value = decodeQuantity(value.get("value"), 0.0, parameter.unit, join(path, "value"), true);
        parameter.bounds = decodeBounds(value.get("bounds"), join(path, "bounds"));
        decodeEvidence(value.get("evidence"), parameter.evidence, join(path, "evidence"));
        validateBounds(parameter.value, parameter.bounds, path);
    }

    void decodeReaction(const json::Value& value, ReactionDefinition& reaction, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        reaction.id = decodeIdentifier(value, "id", join(path, "id"));
        reaction.compartment = optionalString(value, "compartment", "cell", join(path, "compartment"));
        decodeStoichiometry(value.get("reactants"), reaction.reactants, join(path, "reactants"));
        decodeStoichiometry(value.get("products"), reaction.products, join(path, "products"));
        reaction.critical = optionalBool(value, "critical", false, join(path, "critical"));

        if (const auto* delay = value.get("delay")) {
            reaction.delaySeconds = decodeDuration(*delay, join(path, "delay"));
        }
        if (const auto* gate = value.get("gate")) {
            reaction.gate = decodeExpression(*gate, join(path, "gate"), 0);
        }

        const auto* rate = value.get("rate");
        if (rate == nullptr || !rate->isObject()) {
            diagnostics_.error("NVP018", "reaction.rate must be an object.", join(path, "rate"), {});
            return;
        }
        reaction.rate.law = decodeRateLaw(requireString(*rate, "law", join(join(path, "rate"), "law")), join(join(path, "rate"), "law"));
        if (const auto* parameters = rate->get("parameters")) {
            if (!parameters->isArray()) {
                diagnostics_.error("NVP019", "rate.parameters must be an array.", join(join(path, "rate"), "parameters"), {});
            } else {
                for (std::size_t index = 0; index < parameters->asArray().size(); ++index) {
                    const auto& item = parameters->asArray()[index];
                    if (!item.isString()) {
                        diagnostics_.error("NVP020", "Reaction parameter references must be strings.", indexPath(join(join(path, "rate"), "parameters"), index), {});
                    } else {
                        reaction.rate.parameters.emplace_back(item.asString());
                    }
                }
            }
        }
        if (const auto* expression = rate->get("expression")) {
            reaction.rate.expression = decodeExpression(*expression, join(join(path, "rate"), "expression"), 0);
        }

        if (reaction.reactants.empty() && reaction.products.empty()) {
            diagnostics_.error("NVP021", "Reaction must contain at least one reactant or product.", path, {});
        }
        if (reaction.rate.law == RateLawKind::customBytecode && !reaction.rate.expression.has_value()) {
            diagnostics_.error("NVP022", "custom rate law requires rate.expression.", join(join(path, "rate"), "expression"), {});
        }
    }

    void decodeRule(const json::Value& value, RuleDefinition& rule, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        rule.id = decodeIdentifier(value, "id", join(path, "id"));
        const json::Value* condition = value.get("when");
        if (condition == nullptr) {
            condition = value.get("condition");
        }
        if (condition == nullptr) {
            diagnostics_.error("NVP023", "Rule requires a when condition.", join(path, "when"), {});
        } else {
            rule.condition = decodeExpression(*condition, join(path, "when"), 0);
        }
        rule.priority = static_cast<std::int32_t>(optionalInteger(value, "priority", 0, join(path, "priority")));
        if (const auto* refractory = value.get("refractory")) {
            rule.refractorySeconds = decodeDuration(*refractory, join(path, "refractory"));
        }

        const auto* actions = value.get("actions");
        if (actions == nullptr) {
            actions = value.get("then");
        }
        if (actions == nullptr || !actions->isArray()) {
            diagnostics_.error("NVP024", "Rule requires an actions array.", join(path, "actions"), {});
            return;
        }
        if (actions->asArray().empty()) {
            diagnostics_.error("NVP025", "Rule actions cannot be empty.", join(path, "actions"), {});
        }
        for (std::size_t index = 0; index < actions->asArray().size(); ++index) {
            ActionDefinition action;
            decodeAction(actions->asArray()[index], action, indexPath(join(path, "actions"), index));
            rule.actions.push_back(std::move(action));
        }
    }

    void decodeAction(const json::Value& value, ActionDefinition& action, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        action.kind = decodeActionKind(requireString(value, "kind", join(path, "kind")), join(path, "kind"));
        action.target = optionalString(value, "target", {}, join(path, "target"));
        action.unit = optionalString(value, "unit", "1", join(path, "unit"));
        validateUnit(action.unit, join(path, "unit"));
        action.constantValue = optionalNumber(value, "constant", 0.0, join(path, "constant"));
        action.maximumRate = optionalNumber(value, "maximumRate", 0.0, join(path, "maximumRate"));
        if (const auto* expression = value.get("value")) {
            action.value = decodeExpression(*expression, join(path, "value"), 0);
        }

        const bool targetless = action.kind == ActionKind::reversibleShutdown ||
                                action.kind == ActionKind::permanentShutdown;
        if (!targetless && action.target.empty()) {
            diagnostics_.error("NVP026", "Action target is required for this action kind.", join(path, "target"), {});
        }
        if (action.maximumRate < 0.0) {
            diagnostics_.error("NVP027", "Action maximumRate cannot be negative.", join(path, "maximumRate"), {});
        }
    }

    void decodeConstraint(const json::Value& value, ConstraintDefinition& constraint, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        constraint.id = decodeIdentifier(value, "id", join(path, "id"));
        const json::Value* expression = value.get("condition");
        if (expression == nullptr) {
            expression = value.get("expression");
        }
        if (expression == nullptr) {
            diagnostics_.error("NVP028", "Constraint requires an expression.", join(path, "expression"), {});
        } else {
            constraint.condition = decodeExpression(*expression, join(path, "expression"), 0);
        }
        constraint.severity = decodeSeverity(optionalString(value, "severity", "error", join(path, "severity")), join(path, "severity"));
        constraint.response = decodeConstraintResponse(optionalString(value, "response", "reject-step", join(path, "response")), join(path, "response"));
        constraint.message = optionalString(value, "message", {}, join(path, "message"));
    }

    void decodeTermination(const json::Value& value, TerminationDefinition& termination, const std::string& path) {
        if (!requireObject(value, path)) {
            return;
        }
        termination.id = decodeIdentifier(value, "id", join(path, "id"));
        const auto* condition = value.get("when");
        if (condition == nullptr) {
            diagnostics_.error("NVP029", "Termination rule requires a when condition.", join(path, "when"), {});
        } else {
            termination.condition = decodeExpression(*condition, join(path, "when"), 0);
        }
        termination.action = decodeActionKind(optionalString(value, "action", "permanent-shutdown", join(path, "action")), join(path, "action"));
        termination.reason = optionalString(value, "reason", {}, join(path, "reason"));
        if (termination.action != ActionKind::reversibleShutdown &&
            termination.action != ActionKind::permanentShutdown) {
            diagnostics_.error(
                "NVP030",
                "Termination action must be reversible-shutdown or permanent-shutdown.",
                join(path, "action"),
                {}
            );
        }
    }

    Expression decodeExpression(const json::Value& value, const std::string& path, std::size_t depth) {
        if (depth > 64) {
            diagnostics_.error("NVP031", "Expression nesting exceeds 64 levels.", path, {});
            return {};
        }
        if (value.isBool()) {
            Expression result;
            result.kind = ExpressionKind::literal;
            result.value = value.asBool() ? 1.0 : 0.0;
            return result;
        }
        if (value.isNumber()) {
            Expression result;
            result.kind = ExpressionKind::literal;
            result.value = value.asNumber();
            return result;
        }
        if (!value.isObject()) {
            diagnostics_.error("NVP032", "Expression must be a number, boolean, or object.", path, {});
            return {};
        }

        if (const auto* literal = value.get("literal")) {
            Expression result;
            result.kind = ExpressionKind::literal;
            if (literal->isNumber()) {
                result.value = literal->asNumber();
            } else if (literal->isObject()) {
                result.value = optionalNumber(*literal, "value", 0.0, join(path, "literal.value"));
                result.unit = optionalString(*literal, "unit", "1", join(path, "literal.unit"));
                validateUnit(result.unit, join(path, "literal.unit"));
            } else {
                diagnostics_.error("NVP033", "literal must be a number or quantity object.", join(path, "literal"), {});
            }
            return result;
        }

        if (auto reference = decodeReferenceExpression(value, path); reference.has_value()) {
            return std::move(*reference);
        }

        if (const auto* operand = value.get("not")) {
            Expression result;
            result.kind = ExpressionKind::logicalNot;
            result.children.push_back(decodeExpression(*operand, join(path, "not"), depth + 1));
            return result;
        }
        if (const auto* operands = value.get("all")) {
            return decodeVariadic(ExpressionKind::logicalAll, *operands, join(path, "all"), depth);
        }
        if (const auto* operands = value.get("any")) {
            return decodeVariadic(ExpressionKind::logicalAny, *operands, join(path, "any"), depth);
        }

        const struct {
            std::string_view key;
            ExpressionKind kind;
        } binaryOperators[] = {
            {"gt", ExpressionKind::greater},
            {"gte", ExpressionKind::greaterEqual},
            {"lt", ExpressionKind::less},
            {"lte", ExpressionKind::lessEqual},
            {"eq", ExpressionKind::equal},
            {"neq", ExpressionKind::notEqual},
            {"add", ExpressionKind::add},
            {"subtract", ExpressionKind::subtract},
            {"multiply", ExpressionKind::multiply},
            {"divide", ExpressionKind::divide},
            {"min", ExpressionKind::minimum},
            {"max", ExpressionKind::maximum}
        };

        for (const auto& operation : binaryOperators) {
            if (const auto* operands = value.get(operation.key)) {
                return decodeBinary(operation.kind, *operands, join(path, operation.key), depth);
            }
        }

        if (const auto* clamp = value.get("clamp")) {
            Expression result;
            result.kind = ExpressionKind::clamp;
            if (!clamp->isObject()) {
                diagnostics_.error("NVP034", "clamp must be an object.", join(path, "clamp"), {});
                return result;
            }
            const auto* source = clamp->get("value");
            const auto* minimum = clamp->get("min");
            const auto* maximum = clamp->get("max");
            if (source == nullptr || minimum == nullptr || maximum == nullptr) {
                diagnostics_.error("NVP035", "clamp requires value, min, and max operands.", join(path, "clamp"), {});
                return result;
            }
            result.children.push_back(decodeExpression(*source, join(path, "clamp.value"), depth + 1));
            result.children.push_back(decodeExpression(*minimum, join(path, "clamp.min"), depth + 1));
            result.children.push_back(decodeExpression(*maximum, join(path, "clamp.max"), depth + 1));
            return result;
        }

        if (const auto* sustained = value.get("sustained")) {
            return decodeTemporal(ExpressionKind::sustained, *sustained, join(path, "sustained"), "for", depth);
        }
        if (const auto* within = value.get("within")) {
            return decodeTemporal(ExpressionKind::within, *within, join(path, "within"), "window", depth);
        }
        if (const auto* rising = value.get("risingEdge")) {
            Expression result;
            result.kind = ExpressionKind::risingEdge;
            result.children.push_back(decodeExpression(*rising, join(path, "risingEdge"), depth + 1));
            return result;
        }
        if (const auto* falling = value.get("fallingEdge")) {
            Expression result;
            result.kind = ExpressionKind::fallingEdge;
            result.children.push_back(decodeExpression(*falling, join(path, "fallingEdge"), depth + 1));
            return result;
        }

        diagnostics_.error(
            "NVP036",
            "Expression object does not contain a recognized operator or reference.",
            path,
            {}
        );
        return {};
    }

    std::optional<Expression> decodeReferenceExpression(const json::Value& value, const std::string& path) {
        const struct {
            std::string_view key;
            ReferenceKind kind;
        } references[] = {
            {"signal", ReferenceKind::input},
            {"input", ReferenceKind::input},
            {"species", ReferenceKind::species},
            {"state", ReferenceKind::state},
            {"parameter", ReferenceKind::parameter},
            {"host", ReferenceKind::hostChannel}
        };

        for (const auto& candidate : references) {
            if (const auto* reference = value.get(candidate.key)) {
                if (!reference->isString()) {
                    diagnostics_.error("NVP037", "Expression reference must be a string.", join(path, candidate.key), {});
                    return Expression{};
                }
                Expression result;
                result.kind = ExpressionKind::reference;
                result.referenceKind = candidate.kind;
                result.reference = std::string(reference->asString());
                if (const auto* unit = value.get("unit")) {
                    if (!unit->isString()) {
                        diagnostics_.error("NVP038", "Expression unit must be a string.", join(path, "unit"), {});
                    } else {
                        result.unit = std::string(unit->asString());
                        validateUnit(result.unit, join(path, "unit"));
                    }
                }
                return result;
            }
        }

        if (const auto* reference = value.get("ref")) {
            if (!reference->isObject()) {
                diagnostics_.error("NVP039", "ref must be an object.", join(path, "ref"), {});
                return Expression{};
            }
            Expression result;
            result.kind = ExpressionKind::reference;
            result.reference = requireString(*reference, "id", join(path, "ref.id"));
            const std::string kind = optionalString(*reference, "kind", "input", join(path, "ref.kind"));
            if (kind == "input" || kind == "signal") result.referenceKind = ReferenceKind::input;
            else if (kind == "species") result.referenceKind = ReferenceKind::species;
            else if (kind == "state") result.referenceKind = ReferenceKind::state;
            else if (kind == "parameter") result.referenceKind = ReferenceKind::parameter;
            else if (kind == "host") result.referenceKind = ReferenceKind::hostChannel;
            else {
                diagnostics_.error("NVP040", "Unknown expression reference kind '" + kind + "'.", join(path, "ref.kind"), {});
            }
            result.unit = optionalString(*reference, "unit", {}, join(path, "ref.unit"));
            if (!result.unit.empty()) validateUnit(result.unit, join(path, "ref.unit"));
            return result;
        }

        if (value.get("time") != nullptr) {
            Expression result;
            result.kind = ExpressionKind::reference;
            result.referenceKind = ReferenceKind::time;
            result.reference = "time";
            result.unit = "s";
            return result;
        }
        return std::nullopt;
    }

    Expression decodeVariadic(ExpressionKind kind,
                              const json::Value& value,
                              const std::string& path,
                              std::size_t depth) {
        Expression result;
        result.kind = kind;
        if (!value.isArray() || value.asArray().empty()) {
            diagnostics_.error("NVP041", "Variadic expression requires a non-empty array.", path, {});
            return result;
        }
        for (std::size_t index = 0; index < value.asArray().size(); ++index) {
            result.children.push_back(decodeExpression(value.asArray()[index], indexPath(path, index), depth + 1));
        }
        return result;
    }

    Expression decodeBinary(ExpressionKind kind,
                            const json::Value& value,
                            const std::string& path,
                            std::size_t depth) {
        Expression result;
        result.kind = kind;

        if (value.isArray()) {
            if (value.asArray().size() != 2) {
                diagnostics_.error("NVP042", "Binary expression array must contain exactly two operands.", path, {});
                return result;
            }
            result.children.push_back(decodeExpression(value.asArray()[0], indexPath(path, 0), depth + 1));
            result.children.push_back(decodeExpression(value.asArray()[1], indexPath(path, 1), depth + 1));
            return result;
        }

        if (!value.isObject()) {
            diagnostics_.error("NVP043", "Binary expression must be an array or object.", path, {});
            return result;
        }

        const auto* left = value.get("left");
        const auto* right = value.get("right");
        if (left != nullptr && right != nullptr) {
            result.children.push_back(decodeExpression(*left, join(path, "left"), depth + 1));
            result.children.push_back(decodeExpression(*right, join(path, "right"), depth + 1));
            return result;
        }

        auto shorthandReference = decodeReferenceExpression(value, path);
        const auto* threshold = value.get("value");
        if (shorthandReference.has_value() && threshold != nullptr) {
            result.children.push_back(std::move(*shorthandReference));
            Expression literal;
            literal.kind = ExpressionKind::literal;
            if (!threshold->isNumber()) {
                diagnostics_.error("NVP044", "Comparison shorthand value must be numeric.", join(path, "value"), {});
            } else {
                literal.value = threshold->asNumber();
            }
            literal.unit = optionalString(value, "unit", result.children.front().unit, join(path, "unit"));
            if (!literal.unit.empty()) validateUnit(literal.unit, join(path, "unit"));
            result.children.push_back(std::move(literal));
            return result;
        }

        diagnostics_.error("NVP045", "Binary expression requires left/right or reference/value operands.", path, {});
        return result;
    }

    Expression decodeTemporal(ExpressionKind kind,
                              const json::Value& value,
                              const std::string& path,
                              std::string_view durationKey,
                              std::size_t depth) {
        Expression result;
        result.kind = kind;
        if (!value.isObject()) {
            diagnostics_.error("NVP046", "Temporal expression must be an object.", path, {});
            return result;
        }
        const auto* condition = value.get("condition");
        const auto* duration = value.get(durationKey);
        if (condition == nullptr || duration == nullptr) {
            diagnostics_.error(
                "NVP047",
                "Temporal expression requires condition and " + std::string(durationKey) + ".",
                path,
                {}
            );
            return result;
        }
        result.children.push_back(decodeExpression(*condition, join(path, "condition"), depth + 1));
        result.durationSeconds = decodeDuration(*duration, join(path, durationKey));
        if (result.durationSeconds <= 0.0) {
            diagnostics_.error("NVP048", "Temporal duration must be positive.", join(path, durationKey), {});
        }
        return result;
    }

    void decodeStoichiometry(const json::Value* value,
                             std::vector<StoichiometryTerm>& output,
                             const std::string& path) {
        if (value == nullptr) {
            return;
        }
        if (!value->isArray()) {
            diagnostics_.error("NVP049", "Stoichiometry must be an array.", path, {});
            return;
        }
        for (std::size_t index = 0; index < value->asArray().size(); ++index) {
            const auto& item = value->asArray()[index];
            const std::string itemPath = indexPath(path, index);
            StoichiometryTerm term;
            if (item.isString()) {
                term.species = std::string(item.asString());
                term.coefficient = 1;
            } else if (item.isObject()) {
                term.species = requireString(item, "species", join(itemPath, "species"));
                const auto coefficient = optionalInteger(item, "coefficient", 1, join(itemPath, "coefficient"));
                if (coefficient <= 0 || coefficient > std::numeric_limits<std::int16_t>::max()) {
                    diagnostics_.error("NVP050", "Stoichiometric coefficient must be in 1...32767.", join(itemPath, "coefficient"), {});
                } else {
                    term.coefficient = static_cast<std::int16_t>(coefficient);
                }
            } else {
                diagnostics_.error("NVP051", "Stoichiometry term must be a species string or object.", itemPath, {});
            }
            output.push_back(std::move(term));
        }
    }

    Bounds decodeBounds(const json::Value* value, const std::string& path) {
        Bounds bounds;
        if (value == nullptr) {
            return bounds;
        }
        if (!value->isObject()) {
            diagnostics_.error("NVP052", "bounds must be an object.", path, {});
            return bounds;
        }
        if (const auto* minimum = value->get("min")) {
            if (!minimum->isNumber() || !std::isfinite(minimum->asNumber())) {
                diagnostics_.error("NVP053", "bounds.min must be finite and numeric.", join(path, "min"), {});
            } else {
                bounds.minimum = minimum->asNumber();
            }
        }
        if (const auto* maximum = value->get("max")) {
            if (!maximum->isNumber() || !std::isfinite(maximum->asNumber())) {
                diagnostics_.error("NVP054", "bounds.max must be finite and numeric.", join(path, "max"), {});
            } else {
                bounds.maximum = maximum->asNumber();
            }
        }
        if (bounds.minimum.has_value() && bounds.maximum.has_value() && *bounds.minimum > *bounds.maximum) {
            diagnostics_.error("NVP055", "bounds.min cannot exceed bounds.max.", path, {});
        }
        return bounds;
    }

    void decodeEvidence(const json::Value* value, Evidence& evidence, const std::string& path) {
        if (value == nullptr) {
            return;
        }
        if (!value->isObject()) {
            diagnostics_.error("NVP056", "evidence must be an object.", path, {});
            return;
        }
        const std::string classification = optionalString(*value, "class", "assumed", join(path, "class"));
        if (classification == "observed") evidence.classification = EvidenceClass::observed;
        else if (classification == "derived") evidence.classification = EvidenceClass::derived;
        else if (classification == "calibrated") evidence.classification = EvidenceClass::calibrated;
        else if (classification == "inferred") evidence.classification = EvidenceClass::inferred;
        else if (classification == "assumed") evidence.classification = EvidenceClass::assumed;
        else if (classification == "hypothetical") evidence.classification = EvidenceClass::hypothetical;
        else diagnostics_.error("NVP057", "Unknown evidence class '" + classification + "'.", join(path, "class"), {});
        evidence.source = optionalString(*value, "source", {}, join(path, "source"));
        evidence.dataset = optionalString(*value, "dataset", {}, join(path, "dataset"));
        evidence.context = optionalString(*value, "context", {}, join(path, "context"));
        evidence.note = optionalString(*value, "note", {}, join(path, "note"));
    }

    double decodeQuantity(const json::Value* value,
                          double fallback,
                          std::string_view expectedUnit,
                          const std::string& path,
                          bool required = false) {
        if (value == nullptr) {
            if (required) {
                diagnostics_.error("NVP058", "Quantity is required.", path, {});
            }
            return fallback;
        }
        if (value->isNumber()) {
            if (!std::isfinite(value->asNumber())) {
                diagnostics_.error("NVP059", "Quantity must be finite.", path, {});
                return fallback;
            }
            return value->asNumber();
        }
        if (!value->isObject()) {
            diagnostics_.error("NVP060", "Quantity must be a number or {value, unit} object.", path, {});
            return fallback;
        }
        const double numeric = optionalNumber(*value, "value", fallback, join(path, "value"));
        const std::string sourceUnit = optionalString(*value, "unit", std::string(expectedUnit), join(path, "unit"));
        validateUnit(sourceUnit, join(path, "unit"));
        const auto converted = units_.convert(numeric, sourceUnit, expectedUnit);
        if (!converted.has_value()) {
            diagnostics_.error(
                "NVP061",
                "Quantity unit '" + sourceUnit + "' is not compatible with '" + std::string(expectedUnit) + "'.",
                path,
                {}
            );
            return fallback;
        }
        return *converted;
    }

    double decodeDuration(const json::Value& value, const std::string& path) {
        if (value.isNumber()) {
            if (!std::isfinite(value.asNumber()) || value.asNumber() < 0.0) {
                diagnostics_.error("NVP062", "Duration in seconds must be finite and non-negative.", path, {});
                return 0.0;
            }
            return value.asNumber();
        }
        if (value.isString()) {
            const auto seconds = units_.parseDurationSeconds(value.asString());
            if (!seconds.has_value()) {
                diagnostics_.error("NVP063", "Invalid duration '" + std::string(value.asString()) + "'.", path, "Use a value such as '20 min', '6 h', or '21 d'.");
                return 0.0;
            }
            return *seconds;
        }
        diagnostics_.error("NVP064", "Duration must be seconds or a duration string.", path, {});
        return 0.0;
    }

    void validateBounds(double value, const Bounds& bounds, const std::string& path) {
        if (!std::isfinite(value)) {
            diagnostics_.error("NVP065", "Value must be finite.", path, {});
            return;
        }
        if (bounds.minimum.has_value() && value < *bounds.minimum) {
            diagnostics_.error("NVP066", "Initial value is below the declared minimum.", path, {});
        }
        if (bounds.maximum.has_value() && value > *bounds.maximum) {
            diagnostics_.error("NVP067", "Initial value exceeds the declared maximum.", path, {});
        }
    }

    void validateUnit(std::string_view unit, const std::string& path) {
        if (units_.find(unit) == nullptr) {
            diagnostics_.error("NVP068", "Unknown unit '" + std::string(unit) + "'.", path, "Declare a supported canonical unit rather than relying on implicit scaling.");
        }
    }

    FidelityLevel decodeFidelity(const json::Value& value, const std::string& path) {
        if (!value.isString()) {
            diagnostics_.error("NVP069", "minimumFidelity must be a string.", path, {});
            return FidelityLevel::f0Logic;
        }
        const auto text = value.asString();
        if (text == "F0") return FidelityLevel::f0Logic;
        if (text == "F1") return FidelityLevel::f1Deterministic;
        if (text == "F2") return FidelityLevel::f2Stochastic;
        if (text == "F3") return FidelityLevel::f3Spatial;
        if (text == "F4") return FidelityLevel::f4Tissue;
        diagnostics_.error("NVP070", "Unknown fidelity level '" + std::string(text) + "'.", path, {});
        return FidelityLevel::f0Logic;
    }

    SignalSource decodeSignalSource(std::string_view value, const std::string& path) {
        if (value == "intracellular") return SignalSource::intracellular;
        if (value == "extracellular") return SignalSource::extracellular;
        if (value == "membrane") return SignalSource::membrane;
        if (value == "host") return SignalSource::host;
        if (value == "numanx") return SignalSource::numanX;
        if (value == "numitissue") return SignalSource::numiTissue;
        if (value == "numibrain") return SignalSource::numiBrain;
        if (value == "experiment") return SignalSource::experiment;
        diagnostics_.error("NVP071", "Unknown signal source '" + std::string(value) + "'.", path, {});
        return SignalSource::extracellular;
    }

    SpeciesKind decodeSpeciesKind(std::string_view value, const std::string& path) {
        if (value == "count") return SpeciesKind::molecularCount;
        if (value == "concentration") return SpeciesKind::concentration;
        if (value == "activity") return SpeciesKind::activity;
        if (value == "occupancy") return SpeciesKind::occupancy;
        if (value == "output") return SpeciesKind::output;
        if (value == "external-field") return SpeciesKind::externalField;
        if (value == "latent") return SpeciesKind::latent;
        diagnostics_.error("NVP072", "Unknown species kind '" + std::string(value) + "'.", path, {});
        return SpeciesKind::concentration;
    }

    StateKind decodeStateKind(std::string_view value, const std::string& path) {
        if (value == "scalar") return StateKind::scalar;
        if (value == "leaky-integrator") return StateKind::leakyIntegrator;
        if (value == "latch") return StateKind::latch;
        if (value == "counter") return StateKind::counter;
        if (value == "timer") return StateKind::timer;
        if (value == "finite-state") return StateKind::finiteState;
        if (value == "permanent-memory") return StateKind::permanentMemory;
        diagnostics_.error("NVP073", "Unknown state type '" + std::string(value) + "'.", path, {});
        return StateKind::scalar;
    }

    RateLawKind decodeRateLaw(std::string_view value, const std::string& path) {
        if (value == "zero-order") return RateLawKind::zeroOrder;
        if (value == "mass-action") return RateLawKind::massAction;
        if (value == "hill-activation") return RateLawKind::hillActivation;
        if (value == "hill-repression") return RateLawKind::hillRepression;
        if (value == "michaelis-menten") return RateLawKind::michaelisMenten;
        if (value == "reversible-mass-action") return RateLawKind::reversibleMassAction;
        if (value == "passive-transport") return RateLawKind::passiveTransport;
        if (value == "saturable-transport") return RateLawKind::saturableTransport;
        if (value == "degradation") return RateLawKind::degradation;
        if (value == "custom") return RateLawKind::customBytecode;
        diagnostics_.error("NVP074", "Unknown rate law '" + std::string(value) + "'.", path, {});
        return RateLawKind::massAction;
    }

    ActionKind decodeActionKind(std::string_view value, const std::string& path) {
        if (value == "set-output") return ActionKind::setOutput;
        if (value == "add-output") return ActionKind::addOutput;
        if (value == "express") return ActionKind::express;
        if (value == "suppress") return ActionKind::suppress;
        if (value == "degrade") return ActionKind::degrade;
        if (value == "set-state") return ActionKind::setState;
        if (value == "increment-state") return ActionKind::incrementState;
        if (value == "emit-event") return ActionKind::emitEvent;
        if (value == "request-differentiation") return ActionKind::requestDifferentiation;
        if (value == "request-migration") return ActionKind::requestMigration;
        if (value == "reversible-shutdown") return ActionKind::reversibleShutdown;
        if (value == "permanent-shutdown") return ActionKind::permanentShutdown;
        diagnostics_.error("NVP075", "Unknown action kind '" + std::string(value) + "'.", path, {});
        return ActionKind::setOutput;
    }

    ConstraintResponse decodeConstraintResponse(std::string_view value, const std::string& path) {
        if (value == "record") return ConstraintResponse::record;
        if (value == "clamp") return ConstraintResponse::clamp;
        if (value == "reject-step") return ConstraintResponse::rejectStep;
        if (value == "substep") return ConstraintResponse::substep;
        if (value == "reversible-shutdown") return ConstraintResponse::reversibleShutdown;
        if (value == "permanent-shutdown") return ConstraintResponse::permanentShutdown;
        diagnostics_.error("NVP076", "Unknown constraint response '" + std::string(value) + "'.", path, {});
        return ConstraintResponse::rejectStep;
    }

    Severity decodeSeverity(std::string_view value, const std::string& path) {
        if (value == "note") return Severity::note;
        if (value == "warning") return Severity::warning;
        if (value == "error") return Severity::error;
        if (value == "fatal") return Severity::fatal;
        diagnostics_.error("NVP077", "Unknown severity '" + std::string(value) + "'.", path, {});
        return Severity::error;
    }

    bool requireObject(const json::Value& value, const std::string& path) {
        if (!value.isObject()) {
            diagnostics_.error("NVP078", "Expected an object.", path, {});
            return false;
        }
        return true;
    }

    std::string decodeIdentifier(const json::Value& object, std::string_view key, const std::string& path) {
        std::string identifier = requireString(object, key, path);
        if (!identifier.empty() && !isValidIdentifier(identifier)) {
            diagnostics_.error("NVP079", "Invalid identifier '" + identifier + "'.", path, {});
        }
        return identifier;
    }

    std::string requireString(const json::Value& object, std::string_view key, const std::string& path) {
        const auto* value = object.get(key);
        if (value == nullptr) {
            diagnostics_.error("NVP080", "Required string is missing.", path, {});
            return {};
        }
        if (!value->isString()) {
            diagnostics_.error("NVP081", "Expected a string.", path, {});
            return {};
        }
        return std::string(value->asString());
    }

    std::string optionalString(const json::Value& object,
                               std::string_view key,
                               std::string fallback,
                               const std::string& path) {
        const auto* value = object.get(key);
        if (value == nullptr) return fallback;
        if (!value->isString()) {
            diagnostics_.error("NVP082", "Expected a string.", path, {});
            return fallback;
        }
        return std::string(value->asString());
    }

    double optionalNumber(const json::Value& object,
                          std::string_view key,
                          double fallback,
                          const std::string& path) {
        const auto* value = object.get(key);
        if (value == nullptr) return fallback;
        if (!value->isNumber() || !std::isfinite(value->asNumber())) {
            diagnostics_.error("NVP083", "Expected a finite number.", path, {});
            return fallback;
        }
        return value->asNumber();
    }

    std::int64_t optionalInteger(const json::Value& object,
                                 std::string_view key,
                                 std::int64_t fallback,
                                 const std::string& path) {
        const auto* value = object.get(key);
        if (value == nullptr) return fallback;
        if (!value->isNumber() || !std::isfinite(value->asNumber()) || std::floor(value->asNumber()) != value->asNumber()) {
            diagnostics_.error("NVP084", "Expected an integer.", path, {});
            return fallback;
        }
        if (value->asNumber() < static_cast<double>(std::numeric_limits<std::int64_t>::min()) ||
            value->asNumber() > static_cast<double>(std::numeric_limits<std::int64_t>::max())) {
            diagnostics_.error("NVP085", "Integer is outside the supported range.", path, {});
            return fallback;
        }
        return static_cast<std::int64_t>(value->asNumber());
    }

    bool optionalBool(const json::Value& object,
                      std::string_view key,
                      bool fallback,
                      const std::string& path) {
        const auto* value = object.get(key);
        if (value == nullptr) return fallback;
        if (!value->isBool()) {
            diagnostics_.error("NVP086", "Expected a boolean.", path, {});
            return fallback;
        }
        return value->asBool();
    }

    template <typename Callback>
    void decodeArray(const json::Value& object,
                     std::string_view key,
                     const std::string& path,
                     Callback&& callback) {
        const auto* value = object.get(key);
        if (value == nullptr) return;
        if (!value->isArray()) {
            diagnostics_.error("NVP087", "Expected an array.", path, {});
            return;
        }
        for (std::size_t index = 0; index < value->asArray().size(); ++index) {
            callback(value->asArray()[index], indexPath(path, index));
        }
    }

    static std::string join(std::string_view base, std::string_view component) {
        if (component.find('.') != std::string_view::npos) {
            return std::string(base) + "." + std::string(component);
        }
        return std::string(base) + "." + std::string(component);
    }

    static std::string indexPath(std::string_view base, std::size_t index) {
        return std::string(base) + "[" + std::to_string(index) + "]";
    }

    const json::Value& root_;
    const UnitRegistry& units_;
    Diagnostics diagnostics_;
};

} // namespace

DecodeResult decodeProgram(const json::Value& root, const UnitRegistry& units) {
    return ProgramDecoder(root, units).run();
}

} // namespace nvivo
