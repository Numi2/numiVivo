#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <set>

namespace nvivo {

namespace {

void collectReferences(const Expression& expression,
                       std::set<std::string, std::less<>>& all,
                       std::set<std::string, std::less<>>& inputs) {
    if (expression.kind == ExpressionKind::reference) {
        all.insert(expression.reference);
        if (expression.referenceKind == ReferenceKind::input ||
            expression.referenceKind == ReferenceKind::hostChannel) {
            inputs.insert(expression.reference);
        }
    }
    for (const auto& child : expression.children) {
        collectReferences(child, all, inputs);
    }
}

bool isBiologicalOutput(ActionKind kind) {
    switch (kind) {
        case ActionKind::setOutput:
        case ActionKind::addOutput:
        case ActionKind::express:
        case ActionKind::suppress:
        case ActionKind::degrade:
        case ActionKind::requestDifferentiation:
        case ActionKind::requestMigration:
            return true;
        default:
            return false;
    }
}

void addFinding(SafetyReport& report,
                Diagnostics& diagnostics,
                std::string id,
                FindingCategory category,
                Severity severity,
                std::string subject,
                std::string message,
                std::string mitigation) {
    report.findings.push_back({
        id,
        category,
        severity,
        subject,
        message,
        mitigation
    });

    const std::string diagnosticCode = "NVS" + id;
    if (severity == Severity::fatal) {
        diagnostics.fatal(diagnosticCode, message, subject, mitigation);
        report.hasBlockingFinding = true;
    } else if (severity == Severity::error) {
        diagnostics.error(diagnosticCode, message, subject, mitigation);
        report.hasBlockingFinding = true;
    } else if (severity == Severity::warning) {
        diagnostics.warning(diagnosticCode, message, subject, mitigation);
    } else {
        diagnostics.note(diagnosticCode, message, subject, mitigation);
    }
}

const SpeciesDefinition* findSpecies(const Program& program, std::string_view id) {
    const auto iterator = std::find_if(program.species.begin(), program.species.end(), [id](const auto& item) {
        return item.id == id;
    });
    return iterator == program.species.end() ? nullptr : &*iterator;
}

} // namespace

SafetyReport SafetyAnalyzer::analyze(const Program& program,
                                     const CompileOptions& options,
                                     Diagnostics& diagnostics) const {
    SafetyReport report;
    std::set<std::string, std::less<>> monitoredReferences;
    std::set<std::string, std::less<>> unusedInputs;

    for (const auto& constraint : program.constraints) {
        std::set<std::string, std::less<>> inputs;
        collectReferences(constraint.condition, monitoredReferences, inputs);
    }
    for (const auto& termination : program.termination) {
        std::set<std::string, std::less<>> inputs;
        collectReferences(termination.condition, monitoredReferences, inputs);
    }

    std::set<std::string, std::less<>> outputTargets;
    bool hasPermanentShutdown = false;
    bool hasReversibleShutdown = false;

    for (std::size_t ruleIndex = 0; ruleIndex < program.rules.size(); ++ruleIndex) {
        const auto& rule = program.rules[ruleIndex];
        std::set<std::string, std::less<>> conditionReferences;
        std::set<std::string, std::less<>> conditionInputs;
        collectReferences(rule.condition, conditionReferences, conditionInputs);

        bool ruleHasBiologicalOutput = false;
        for (const auto& action : rule.actions) {
            if (isBiologicalOutput(action.kind)) {
                ruleHasBiologicalOutput = true;
                if (!action.target.empty()) {
                    outputTargets.insert(action.target);
                }
            }
            if (action.kind == ActionKind::permanentShutdown) {
                hasPermanentShutdown = true;
                ++report.irreversibleActionCount;
            }
            if (action.kind == ActionKind::reversibleShutdown) {
                hasReversibleShutdown = true;
            }
            if (action.kind == ActionKind::requestDifferentiation ||
                action.kind == ActionKind::permanentShutdown) {
                ++report.irreversibleActionCount;
            }

            if ((action.kind == ActionKind::addOutput || action.kind == ActionKind::express) &&
                action.maximumRate <= 0.0) {
                const auto* species = findSpecies(program, action.target);
                const bool boundedSpecies = species != nullptr && species->bounds.maximum.has_value();
                if (!boundedSpecies) {
                    addFinding(
                        report,
                        diagnostics,
                        "008-" + std::to_string(ruleIndex),
                        FindingCategory::unboundedAccumulation,
                        options.strictSafety ? Severity::error : Severity::warning,
                        "$.spec.rules[" + std::to_string(ruleIndex) + "]",
                        "Accumulating action '" + action.target + "' has neither a maximum rate nor a bounded target state.",
                        "Declare maximumRate and a finite target bound, then monitor cumulative output."
                    );
                }
            }
        }

        if (ruleHasBiologicalOutput && conditionInputs.size() < 2) {
            addFinding(
                report,
                diagnostics,
                "003-" + std::to_string(ruleIndex),
                FindingCategory::singleSignalActivation,
                Severity::warning,
                "$.spec.rules[" + std::to_string(ruleIndex) + "]",
                "Biological output is gated by fewer than two independent input or host signals.",
                "Document why a single-signal decision is acceptable or add an orthogonal context gate."
            );
        }
    }

    for (std::size_t stateIndex = 0; stateIndex < program.state.size(); ++stateIndex) {
        if (program.state[stateIndex].kind == StateKind::permanentMemory) {
            ++report.irreversibleActionCount;
            addFinding(
                report,
                diagnostics,
                "004-state-" + std::to_string(stateIndex),
                FindingCategory::irreversibleAction,
                Severity::warning,
                "$.spec.state[" + std::to_string(stateIndex) + "]",
                "Program contains persistent molecular memory whose state is not automatically reversible.",
                "Provide a reset, retrieval, containment, and failure-observation strategy for the intended context."
            );
        }
    }

    for (std::size_t terminationIndex = 0; terminationIndex < program.termination.size(); ++terminationIndex) {
        const auto action = program.termination[terminationIndex].action;
        if (action == ActionKind::permanentShutdown) {
            hasPermanentShutdown = true;
        } else if (action == ActionKind::reversibleShutdown) {
            hasReversibleShutdown = true;
        }
    }

    if (options.requireTermination && !hasPermanentShutdown && !hasReversibleShutdown) {
        addFinding(
            report,
            diagnostics,
            "002",
            FindingCategory::missingTermination,
            options.strictSafety ? Severity::error : Severity::warning,
            "$.spec.termination",
            "Program does not define a bounded shutdown path.",
            "Add a termination rule driven by time, abnormal cell state, or an independently controlled signal."
        );
    }

    for (const auto& target : outputTargets) {
        if (monitoredReferences.contains(target)) {
            ++report.monitoredOutputCount;
            continue;
        }
        addFinding(
            report,
            diagnostics,
            "007-" + target,
            FindingCategory::unmonitoredOutput,
            options.strictSafety ? Severity::error : Severity::warning,
            target,
            "Biological output '" + target + "' is not referenced by a constraint or termination monitor.",
            "Add a direct output, rate, cumulative-exposure, or downstream-state monitor."
        );
    }

    for (std::size_t parameterIndex = 0; parameterIndex < program.parameters.size(); ++parameterIndex) {
        const auto& parameter = program.parameters[parameterIndex];
        if (parameter.evidence.classification == EvidenceClass::hypothetical) {
            addFinding(
                report,
                diagnostics,
                "005-" + std::to_string(parameterIndex),
                FindingCategory::hypotheticalParameter,
                options.permitHypotheticalParameters ? Severity::note : Severity::error,
                "$.spec.parameters[" + std::to_string(parameterIndex) + "]",
                "Parameter '" + parameter.id + "' is explicitly hypothetical.",
                "Retain this classification in result certificates and calibrate before higher-fidelity claims."
            );
        }
        if (!parameter.bounds.minimum.has_value() || !parameter.bounds.maximum.has_value()) {
            addFinding(
                report,
                diagnostics,
                "001-param-" + std::to_string(parameterIndex),
                FindingCategory::missingBound,
                Severity::warning,
                "$.spec.parameters[" + std::to_string(parameterIndex) + "]",
                "Parameter '" + parameter.id + "' does not have a complete finite uncertainty interval.",
                "Declare lower and upper bounds or explicitly record why the parameter is fixed."
            );
        }
    }

    if (program.target.species.empty() || program.target.deliveryMode.empty()) {
        addFinding(
            report,
            diagnostics,
            "004-context",
            FindingCategory::unsupportedContext,
            Severity::warning,
            "$.spec.target",
            "Target context omits species or deliveryMode, limiting interpretation beyond abstract circuit behavior.",
            "Declare the host species and delivery mode before F2-F4 contextual simulation."
        );
    }

    if (report.irreversibleActionCount > 0 && !hasPermanentShutdown) {
        addFinding(
            report,
            diagnostics,
            "004-irreversible",
            FindingCategory::irreversibleAction,
            options.strictSafety ? Severity::error : Severity::warning,
            "$.spec",
            "Program contains irreversible state transitions without a declared permanent shutdown path.",
            "Add an independent terminal containment mechanism and model its failure probability."
        );
    }

    return report;
}

} // namespace nvivo
