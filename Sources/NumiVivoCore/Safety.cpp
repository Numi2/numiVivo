#include "NumiVivoCore/Core.hpp"
#include "SourceSemantics.hpp"
#include <algorithm>
#include <set>

namespace nvivo {
namespace {
void collect(const Expression& expression, std::set<std::string,std::less<>>& all,
             std::set<std::string,std::less<>>& inputs) {
    if (expression.kind == ExpressionKind::reference) {
        all.insert(expression.reference);
        if (expression.referenceKind == ReferenceKind::input || expression.referenceKind == ReferenceKind::hostChannel)
            inputs.insert(expression.reference);
    }
    for (const auto& child : expression.children) collect(child,all,inputs);
}
bool output(ActionKind kind) {
    return kind <= ActionKind::degrade || kind == ActionKind::requestDifferentiation || kind == ActionKind::requestMigration;
}
}
SafetyReport SafetyAnalyzer::analyze(const Program& program, const CompileOptions& options, Diagnostics& diagnostics) const {
    SafetyReport report;
    // Compiler::compile invokes this before any floating-point narrowing or
    // information-losing lowering. Disabling risk heuristics cannot disable
    // physical units, reference categories or source-representability checks.
    validateExecutableSource(program,diagnostics);
    if (diagnostics.hasErrors()) return report;
    auto finding = [&](std::string id, FindingCategory category, Severity severity,
                       std::string subject, std::string message, std::string mitigation) {
        report.findings.push_back({id,category,severity,subject,message,mitigation});
        diagnostics.add({severity,"NVS" + id,message,subject,mitigation,{}});
        report.hasBlockingFinding |= severity == Severity::error || severity == Severity::fatal;
    };
    const auto requiredSeverity = options.strictSafety ? Severity::error : Severity::warning;
    std::set<std::string,std::less<>> monitored, outputs;
    for (const auto& constraint : program.constraints) {
        std::set<std::string,std::less<>> inputs;
        collect(constraint.condition,monitored,inputs);
    }
    for (const auto& termination : program.termination) {
        std::set<std::string,std::less<>> inputs;
        collect(termination.condition,monitored,inputs);
    }
    bool permanent = false, reversible = false;
    for (std::size_t i = 0; i < program.rules.size(); ++i) {
        const auto& rule = program.rules[i];
        const std::string path = "$.spec.rules[" + std::to_string(i) + "]";
        std::set<std::string,std::less<>> references, inputs;
        collect(rule.condition,references,inputs);
        bool biologicalOutput = false;
        for (const auto& action : rule.actions) {
            if (output(action.kind)) { biologicalOutput = true; outputs.insert(action.target); }
            if (action.kind == ActionKind::permanentShutdown) permanent = true;
            if (action.kind == ActionKind::reversibleShutdown) reversible = true;
            if (action.kind == ActionKind::permanentShutdown || action.kind == ActionKind::requestDifferentiation)
                ++report.irreversibleActionCount;
            if ((action.kind == ActionKind::addOutput || action.kind == ActionKind::express) && action.maximumRate <= 0) {
                const auto found = std::find_if(program.species.begin(),program.species.end(),[&](const auto& s){ return s.id == action.target; });
                if (found == program.species.end() || !found->bounds.maximum)
                    finding("008-" + std::to_string(i),FindingCategory::unboundedAccumulation,requiredSeverity,path,
                            "Accumulating output '" + action.target + "' has no rate or target upper bound.",
                            "Declare a finite bound and monitor cumulative exposure explicitly.");
            }
        }
        if (biologicalOutput && inputs.size() < 2)
            finding("003-" + std::to_string(i),FindingCategory::singleSignalActivation,Severity::warning,path,
                    "Output is gated by fewer than two distinct named input/host references.",
                    "Document context specificity. Distinct identifiers alone do not establish biological independence.");
        if (rule.condition.kind == ExpressionKind::literal && rule.condition.value <= 0.5)
            finding("010-" + std::to_string(i),FindingCategory::unreachableRule,Severity::warning,path,
                    "Rule condition is a constant false value.","Remove the rule or provide an executable condition.");
    }
    for (const auto& termination : program.termination) {
        if (termination.action == ActionKind::permanentShutdown) permanent = true;
        if (termination.action == ActionKind::reversibleShutdown) reversible = true;
        if (termination.condition.kind == ExpressionKind::literal && termination.condition.value <= 0.5)
            finding("010-" + termination.id,FindingCategory::unreachableRule,Severity::warning,"$.spec.termination",
                    "A declared termination condition is constant false.",
                    "A syntactically present shutdown rule is not a reachable or bounded shutdown guarantee.");
    }
    if (options.requireTermination && !permanent && !reversible)
        finding("002",FindingCategory::missingTermination,requiredSeverity,"$.spec.termination",
                "No shutdown action is declared.","Add a reachable termination rule and verify its timing independently.");
    for (const auto& id : outputs) {
        if (monitored.contains(id)) { ++report.monitoredOutputCount; continue; }
        finding("007-" + id,FindingCategory::unmonitoredOutput,requiredSeverity,id,
                "Output is not referenced by any constraint or termination expression.",
                "Add an appropriate output monitor. Reference coverage alone does not prove a bound.");
    }
    for (std::size_t i = 0; i < program.parameters.size(); ++i) {
        const auto& parameter = program.parameters[i];
        const std::string path = "$.spec.parameters[" + std::to_string(i) + "]";
        if (parameter.evidence.classification == EvidenceClass::hypothetical)
            finding("005-" + std::to_string(i),FindingCategory::hypotheticalParameter,
                    options.permitHypotheticalParameters ? Severity::note : Severity::error,path,
                    "Parameter '" + parameter.id + "' is hypothetical.","Retain its evidence class in every derived result.");
        if (!parameter.bounds.minimum || !parameter.bounds.maximum)
            finding("001-param-" + std::to_string(i),FindingCategory::missingBound,Severity::warning,path,
                    "Parameter lacks a complete finite uncertainty interval.","Declare a justified interval or identify it explicitly as fixed.");
    }
    if (program.target.species.empty() || program.target.deliveryMode.empty())
        finding("004-context",FindingCategory::unsupportedContext,Severity::warning,"$.spec.target",
                "Species or delivery context is missing.","Limit conclusions to the stated abstract model until context is specified.");
    if (report.irreversibleActionCount > 0 && !permanent)
        finding("004-irreversible",FindingCategory::irreversibleAction,requiredSeverity,"$.spec",
                "Irreversible actions are declared without permanent shutdown.",
                "Specify containment, observability and a terminal control mechanism; these heuristics do not establish physical safety.");
    return report;
}
}
