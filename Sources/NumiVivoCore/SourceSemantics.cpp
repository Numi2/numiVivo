#include "SourceSemantics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace nvivo {
namespace {
struct Failure : std::runtime_error { using std::runtime_error::runtime_error; };
void check(bool value, const std::string& message) { if (!value) throw Failure(message); }
struct Unit {
    std::array<int,6> powers{};
    double scale = 1;
};
Unit unit(const UnitRegistry& registry, const std::string& name) {
    const auto* definition = registry.find(name.empty() ? "1" : name);
    check(definition != nullptr, "Unknown source unit '" + name + "'.");
    check(definition->offsetToSI == 0 && std::isfinite(definition->scaleToSI) && definition->scaleToSI > 0,
          "Affine units need explicit conversion before compilation: " + name);
    const auto& d = definition->dimension;
    return {{d.length,d.mass,d.time,d.amount,d.temperature,d.count},definition->scaleToSI};
}
Unit combine(Unit a, const Unit& b, int sign) {
    for (std::size_t i = 0; i < a.powers.size(); ++i) {
        a.powers[i] += sign * b.powers[i];
        check(std::abs(a.powers[i]) <= 127, "Derived unit exponent exceeds its bounded representation.");
    }
    a.scale = sign > 0 ? a.scale * b.scale : a.scale / b.scale;
    check(std::isfinite(a.scale) && a.scale > 0, "Derived unit scale overflows or underflows.");
    return a;
}
bool basis(const Unit& a, const Unit& b) {
    return a.powers == b.powers && std::abs(a.scale / b.scale - 1) <= 1e-10;
}
void representable(double value, const std::string& subject) {
    check(std::isfinite(value) && std::abs(value) <= std::numeric_limits<float>::max(), subject + " exceeds finite FP32 range.");
    check(value == 0 || static_cast<float>(value) != 0, subject + " underflows FP32; no implicit zero substitution is permitted.");
}
void countValue(double value, const std::string& subject) {
    check(std::isfinite(value) && value >= 0 && value <= 16'777'216 && std::trunc(value) == value,
          subject + " is not an exactly representable nonnegative FP32 count; use the UInt32 hybrid runtime for larger populations.");
}
void bounds(const Bounds& value, const std::string& subject, bool runtimeFP32) {
    if (value.minimum) { check(std::isfinite(*value.minimum), subject + " has a nonfinite lower bound."); if(runtimeFP32) representable(*value.minimum,subject); }
    if (value.maximum) { check(std::isfinite(*value.maximum), subject + " has a nonfinite upper bound."); if(runtimeFP32) representable(*value.maximum,subject); }
    check(!value.minimum || !value.maximum || *value.minimum <= *value.maximum, subject + " has inverted bounds.");
}
struct Symbol { ReferenceKind kind; Unit units; bool external; double value; };
class Validator {
    const Program& program;
    UnitRegistry registry;
    std::map<std::string,Symbol,std::less<>> symbols;
    std::size_t expressionNodes = 0;
    Unit timeUnit = unit(registry,"s");
    const Symbol& symbol(const std::string& id) const {
        auto found = symbols.find(id);
        check(found != symbols.end(), "Unresolved source symbol '" + id + "'.");
        return found->second;
    }
    void add(const std::string& id, ReferenceKind kind, const std::string& units, bool external, double value) {
        check(symbols.emplace(id,Symbol{kind,unit(registry,units),external,value}).second, "Duplicate source symbol '" + id + "'.");
    }
    Unit expression(const Expression& e, unsigned depth = 0) {
        check(depth <= 64 && ++expressionNodes <= 1'048'576, "Source expression work bound exceeded.");
        representable(e.value,"Expression literal");
        representable(e.durationSeconds,"Temporal duration");
        if (e.kind == ExpressionKind::literal) return unit(registry,e.unit);
        if (e.kind == ExpressionKind::reference) {
            if (e.referenceKind == ReferenceKind::time) return timeUnit;
            const auto& s = symbol(e.reference);
            check(s.kind == e.referenceKind || (e.referenceKind == ReferenceKind::hostChannel && s.kind == ReferenceKind::input && s.external),
                  "Reference category disagrees with declaration of '" + e.reference + "'.");
            if (e.unit.empty()) return s.units;
            const auto requested = unit(registry,e.unit);
            check(requested.powers == s.units.powers, "Reference unit dimension mismatch.");
            representable(s.units.scale / requested.scale,"Reference unit conversion factor");
            return requested;
        }
        const bool fold = e.kind == ExpressionKind::logicalAll || e.kind == ExpressionKind::logicalAny;
        const bool unary = e.kind == ExpressionKind::logicalNot || e.kind == ExpressionKind::sustained ||
            e.kind == ExpressionKind::within || e.kind == ExpressionKind::risingEdge || e.kind == ExpressionKind::fallingEdge;
        const std::size_t required = e.kind == ExpressionKind::clamp ? 3 : unary ? 1 : 2;
        check(fold ? !e.children.empty() : e.children.size() == required, "Expression operand arity mismatch.");
        std::vector<Unit> children;
        for (const auto& child : e.children) children.push_back(expression(child,depth + 1));
        const Unit scalar{};
        if (fold || unary) {
            for (const auto& child : children) check(child.powers == scalar.powers, "Logical operand has physical units.");
            return scalar;
        }
        if (e.kind == ExpressionKind::multiply) return combine(children[0],children[1],1);
        if (e.kind == ExpressionKind::divide) return combine(children[0],children[1],-1);
        for (const auto& child : children) {
            check(child.powers == children[0].powers, "Expression combines incompatible dimensions.");
            representable(child.scale / children[0].scale,"Expression unit conversion factor");
        }
        if (e.kind >= ExpressionKind::greater && e.kind <= ExpressionKind::notEqual) return scalar;
        return children[0];
    }
    void condition(const Expression& e) { check(expression(e).powers == Unit{}.powers,"Condition is not dimensionless."); }
    void kinetic(const ReactionDefinition& r) {
        std::map<std::string,long long,std::less<>> net;
        auto product = [&](const std::vector<StoichiometryTerm>& terms, int sign) {
            Unit result;
            for (const auto& term : terms) {
                const auto& s = symbol(term.species);
                check(s.kind != ReferenceKind::parameter && term.coefficient > 0 && term.coefficient <= 32,
                      "Reaction species or stoichiometric order is unsupported.");
                net[term.species] += sign * static_cast<long long>(term.coefficient);
                for (int i = 0; i < term.coefficient; ++i) result = combine(result,s.units,1);
            }
            return result;
        };
        const auto reactants = product(r.reactants,-1), products = product(r.products,1);
        std::optional<Unit> extent;
        for (const auto& [id,change] : net) {
            const auto& s = symbol(id);
            if (change == 0 || s.external) continue;
            if (!extent) extent = s.units;
            else check(basis(*extent,s.units), "Reaction '" + r.id + "' mixes numeric species bases; explicitly normalize all changing species units.");
        }
        check(extent.has_value(), "Reaction '" + r.id + "' has no internally owned net state change.");
        auto requireParameter = [&](std::size_t index, Unit expected) {
            check(index < r.rate.parameters.size(), "Missing kinetic parameter in '" + r.id + "'.");
            const auto& p = symbol(r.rate.parameters[index]);
            check(p.kind == ReferenceKind::parameter && p.value >= 0, "Kinetic parameter is not a nonnegative parameter declaration.");
            check(basis(p.units,expected), "Kinetic parameter '" + r.rate.parameters[index] +
                  "' needs an explicit unit/value conversion to the species numeric basis and seconds; runtime parameters are not automatically rescaled.");
        };
        const auto rate = combine(*extent,timeUnit,-1);
        switch (r.rate.law) {
        case RateLawKind::zeroOrder: requireParameter(0,rate); break;
        case RateLawKind::massAction: case RateLawKind::degradation:
            requireParameter(0,combine(rate,reactants,-1)); break;
        case RateLawKind::reversibleMassAction:
            requireParameter(0,combine(rate,reactants,-1)); requireParameter(1,combine(rate,products,-1)); break;
        case RateLawKind::hillActivation: case RateLawKind::hillRepression: case RateLawKind::michaelisMenten:
            check(!r.reactants.empty(), "Saturating kinetics require an input species.");
            requireParameter(0,rate); requireParameter(1,symbol(r.reactants[0].species).units);
            if (r.rate.law != RateLawKind::michaelisMenten) requireParameter(2,Unit{});
            break;
        case RateLawKind::passiveTransport: case RateLawKind::saturableTransport:
            check(!r.reactants.empty() && !r.products.empty(), "Exchange kinetics require both compartments.");
            check(basis(symbol(r.reactants[0].species).units,*extent) && basis(symbol(r.products[0].species).units,*extent),
                  "Exchange source/destination numeric units disagree.");
            if (r.rate.law == RateLawKind::passiveTransport) requireParameter(0,combine(Unit{},timeUnit,-1));
            else { requireParameter(0,rate); requireParameter(1,*extent); }
            break;
        case RateLawKind::customBytecode:
            check(r.rate.expression.has_value(), "Custom kinetics require an expression.");
            check(basis(expression(*r.rate.expression),rate), "Custom rate expression result has incompatible units.");
            break;
        }
        if (r.rate.expression && r.rate.law != RateLawKind::customBytecode)
            check(false,"A built-in kinetic law cannot silently ignore a custom rate expression.");
        if (r.gate) condition(*r.gate);
        representable(r.delaySeconds,"Reaction delay");
    }
public:
    explicit Validator(const Program& p) : program(p) {}
    void run() {
        for (const auto& s : program.inputs) {
            representable(s.defaultValue,s.id); bounds(s.bounds,s.id,true);
            if (unit(registry,s.unit).powers == unit(registry,"count").powers) countValue(s.defaultValue,s.id);
            add(s.id,ReferenceKind::input,s.unit,true,s.defaultValue);
        }
        for (const auto& s : program.species) {
            representable(s.initialValue,s.id); bounds(s.bounds,s.id,true);
            if (s.kind == SpeciesKind::molecularCount) countValue(s.initialValue,s.id);
            add(s.id,ReferenceKind::species,s.unit,s.externallyOwned,s.initialValue);
        }
        for (const auto& s : program.state) {
            check(s.kind == StateKind::scalar || s.kind == StateKind::counter,
                  "State '" + s.id + "' has a behavior not lowered by ProgramPack v1; express it explicitly with reactions/rules instead of losing its semantics.");
            check(s.halfLifeSeconds == 0 && s.states.empty(), "Scalar/counter state has ignored half-life or finite-state labels.");
            representable(s.initialValue,s.id); bounds(s.bounds,s.id,true);
            if (s.kind == StateKind::counter) countValue(s.initialValue,s.id);
            add(s.id,ReferenceKind::state,s.unit,false,s.initialValue);
        }
        for (const auto& p : program.parameters) {
            representable(p.value,p.id); bounds(p.bounds,p.id,false);
            add(p.id,ReferenceKind::parameter,p.unit,false,p.value);
        }
        for (const auto& r : program.reactions) kinetic(r);
        for (const auto& r : program.rules) {
            condition(r.condition); representable(r.refractorySeconds,"Rule refractory duration");
            for (const auto& a : r.actions) {
                representable(a.constantValue,"Action constant"); representable(a.maximumRate,"Action rate limit");
                std::optional<Unit> value;
                if (a.value) value = expression(*a.value);
                if (a.kind >= ActionKind::emitEvent) continue;
                const auto& target = symbol(a.target);
                check(!target.external && target.kind != ReferenceKind::parameter,"Action target is not internally owned state.");
                if (a.kind == ActionKind::setState || a.kind == ActionKind::incrementState)
                    check(target.kind == ReferenceKind::state,"State action targets a non-state symbol.");
                Unit expected = target.units;
                if (a.kind == ActionKind::express || a.kind == ActionKind::suppress || a.kind == ActionKind::degrade)
                    expected = combine(expected,timeUnit,-1);
                check(basis(value.value_or(unit(registry,a.unit)),expected),"Action result/constant unit differs from the target numeric basis.");
                if ((a.kind == ActionKind::addOutput || a.kind == ActionKind::incrementState) && a.maximumRate != 0)
                    check(false,"Add/increment actions do not apply maximumRate; use an explicitly bounded rate action.");
            }
        }
        for (const auto& c : program.constraints) condition(c.condition);
        for (const auto& t : program.termination) condition(t.condition);
    }
};
}
void validateExecutableSource(const Program& program, Diagnostics& diagnostics) {
    try { Validator(program).run(); }
    catch (const Failure& failure) { diagnostics.error("NVC100",failure.what(),"$.spec","Normalize units or use explicitly supported source semantics before compilation."); }
}
}
