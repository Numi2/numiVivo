#include "NumiVivoCore/Core.hpp"
#include "SourceSemantics.hpp"

#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace nvivo;
namespace {
unsigned checks = 0;
void expect(bool value, const std::string& label) {
    if (!value) throw std::runtime_error(label);
    ++checks;
}
template<class T> T load(const Bytes& bytes, std::size_t offset) {
    if (offset > bytes.size() || sizeof(T) > bytes.size() - offset) throw std::runtime_error("fixture read bounds");
    T result{}; std::memcpy(&result,bytes.data() + offset,sizeof(T)); return result;
}
template<class T> void store(Bytes& bytes, std::size_t offset, const T& value) {
    if (offset > bytes.size() || sizeof(T) > bytes.size() - offset) throw std::runtime_error("fixture write bounds");
    std::memcpy(bytes.data() + offset,&value,sizeof(T));
}
PackSectionDescriptor section(const Bytes& bytes, PackSectionType kind) {
    const auto header = load<PackHeader>(bytes,0);
    for (unsigned i = 0; i < header.sectionCount; ++i) {
        const auto s = load<PackSectionDescriptor>(bytes,sizeof(PackHeader) + i * sizeof(PackSectionDescriptor));
        if (s.type == static_cast<unsigned>(kind)) return s;
    }
    throw std::runtime_error("fixture section absent");
}
void rehash(Bytes& bytes) {
    auto header = load<PackHeader>(bytes,0);
    for (unsigned i = 0; i < header.sectionCount; ++i) {
        const auto offset = sizeof(PackHeader) + i * sizeof(PackSectionDescriptor);
        auto s = load<PackSectionDescriptor>(bytes,offset);
        if (s.offset <= bytes.size() && s.size <= bytes.size() - s.offset)
            s.fingerprint = sha256(std::span<const std::byte>(bytes).subspan(s.offset,s.size));
        store(bytes,offset,s);
    }
    header.contentFingerprint = sha256(std::span<const std::byte>(bytes).subspan(sizeof(PackHeader)));
    store(bytes,0,header);
}
ProgramIR fixture() {
    ProgramIR ir;
    ir.fidelity = FidelityLevel::f1Deterministic;
    ir.featureFlags = 1;
    ir.sourceFingerprint = sha256({});
    ir.strings.push_back('\0');
    auto text = [&](const std::string& value) {
        auto offset = static_cast<std::uint32_t>(ir.strings.size());
        ir.strings.insert(ir.strings.end(),value.begin(),value.end());
        ir.strings.push_back('\0');
        return offset;
    };
    const auto x = text("X"), cell = text("cell"), count = text("count"), k = text("k"),
               inverseTime = text("1/s"), reaction = text("decay"), monitor = text("bound"),
               rule = text("noop"), message = text("synthetic fixture"), scalar = text("1");
    // A valid alternate unit string is deliberately present for a mutation that
    // changes kinetics without changing the length or validity of the string table.
    text("1/min");
    ir.species.push_back({x,cell,count,speciesCountValued,10,0,100,0});
    ir.parameters.push_back({k,inverseTime,0,0,0.1,0,1,4,0});
    ir.reactionParameterIndices = {0};
    ir.stoichiometry.push_back({0,1,0});
    ReactionRecord r{};
    r.nameOffset = reaction; r.compartmentOffset = cell;
    r.reactantOffset = 0; r.reactantCount = 1; r.productOffset = 1;
    r.parameterCount = 1; r.rateLaw = 8; r.flags = reactionStochasticEligible;
    r.characteristicRate = 0.1f; r.reserved = UINT32_MAX;
    ir.reactions.push_back(r);
    ir.expressions = {
        {0,0,0,1,0},{255,0,0,0,0}, // rule: true
        {0,0,0,10,0},{255,0,0,0,0}, // set X = 10
        {0,0,0,1,0},{255,0,0,0,0}   // monitor: true
    };
    ir.maximumExpressionStack = 1;
    ActionRecord a{};
    a.kind = 0; a.expressionOffset = 2; a.expressionCount = 2; a.unitOffset = count;
    ir.actions.push_back(a);
    RuleRecord q{};
    q.nameOffset = rule; q.conditionCount = 2; q.actionCount = 1;
    ir.rules.push_back(q);
    MonitorRecord m{};
    m.nameOffset = monitor; m.expressionOffset = 4; m.expressionCount = 2;
    m.messageOffset = message; m.severity = 2; m.response = 2;
    ir.monitors.push_back(m);
    ir.cohorts.push_back({0,1,8,0,1,0.1f,256,0});
    ir.speciesIncidenceOffsets = {0,1};
    ir.speciesIncidence.push_back({0,-1,0,0,0});
    text("{\"program\":{\"sourceFingerprint\":\"" + hexFingerprint(ir.sourceFingerprint) + "\"},\"fidelity\":\"F1\"}");
    (void)scalar;
    return ir;
}
Program source() {
    Program p;
    p.apiVersion = "numivivo.org/v1alpha1"; p.kind = "VivoProgram";
    p.metadata.name = "audit";
    SpeciesDefinition s;
    s.id = "X"; s.kind = SpeciesKind::molecularCount; s.unit = "count"; s.initialValue = 10;
    s.bounds = {0,100}; p.species.push_back(s);
    ParameterDefinition k;
    k.id = "k"; k.unit = "1/s"; k.value = 0.1; k.bounds = {0,1}; p.parameters.push_back(k);
    ReactionDefinition r;
    r.id = "decay"; r.reactants.push_back({"X",1}); r.rate.law = RateLawKind::degradation; r.rate.parameters = {"k"};
    p.reactions.push_back(r);
    return p;
}
}
int main() {
    try {
        const auto built = serializeProgramPack(fixture());
        if (built.bytes.empty()) throw std::runtime_error("valid fixture did not serialize: " + built.diagnostics.toJson());
        const auto& baseline = built.bytes;
        expect(inspectProgramPack(baseline).valid,"baseline inspection");
        expect(inspectProgramPack(baseline,false).valid,"baseline inspection without redundant section hashing");
        const auto emptyHash = hexFingerprint(sha256({}));
        expect(emptyHash == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","SHA256 empty vector");
        auto mutation = [&](const std::string& label, const std::function<void(Bytes&)>& modify, bool recompute = true) {
            Bytes bytes = baseline;
            modify(bytes);
            if (recompute) rehash(bytes);
            const auto inspection = inspectProgramPack(bytes);
            expect(!inspection.valid && inspection.diagnostics.hasErrors(),label);
        };
        auto word = [&](PackSectionType kind, std::size_t offset, std::uint32_t value) {
            return [=](Bytes& bytes) { store(bytes,section(bytes,kind).offset + offset,value); };
        };
        using S = PackSectionType;
        mutation("reactant species index",word(S::stoichiometry,0,42));
        mutation("parameter index",word(S::reactionParameterIndices,0,9));
        mutation("reaction range overflow",word(S::reactions,8,UINT32_MAX));
        mutation("reaction parameter arity",word(S::reactions,28,0));
        mutation("reaction cohort index",word(S::reactions,56,6));
        mutation("ungated sentinel",word(S::reactions,60,0));
        mutation("invalid rate law",word(S::reactions,40,99));
        mutation("cohort class mismatch",word(S::cohorts,8,1));
        mutation("CSR coefficient mismatch",[](Bytes& b) { store(b,section(b,S::speciesIncidence).offset + 4,std::int16_t(-2)); });
        mutation("CSR terminal out of bounds",word(S::speciesIncidenceOffsets,4,2));
        mutation("CSR start nonzero",word(S::speciesIncidenceOffsets,0,1));
        mutation("CSR reaction out of bounds",word(S::speciesIncidence,0,3));
        mutation("bytecode stack underflow",[](Bytes& b) { store(b,section(b,S::expressions).offset,std::uint16_t(12)); });
        mutation("bytecode unknown opcode",[](Bytes& b) { store(b,section(b,S::expressions).offset,std::uint16_t(99)); });
        mutation("bytecode bad species operand",[](Bytes& b) {
            const auto offset = section(b,S::expressions).offset;
            store(b,offset,std::uint16_t(1)); store(b,offset + 4,std::uint32_t(77));
        });
        mutation("bytecode missing terminator",[](Bytes& b) { store(b,section(b,S::expressions).offset + 16,std::uint16_t(3)); });
        mutation("bytecode boundary reference",word(S::rules,4,1));
        mutation("bytecode declared stack mismatch",word(S::runtimeContract,28,9));
        mutation("action out of range",word(S::actions,0,2));
        mutation("unreferenced action",word(S::rules,16,0));
        mutation("monitor invalid response",word(S::monitors,20,10));
        mutation("species nonfinite",word(S::species,16,0x7fc00000));
        mutation("species fractional count",[](Bytes& b) { store(b,section(b,S::species).offset + 16,0.5f); });
        mutation("species bad bound",[](Bytes& b) { store(b,section(b,S::species).offset + 24,5.0f); });
        mutation("parameter underflow",[](Bytes& b) { store(b,section(b,S::parameters).offset + 16,1e-100); });
        mutation("parameter nonfinite",[](Bytes& b) { store(b,section(b,S::parameters).offset + 16,std::numeric_limits<double>::quiet_NaN()); });
        mutation("string interior reference",word(S::species,0,2));
        mutation("invalid UTF8",[](Bytes& b) { b[section(b,S::strings).offset + 1] = std::byte{0xff}; });
        mutation("source header not bound",[](Bytes& b) { b[40] ^= std::byte{1}; });
        mutation("fidelity header not bound",[](Bytes& b) { store(b,24,std::uint32_t(2)); });
        mutation("flags header not bound",[](Bytes& b) { store(b,20,std::uint32_t(3)); });
        mutation("unsafe record alignment",[](Bytes& b) {
            // Descriptor 1 is species.
            store(b,sizeof(PackHeader) + sizeof(PackSectionDescriptor) + 32,std::uint32_t(1));
        });
        mutation("descriptor zero alignment",[](Bytes& b) { store(b,sizeof(PackHeader) + 32,std::uint32_t(0)); });
        mutation("unknown required table",[](Bytes& b) { store(b,sizeof(PackHeader),std::uint32_t(77)); });
        mutation("unknown header reserved",[](Bytes& b) { b[104] = std::byte{1}; });
        mutation("section hash mismatch",[](Bytes& b) { b[section(b,S::species).offset + 16] ^= std::byte{1}; },false);
        mutation("cached noncanonical rate units",[](Bytes& b) {
            const auto strings = section(b,S::strings);
            const auto* raw = reinterpret_cast<const char*>(b.data() + strings.offset);
            std::string_view text(raw,static_cast<std::size_t>(strings.size));
            const auto offset = text.find("1/min");
            if (offset == std::string_view::npos) throw std::runtime_error("fixture unit absent");
            store(b,section(b,S::parameters).offset + 4,static_cast<std::uint32_t>(offset));
        });
        Bytes shortPack(7);
        expect(!inspectProgramPack(shortPack).valid,"truncated pack");
        auto sourceCheck = [&](const std::string& label,const std::function<void(Program&)>& mutate,bool accepted) {
            auto program = source(); mutate(program); Diagnostics diagnostics;
            validateExecutableSource(program,diagnostics);
            expect(!diagnostics.hasErrors() == accepted,label);
        };
        sourceCheck("valid source",[](Program&){},true);
        sourceCheck("source per-minute scale",[](Program& p){ p.parameters[0].unit = "1/min"; },false);
        sourceCheck("source FP32 overflow",[](Program& p){ p.parameters[0].value = 1e100; },false);
        sourceCheck("source FP32 underflow",[](Program& p){ p.parameters[0].value = 1e-100; },false);
        sourceCheck("ignored leaky-integrator behavior",[](Program& p){ StateDefinition s; s.id="memory"; s.kind=StateKind::leakyIntegrator; s.halfLifeSeconds=30; p.state.push_back(s); },false);
        sourceCheck("mixed reaction bases",[](Program& p){ p.species[0].unit="nM"; auto s=p.species[0];s.id="Y";s.unit="uM";p.species.push_back(s);p.reactions[0].products.push_back({"Y",1});},false);
        sourceCheck("reference category mismatch",[](Program& p){ RuleDefinition r;r.id="r";r.condition.kind=ExpressionKind::reference;r.condition.reference="X";r.condition.referenceKind=ReferenceKind::parameter;p.rules.push_back(r); },false);
        sourceCheck("ignored built-in custom expression",[](Program& p){p.reactions[0].rate.expression=Expression{};},false);
        sourceCheck("negative kinetic parameter",[](Program& p){p.parameters[0].value=-1;},false);
        expect(inspectProgramPack(baseline).valid,"baseline preserved after mutations");
        std::cout << "Native pack/semantic checks passed: " << checks << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Native audit check failed: " << error.what() << "\n";
        return 1;
    }
}
