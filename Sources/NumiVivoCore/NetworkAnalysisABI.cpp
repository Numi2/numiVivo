#include "NumiVivoCore/NumiVivoNetworkAnalysis.h"
#include "NumiVivoCore/Core.hpp"
#include "NumiVivoCore/ReactionNetworkAnalysis.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <span>
#include <string>
#include <string_view>

namespace {

using namespace nvivo;

void clearBuffer(NVivoByteBuffer* buffer) {
    if (buffer == nullptr) return;
    buffer->data = nullptr;
    buffer->size = 0;
}

bool copyBuffer(std::string_view value, NVivoByteBuffer* output) {
    if (output == nullptr) return true;
    clearBuffer(output);
    if (value.empty()) return true;
    void* memory = std::malloc(value.size());
    if (memory == nullptr) return false;
    std::memcpy(memory, value.data(), value.size());
    output->data = static_cast<std::uint8_t*>(memory);
    output->size = value.size();
    return true;
}

bool decodeCompileOptions(const NVivoCompileOptions* source,
                          CompileOptions& destination,
                          Diagnostics& diagnostics) {
    NVivoCompileOptions defaults{};
    if (source == nullptr) {
        nvivo_default_compile_options(&defaults);
        source = &defaults;
    }
    if (source->struct_size < sizeof(NVivoCompileOptions) ||
        source->abi_version != NVIVO_COMPILER_ABI_VERSION ||
        source->requested_fidelity > NVIVO_FIDELITY_F4_TISSUE) {
        diagnostics.error("NVNA001", "Compile options are incompatible with the network-analysis ABI.");
        return false;
    }
    destination.requestedFidelity = static_cast<FidelityLevel>(source->requested_fidelity);
    destination.strictUnits = (source->flags & NVIVO_COMPILE_STRICT_UNITS) != 0;
    destination.strictSafety = (source->flags & NVIVO_COMPILE_STRICT_SAFETY) != 0;
    destination.deterministicPack = (source->flags & NVIVO_COMPILE_DETERMINISTIC_PACK) != 0;
    destination.permitHypotheticalParameters =
        (source->flags & NVIVO_COMPILE_PERMIT_HYPOTHETICAL_PARAMETERS) != 0;
    destination.requireTermination = (source->flags & NVIVO_COMPILE_REQUIRE_TERMINATION) != 0;
    destination.limits.maximumSpecies = source->limits.maximum_species;
    destination.limits.maximumParameters = source->limits.maximum_parameters;
    destination.limits.maximumReactions = source->limits.maximum_reactions;
    destination.limits.maximumRules = source->limits.maximum_rules;
    destination.limits.maximumConstraints = source->limits.maximum_constraints;
    destination.limits.maximumExpressionInstructions = source->limits.maximum_expression_instructions;
    destination.limits.maximumStoichiometryTerms = source->limits.maximum_stoichiometry_terms;
    destination.limits.maximumTemporalStates = source->limits.maximum_temporal_states;
    if (destination.limits.maximumSpecies == 0 ||
        destination.limits.maximumParameters == 0 ||
        destination.limits.maximumReactions == 0 ||
        destination.limits.maximumExpressionInstructions == 0 ||
        destination.limits.maximumStoichiometryTerms == 0 ||
        destination.limits.maximumTemporalStates == 0) {
        diagnostics.error("NVNA002", "Compile resource limits must be non-zero.");
        return false;
    }
    return true;
}

bool decodeAnalysisOptions(const NVivoNetworkAnalysisOptions* source,
                           network::Options& destination,
                           Diagnostics& diagnostics) {
    NVivoNetworkAnalysisOptions defaults{};
    if (source == nullptr) {
        nvivo_default_network_analysis_options(&defaults);
        source = &defaults;
    }
    if (source->struct_size < sizeof(NVivoNetworkAnalysisOptions) ||
        source->abi_version != NVIVO_NETWORK_ANALYSIS_ABI_VERSION ||
        source->maximum_dense_species == 0 ||
        source->maximum_dense_reactions == 0 ||
        source->maximum_dense_elements == 0 ||
        source->maximum_conservation_laws == 0 ||
        !std::isfinite(source->rank_tolerance) ||
        source->rank_tolerance <= 0.0 ||
        !std::isfinite(source->conservation_tolerance) ||
        source->conservation_tolerance <= 0.0) {
        diagnostics.error("NVNA003", "Network-analysis options are invalid.");
        return false;
    }
    destination.rankTolerance = source->rank_tolerance;
    destination.conservationTolerance = source->conservation_tolerance;
    destination.maximumDenseSpecies = source->maximum_dense_species;
    destination.maximumDenseReactions = source->maximum_dense_reactions;
    destination.maximumDenseElements = source->maximum_dense_elements;
    destination.maximumConservationLaws = source->maximum_conservation_laws;
    destination.inferConservationLaws =
        (source->flags & NVIVO_NETWORK_INFER_CONSERVATION) != 0;
    destination.computeComplexDeficiency =
        (source->flags & NVIVO_NETWORK_COMPUTE_COMPLEX_DEFICIENCY) != 0;
    return true;
}

NVivoStatus publishDiagnostics(const Diagnostics& diagnostics,
                               NVivoByteBuffer* output,
                               NVivoStatus status) {
    const std::string report = diagnostics.toJson();
    return copyBuffer(report, output) ? status : NVIVO_STATUS_OUT_OF_MEMORY;
}

} // namespace

extern "C" {

void nvivo_default_network_analysis_options(
    NVivoNetworkAnalysisOptions* options) {
    if (options == nullptr) return;
    std::memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->abi_version = NVIVO_NETWORK_ANALYSIS_ABI_VERSION;
    options->flags = NVIVO_NETWORK_INFER_CONSERVATION |
                     NVIVO_NETWORK_COMPUTE_COMPLEX_DEFICIENCY;
    options->maximum_dense_species = 2'048;
    options->maximum_dense_reactions = 2'048;
    options->maximum_conservation_laws = 256;
    options->maximum_dense_elements = 16'777'216;
    options->rank_tolerance = 1.0e-10;
    options->conservation_tolerance = 1.0e-9;
}

NVivoStatus nvivo_analyze_program_json(
    const std::uint8_t* source,
    std::size_t sourceSize,
    const NVivoCompileOptions* compileOptions,
    const NVivoNetworkAnalysisOptions* analysisOptions,
    NVivoByteBuffer* analysisJson,
    NVivoByteBuffer* diagnosticsJson) {
    clearBuffer(analysisJson);
    clearBuffer(diagnosticsJson);
    if (analysisJson == nullptr || (source == nullptr && sourceSize != 0)) {
        return NVIVO_STATUS_INVALID_ARGUMENT;
    }

    try {
        Diagnostics diagnostics;
        CompileOptions compilerConfiguration;
        network::Options analysisConfiguration;
        if (!decodeCompileOptions(compileOptions, compilerConfiguration, diagnostics) ||
            !decodeAnalysisOptions(analysisOptions, analysisConfiguration, diagnostics)) {
            return publishDiagnostics(
                diagnostics,
                diagnosticsJson,
                NVIVO_STATUS_INVALID_ARGUMENT
            );
        }

        const std::string_view sourceText(
            source == nullptr ? "" : reinterpret_cast<const char*>(source),
            sourceSize
        );
        auto parsed = json::parse(sourceText);
        diagnostics.append(parsed.diagnostics);
        if (!parsed.root.has_value()) {
            return publishDiagnostics(
                diagnostics,
                diagnosticsJson,
                NVIVO_STATUS_PARSE_ERROR
            );
        }

        UnitRegistry units;
        auto decoded = decodeProgram(*parsed.root, units);
        diagnostics.append(decoded.diagnostics);
        if (!decoded.program.has_value()) {
            return publishDiagnostics(
                diagnostics,
                diagnosticsJson,
                NVIVO_STATUS_VALIDATION_ERROR
            );
        }

        const auto sourceBytes = std::span<const std::byte>(
            reinterpret_cast<const std::byte*>(sourceText.data()),
            sourceText.size()
        );
        Compiler compiler(std::move(units));
        auto compiled = compiler.compile(
            *decoded.program,
            sourceBytes,
            compilerConfiguration
        );
        diagnostics.append(compiled.diagnostics);
        if (!compiled.ir.has_value()) {
            const NVivoStatus status = compiled.safety.hasBlockingFinding
                ? NVIVO_STATUS_SAFETY_REJECTED
                : NVIVO_STATUS_COMPILE_ERROR;
            return publishDiagnostics(diagnostics, diagnosticsJson, status);
        }

        auto analysis = network::Analyzer().analyze(
            *compiled.ir,
            analysisConfiguration
        );
        diagnostics.append(analysis.diagnostics);
        const std::string analysisText = analysis.toJson();
        const std::string diagnosticsText = diagnostics.toJson();
        if (!copyBuffer(analysisText, analysisJson)) {
            return NVIVO_STATUS_OUT_OF_MEMORY;
        }
        if (!copyBuffer(diagnosticsText, diagnosticsJson)) {
            nvivo_buffer_release(analysisJson);
            return NVIVO_STATUS_OUT_OF_MEMORY;
        }
        return diagnostics.hasErrors()
            ? NVIVO_STATUS_VALIDATION_ERROR
            : NVIVO_STATUS_OK;
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& error) {
        const std::string report = std::string("{\"status\":\"internal-error\",\"message\":\"") +
                                   json::escape(error.what()) + "\"}";
        return copyBuffer(report, diagnosticsJson)
            ? NVIVO_STATUS_INTERNAL_ERROR
            : NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (...) {
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}

} // extern "C"
