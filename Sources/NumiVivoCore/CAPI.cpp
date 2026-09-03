#include "NumiVivoCore/NumiVivoCore.h"
#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <optional>
#include <sstream>
#include <string_view>

namespace {

using namespace nvivo;

struct PipelineResult {
    std::optional<Program> program;
    std::optional<ProgramIR> ir;
    SafetyReport safety;
    Diagnostics diagnostics;
    NVivoStatus status = NVIVO_STATUS_OK;
};

void clearBuffer(NVivoByteBuffer* buffer) {
    if (buffer == nullptr) return;
    buffer->data = nullptr;
    buffer->size = 0;
}

bool copyBuffer(std::span<const std::byte> bytes, NVivoByteBuffer* output) {
    if (output == nullptr) return true;
    clearBuffer(output);
    if (bytes.empty()) return true;
    void* memory = std::malloc(bytes.size());
    if (memory == nullptr) return false;
    std::memcpy(memory, bytes.data(), bytes.size());
    output->data = static_cast<std::uint8_t*>(memory);
    output->size = bytes.size();
    return true;
}

bool copyBuffer(std::string_view string, NVivoByteBuffer* output) {
    return copyBuffer(
        {reinterpret_cast<const std::byte*>(string.data()), string.size()},
        output
    );
}

CompileOptions defaultOptions() {
    CompileOptions options;
    options.requestedFidelity = FidelityLevel::f2Stochastic;
    options.strictUnits = true;
    options.strictSafety = true;
    options.deterministicPack = true;
    options.permitHypotheticalParameters = true;
    options.requireTermination = true;
    return options;
}

bool decodeOptions(const NVivoCompileOptions* source,
                   CompileOptions& destination,
                   Diagnostics& diagnostics) {
    destination = defaultOptions();
    if (source == nullptr) return true;

    if (source->struct_size < sizeof(NVivoCompileOptions)) {
        diagnostics.error(
            "NVA001",
            "NVivoCompileOptions.struct_size is smaller than the current ABI structure."
        );
        return false;
    }
    if (source->abi_version != NVIVO_COMPILER_ABI_VERSION) {
        diagnostics.error(
            "NVA002",
            "NVivoCompileOptions.abi_version is unsupported."
        );
        return false;
    }
    if (source->requested_fidelity > NVIVO_FIDELITY_F4_TISSUE) {
        diagnostics.error("NVA003", "Requested fidelity value is invalid.");
        return false;
    }

    destination.requestedFidelity = static_cast<FidelityLevel>(source->requested_fidelity);
    destination.strictUnits = (source->flags & NVIVO_COMPILE_STRICT_UNITS) != 0;
    destination.strictSafety = (source->flags & NVIVO_COMPILE_STRICT_SAFETY) != 0;
    destination.deterministicPack = (source->flags & NVIVO_COMPILE_DETERMINISTIC_PACK) != 0;
    destination.permitHypotheticalParameters =
        (source->flags & NVIVO_COMPILE_PERMIT_HYPOTHETICAL_PARAMETERS) != 0;
    destination.requireTermination = (source->flags & NVIVO_COMPILE_REQUIRE_TERMINATION) != 0;

    const auto& limits = source->limits;
    destination.limits.maximumSpecies = limits.maximum_species;
    destination.limits.maximumParameters = limits.maximum_parameters;
    destination.limits.maximumReactions = limits.maximum_reactions;
    destination.limits.maximumRules = limits.maximum_rules;
    destination.limits.maximumConstraints = limits.maximum_constraints;
    destination.limits.maximumExpressionInstructions = limits.maximum_expression_instructions;
    destination.limits.maximumStoichiometryTerms = limits.maximum_stoichiometry_terms;
    destination.limits.maximumTemporalStates = limits.maximum_temporal_states;

    if (destination.limits.maximumSpecies == 0 ||
        destination.limits.maximumParameters == 0 ||
        destination.limits.maximumReactions == 0 ||
        destination.limits.maximumExpressionInstructions == 0 ||
        destination.limits.maximumStoichiometryTerms == 0 ||
        destination.limits.maximumTemporalStates == 0) {
        diagnostics.error("NVA004", "Compiler resource limits must be non-zero.");
        return false;
    }
    return true;
}

PipelineResult compilePipeline(const std::uint8_t* source,
                               std::size_t sourceSize,
                               const NVivoCompileOptions* abiOptions) {
    PipelineResult pipeline;
    CompileOptions options;
    if (!decodeOptions(abiOptions, options, pipeline.diagnostics)) {
        pipeline.status = NVIVO_STATUS_INVALID_ARGUMENT;
        return pipeline;
    }

    if (source == nullptr && sourceSize != 0) {
        pipeline.diagnostics.error("NVA005", "Source pointer is null while source_size is non-zero.");
        pipeline.status = NVIVO_STATUS_INVALID_ARGUMENT;
        return pipeline;
    }

    const std::string_view sourceText(
        source == nullptr ? "" : reinterpret_cast<const char*>(source),
        sourceSize
    );
    auto parsed = json::parse(sourceText);
    pipeline.diagnostics.append(parsed.diagnostics);
    if (!parsed.root.has_value()) {
        pipeline.status = NVIVO_STATUS_PARSE_ERROR;
        return pipeline;
    }

    UnitRegistry units;
    auto decoded = decodeProgram(*parsed.root, units);
    pipeline.diagnostics.append(decoded.diagnostics);
    if (!decoded.program.has_value()) {
        pipeline.status = NVIVO_STATUS_VALIDATION_ERROR;
        return pipeline;
    }
    pipeline.program = std::move(decoded.program);

    const auto sourceBytes = std::span<const std::byte>(
        reinterpret_cast<const std::byte*>(sourceText.data()),
        sourceText.size()
    );
    Compiler compiler(std::move(units));
    auto compiled = compiler.compile(*pipeline.program, sourceBytes, options);
    pipeline.safety = std::move(compiled.safety);
    pipeline.diagnostics.append(compiled.diagnostics);
    pipeline.ir = std::move(compiled.ir);

    if (!pipeline.ir.has_value()) {
        pipeline.status = pipeline.safety.hasBlockingFinding
            ? NVIVO_STATUS_SAFETY_REJECTED
            : NVIVO_STATUS_COMPILE_ERROR;
        return pipeline;
    }
    pipeline.status = NVIVO_STATUS_OK;
    return pipeline;
}

std::string_view statusName(NVivoStatus status) {
    switch (status) {
        case NVIVO_STATUS_OK: return "ok";
        case NVIVO_STATUS_INVALID_ARGUMENT: return "invalid-argument";
        case NVIVO_STATUS_PARSE_ERROR: return "parse-error";
        case NVIVO_STATUS_VALIDATION_ERROR: return "validation-error";
        case NVIVO_STATUS_SAFETY_REJECTED: return "safety-rejected";
        case NVIVO_STATUS_COMPILE_ERROR: return "compile-error";
        case NVIVO_STATUS_INVALID_PACK: return "invalid-pack";
        case NVIVO_STATUS_RESOURCE_LIMIT: return "resource-limit";
        case NVIVO_STATUS_OUT_OF_MEMORY: return "out-of-memory";
        case NVIVO_STATUS_INTERNAL_ERROR: return "internal-error";
    }
    return "unknown";
}

std::string pipelineReport(const PipelineResult& pipeline,
                           const std::array<std::uint8_t, kFingerprintBytes>* packFingerprint = nullptr,
                           std::size_t packBytes = 0) {
    std::ostringstream stream;
    stream << "{\"status\":\"" << statusName(pipeline.status) << "\""
           << ",\"compiler\":" << pipeline.diagnostics.toJson()
           << ",\"safety\":" << pipeline.safety.toJson();
    if (pipeline.ir.has_value()) {
        stream << ",\"program\":{\"fidelity\":\"" << fidelityName(pipeline.ir->fidelity)
               << "\",\"sourceFingerprint\":\"" << hexFingerprint(pipeline.ir->sourceFingerprint)
               << "\",\"species\":" << pipeline.ir->species.size()
               << ",\"parameters\":" << pipeline.ir->parameters.size()
               << ",\"reactions\":" << pipeline.ir->reactions.size()
               << ",\"rules\":" << pipeline.ir->rules.size()
               << ",\"monitors\":" << pipeline.ir->monitors.size()
               << ",\"cohorts\":" << pipeline.ir->cohorts.size() << '}';
    }
    if (packFingerprint != nullptr) {
        stream << ",\"pack\":{\"bytes\":" << packBytes
               << ",\"fingerprint\":\"" << hexFingerprint(*packFingerprint) << "\"}";
    }
    stream << '}';
    return stream.str();
}

NVivoStatus copyReportOrOOM(std::string_view report,
                            NVivoByteBuffer* output,
                            NVivoStatus intendedStatus) {
    if (!copyBuffer(report, output)) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    }
    return intendedStatus;
}

} // namespace

extern "C" {

void nvivo_default_compile_options(NVivoCompileOptions* options) {
    if (options == nullptr) return;
    std::memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->abi_version = NVIVO_COMPILER_ABI_VERSION;
    options->flags = NVIVO_COMPILE_STRICT_UNITS |
                     NVIVO_COMPILE_STRICT_SAFETY |
                     NVIVO_COMPILE_DETERMINISTIC_PACK |
                     NVIVO_COMPILE_PERMIT_HYPOTHETICAL_PARAMETERS |
                     NVIVO_COMPILE_REQUIRE_TERMINATION;
    options->requested_fidelity = NVIVO_FIDELITY_F2_STOCHASTIC;
    options->limits.maximum_species = 16'384;
    options->limits.maximum_parameters = 65'536;
    options->limits.maximum_reactions = 65'536;
    options->limits.maximum_rules = 16'384;
    options->limits.maximum_constraints = 16'384;
    options->limits.maximum_expression_instructions = 1'048'576;
    options->limits.maximum_stoichiometry_terms = 1'048'576;
    options->limits.maximum_temporal_states = 262'144;
}

NVivoStatus nvivo_compile_program_json(const std::uint8_t* source,
                                       std::size_t sourceSize,
                                       const NVivoCompileOptions* options,
                                       NVivoByteBuffer* programPack,
                                       NVivoByteBuffer* diagnosticsJson) {
    clearBuffer(programPack);
    clearBuffer(diagnosticsJson);
    if (programPack == nullptr) return NVIVO_STATUS_INVALID_ARGUMENT;

    try {
        PipelineResult pipeline = compilePipeline(source, sourceSize, options);
        if (pipeline.status != NVIVO_STATUS_OK || !pipeline.ir.has_value()) {
            return copyReportOrOOM(
                pipelineReport(pipeline),
                diagnosticsJson,
                pipeline.status
            );
        }

        auto packed = serializeProgramPack(*pipeline.ir);
        pipeline.diagnostics.append(packed.diagnostics);
        if (packed.diagnostics.hasErrors() || packed.bytes.empty()) {
            pipeline.status = NVIVO_STATUS_COMPILE_ERROR;
            return copyReportOrOOM(
                pipelineReport(pipeline),
                diagnosticsJson,
                pipeline.status
            );
        }

        const std::string report = pipelineReport(
            pipeline,
            &packed.fingerprint,
            packed.bytes.size()
        );
        if (!copyBuffer(packed.bytes, programPack)) {
            return NVIVO_STATUS_OUT_OF_MEMORY;
        }
        if (!copyBuffer(report, diagnosticsJson)) {
            nvivo_buffer_release(programPack);
            return NVIVO_STATUS_OUT_OF_MEMORY;
        }
        return NVIVO_STATUS_OK;
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& error) {
        const std::string report = std::string("{\"status\":\"internal-error\",\"message\":\"") +
                                   json::escape(error.what()) + "\"}";
        return copyReportOrOOM(report, diagnosticsJson, NVIVO_STATUS_INTERNAL_ERROR);
    } catch (...) {
        return copyReportOrOOM(
            "{\"status\":\"internal-error\",\"message\":\"unknown exception\"}",
            diagnosticsJson,
            NVIVO_STATUS_INTERNAL_ERROR
        );
    }
}

NVivoStatus nvivo_validate_program_json(const std::uint8_t* source,
                                        std::size_t sourceSize,
                                        const NVivoCompileOptions* options,
                                        NVivoByteBuffer* diagnosticsJson) {
    clearBuffer(diagnosticsJson);
    try {
        PipelineResult pipeline = compilePipeline(source, sourceSize, options);
        return copyReportOrOOM(
            pipelineReport(pipeline),
            diagnosticsJson,
            pipeline.status
        );
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (const std::exception& error) {
        const std::string report = std::string("{\"status\":\"internal-error\",\"message\":\"") +
                                   json::escape(error.what()) + "\"}";
        return copyReportOrOOM(report, diagnosticsJson, NVIVO_STATUS_INTERNAL_ERROR);
    } catch (...) {
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}

NVivoStatus nvivo_inspect_program_pack(const std::uint8_t* pack,
                                       std::size_t packSize,
                                       std::uint32_t verifySectionHashes,
                                       NVivoPackSummary* summary,
                                       NVivoByteBuffer* inspectionJson) {
    clearBuffer(inspectionJson);
    if (summary != nullptr) {
        std::memset(summary, 0, sizeof(*summary));
        summary->struct_size = sizeof(*summary);
    }
    if (pack == nullptr && packSize != 0) return NVIVO_STATUS_INVALID_ARGUMENT;

    try {
        const auto bytes = std::span<const std::byte>(
            reinterpret_cast<const std::byte*>(pack),
            packSize
        );
        const auto inspection = inspectProgramPack(bytes, verifySectionHashes != 0);
        if (summary != nullptr && packSize >= sizeof(PackHeader)) {
            summary->pack_major = inspection.header.major;
            summary->pack_minor = inspection.header.minor;
            summary->compiler_abi = inspection.header.compilerABI;
            summary->fidelity = inspection.header.fidelity;
            summary->section_count = inspection.header.sectionCount;
            summary->total_bytes = inspection.header.totalBytes;
            std::copy(
                inspection.header.sourceFingerprint.begin(),
                inspection.header.sourceFingerprint.end(),
                summary->source_fingerprint
            );
            std::copy(
                inspection.header.contentFingerprint.begin(),
                inspection.header.contentFingerprint.end(),
                summary->content_fingerprint
            );
        }
        const std::string report = inspection.toJson();
        const NVivoStatus status = inspection.valid ? NVIVO_STATUS_OK : NVIVO_STATUS_INVALID_PACK;
        return copyReportOrOOM(report, inspectionJson, status);
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (...) {
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}

void nvivo_buffer_release(NVivoByteBuffer* buffer) {
    if (buffer == nullptr) return;
    std::free(buffer->data);
    buffer->data = nullptr;
    buffer->size = 0;
}

} // extern "C"
