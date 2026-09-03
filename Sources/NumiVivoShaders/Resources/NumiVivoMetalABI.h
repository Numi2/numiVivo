#ifndef NUMIVIVO_METAL_ABI_H
#define NUMIVIVO_METAL_ABI_H

#include <metal_stdlib>
using namespace metal;

constant uint NVIVO_METAL_ABI_VERSION = 1u;
constant uint NVIVO_INVALID_INDEX = 0xffffffffu;
constant uint NVIVO_MAX_EXPRESSION_STACK = 64u;
constant uint NVIVO_MAX_POISSON_ATTEMPTS = 16u;

constant uint NVIVO_SPECIES_CONSERVED = 1u << 0u;
constant uint NVIVO_SPECIES_EXTERNALLY_OWNED = 1u << 1u;
constant uint NVIVO_SPECIES_INPUT = 1u << 2u;
constant uint NVIVO_SPECIES_STATE = 1u << 3u;
constant uint NVIVO_SPECIES_OUTPUT = 1u << 4u;
constant uint NVIVO_SPECIES_COUNT_VALUED = 1u << 5u;

constant uint NVIVO_REACTION_CRITICAL = 1u << 0u;
constant uint NVIVO_REACTION_HAS_GATE = 1u << 1u;
constant uint NVIVO_REACTION_DELAYED = 1u << 2u;
constant uint NVIVO_REACTION_STOCHASTIC_ELIGIBLE = 1u << 3u;
constant uint NVIVO_REACTION_SPATIAL = 1u << 4u;

constant uint NVIVO_ACTION_TARGET_IS_STRING = 1u << 0u;
constant uint NVIVO_MONITOR_IS_TERMINATION = 1u << 0u;
constant ushort NVIVO_EXPRESSION_REFERENCE_IS_TIME = 1u << 0u;

constant uint NVIVO_MODE_F1_DETERMINISTIC = 1u;
constant uint NVIVO_MODE_F2_STOCHASTIC = 2u;

constant uint NVIVO_DIAGNOSTIC_NONFINITE = 1u << 0u;
constant uint NVIVO_DIAGNOSTIC_BOUND_VIOLATION = 1u << 1u;
constant uint NVIVO_DIAGNOSTIC_MONITOR = 1u << 2u;
constant uint NVIVO_DIAGNOSTIC_REJECT = 1u << 3u;
constant uint NVIVO_DIAGNOSTIC_SUBSTEP = 1u << 4u;
constant uint NVIVO_DIAGNOSTIC_REVERSIBLE_SHUTDOWN = 1u << 5u;
constant uint NVIVO_DIAGNOSTIC_PERMANENT_SHUTDOWN = 1u << 6u;
constant uint NVIVO_DIAGNOSTIC_EXPRESSION_FAULT = 1u << 7u;
constant uint NVIVO_DIAGNOSTIC_STOCHASTIC_FALLBACK = 1u << 8u;
constant uint NVIVO_DIAGNOSTIC_FLUX_TRUNCATION = 1u << 9u;
constant uint NVIVO_DIAGNOSTIC_EVENT_OVERFLOW = 1u << 10u;

constant uint NVIVO_PUBLICATION_COMMITTED = 0u;
constant uint NVIVO_PUBLICATION_REJECTED = 1u;
constant uint NVIVO_PUBLICATION_SUBSTEP_REQUIRED = 2u;
constant uint NVIVO_PUBLICATION_REVERSIBLE_SHUTDOWN = 3u;
constant uint NVIVO_PUBLICATION_PERMANENT_SHUTDOWN = 4u;

struct NVivoSpeciesRecord {
    uint nameOffset;
    uint compartmentOffset;
    uint unitOffset;
    uint flags;
    float initialValue;
    float minimum;
    float maximum;
    float reserved;
};

struct NVivoGPUParameterRecord {
    float value;
    float minimum;
    float maximum;
    uint flags;
};

struct NVivoStoichiometryRecord {
    uint speciesIndex;
    short coefficient;
    ushort role;
};

struct NVivoReactionRecord {
    uint nameOffset;
    uint compartmentOffset;
    uint reactantOffset;
    uint productOffset;
    uint reactantCount;
    uint productCount;
    uint parameterOffset;
    uint parameterCount;
    uint expressionOffset;
    uint expressionCount;
    uint rateLaw;
    uint flags;
    float delaySeconds;
    float characteristicRate;
    uint cohortIndex;
    uint gateExpressionOffset;
};

struct NVivoExpressionInstruction {
    ushort opcode;
    ushort flags;
    uint operand;
    float immediate;
    uint auxiliary;
};

struct NVivoActionRecord {
    uint targetIndex;
    uint expressionOffset;
    uint expressionCount;
    uint kind;
    float constantValue;
    float maximumRate;
    uint unitOffset;
    uint flags;
};

struct NVivoRuleRecord {
    uint nameOffset;
    uint conditionOffset;
    uint conditionCount;
    uint actionOffset;
    uint actionCount;
    int priority;
    float refractorySeconds;
    uint temporalStateOffset;
};

struct NVivoMonitorRecord {
    uint nameOffset;
    uint expressionOffset;
    uint expressionCount;
    uint messageOffset;
    uint severity;
    uint response;
    uint temporalStateOffset;
    uint flags;
};

struct NVivoCohortRecord {
    uint reactionOffset;
    uint reactionCount;
    uint rateLaw;
    uint flags;
    float maximumStableStep;
    float stiffnessEstimate;
    uint preferredThreads;
    uint reserved;
};

struct NVivoIncidenceRecord {
    uint reactionIndex;
    short netCoefficient;
    ushort reserved16;
    uint reserved0;
    uint reserved1;
};

struct NVivoStepUniforms {
    uint activeCellCount;
    uint cellCapacity;
    uint speciesCount;
    uint parameterCount;

    uint reactionCount;
    uint ruleCount;
    uint monitorCount;
    uint temporalStateCount;

    float deltaTime;
    float absoluteTime;
    uint logicalStepLow;
    uint logicalStepHigh;

    uint seedLow;
    uint seedHigh;
    uint mode;
    uint substepIndex;

    uint eventCapacity;
    uint featureFlags;
    uint delaySlotCount;
    uint delayWriteSlot;
};

struct NVivoCohortUniforms {
    uint reactionOffset;
    uint reactionCount;
    uint dispatchCellCount;
    uint reserved;
};

struct NVivoRuntimeDiagnostics {
    atomic_uint flags;
    atomic_uint nonFiniteCount;
    atomic_uint boundViolationCount;
    atomic_uint monitorViolationCount;
    atomic_uint shutdownCount;
    atomic_uint expressionFaultCount;
    atomic_uint stochasticFallbackCount;
    atomic_uint fluxTruncationCount;
    atomic_uint eventOverflowCount;
    atomic_uint firstCell;
    atomic_uint firstSubject;
    atomic_uint requestedResponse;
    atomic_uint maximumSeverity;
    atomic_uint reserved0;
    atomic_uint reserved1;
    atomic_uint reserved2;
};

struct NVivoPublicationRecord {
    uint committedStepLow;
    uint committedStepHigh;
    uint stateVersion;
    uint status;
    uint diagnosticFlags;
    uint shutdownState;
    uint activeCellCount;
    uint eventCount;
};

struct NVivoEventRecord {
    uint cellIndex;
    uint kind;
    uint subject;
    uint logicalStepLow;
    float value0;
    float value1;
    uint flags;
    uint reserved;
};

struct NVivoProgramArguments {
    device float *committedState [[id(0)]];
    device float *shadowState [[id(1)]];
    device float *committedTemporalState [[id(2)]];
    device float *shadowTemporalState [[id(3)]];
    device const NVivoGPUParameterRecord *parameters [[id(4)]];
    device const NVivoSpeciesRecord *species [[id(5)]];
    device const uint *reactionParameterIndices [[id(6)]];
    device const NVivoStoichiometryRecord *stoichiometry [[id(7)]];
    device const NVivoReactionRecord *reactions [[id(8)]];
    device const NVivoExpressionInstruction *expressions [[id(9)]];
    device const NVivoActionRecord *actions [[id(10)]];
    device const NVivoRuleRecord *rules [[id(11)]];
    device const NVivoMonitorRecord *monitors [[id(12)]];
    device const NVivoCohortRecord *cohorts [[id(13)]];
    device const uint *speciesIncidenceOffsets [[id(14)]];
    device const NVivoIncidenceRecord *speciesIncidence [[id(15)]];
    device float *reactionFlux [[id(16)]];
    device NVivoRuntimeDiagnostics *diagnostics [[id(17)]];
    device NVivoPublicationRecord *publication [[id(18)]];
    device NVivoEventRecord *events [[id(19)]];
    device atomic_uint *eventCount [[id(20)]];
    device const uint *cellActiveMask [[id(21)]];
    device float *committedDelayedFlux [[id(22)]];
    device float *shadowDelayedFlux [[id(23)]];
    device float *committedRuleRefractory [[id(24)]];
    device float *shadowRuleRefractory [[id(25)]];
};

struct NVivoSpatialNeighbor {
    uint nodeIndex;
    uint speciesIndex;
    float conductance;
    float partitionCoefficient;
};

struct NVivoSpatialUniforms {
    uint nodeCount;
    uint nodeCapacity;
    uint speciesCount;
    uint neighborRecordCount;
    float deltaTime;
    float minimumValue;
    uint flags;
    uint reserved;
};

struct NVivoReductionUniforms {
    uint valueCount;
    uint stride;
    uint outputOffset;
    uint flags;
};

#endif
