#pragma once
#include "NumiVivoCore/Core.hpp"
namespace nvivo {
// Mandatory semantic preflight, independent of optional safety-report severity.
// It runs before any source double is narrowed into a ProgramPack FP32 record.
void validateExecutableSource(const Program&, Diagnostics&);
}
