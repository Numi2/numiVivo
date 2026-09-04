#pragma once
#include "NumiVivoCore/Core.hpp"
#include <iterator>

namespace nvivo {
// Called only after the directory and hashes have passed structural inspection.
// Uses bounded copies, never unaligned casts into an untrusted payload.
void validateProgramPackTables(std::span<const std::byte> bytes,
                               const PackHeader& header,
                               const std::vector<PackSectionDescriptor>& sections,
                               Diagnostics& diagnostics);
}
