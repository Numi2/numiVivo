#include "NumiVivoCore/NumiVivoPackAccess.h"
#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <span>
#include <vector>

namespace {

float narrowFinite(double value) {
    const double limit = static_cast<double>(std::numeric_limits<float>::max());
    if (!std::isfinite(value)) return 0.0F;
    return static_cast<float>(std::clamp(value, -limit, limit));
}

std::span<const std::byte> bytesOf(const std::uint8_t* data, std::size_t size) {
    return {
        reinterpret_cast<const std::byte*>(data),
        size
    };
}

const nvivo::PackSectionDescriptor* findSection(
    const nvivo::PackInspection& inspection,
    std::uint32_t type) {
    const auto iterator = std::find_if(
        inspection.sections.begin(),
        inspection.sections.end(),
        [type](const auto& section) { return section.type == type; }
    );
    return iterator == inspection.sections.end() ? nullptr : &*iterator;
}

} // namespace

extern "C" {

NVivoStatus nvivo_program_pack_section(const std::uint8_t* pack,
                                       std::size_t packSize,
                                       std::uint32_t sectionType,
                                       NVivoPackSectionView* view) {
    if (view == nullptr) return NVIVO_STATUS_INVALID_ARGUMENT;
    std::memset(view, 0, sizeof(*view));
    if (pack == nullptr && packSize != 0) return NVIVO_STATUS_INVALID_ARGUMENT;

    try {
        const auto inspection = nvivo::inspectProgramPack(bytesOf(pack, packSize), true);
        if (!inspection.valid) return NVIVO_STATUS_INVALID_PACK;

        const auto* section = findSection(inspection, sectionType);
        if (section == nullptr) return NVIVO_STATUS_INVALID_ARGUMENT;
        if (section->offset > packSize || section->size > packSize - section->offset) {
            return NVIVO_STATUS_INVALID_PACK;
        }

        view->data = pack + static_cast<std::size_t>(section->offset);
        view->size = static_cast<std::size_t>(section->size);
        view->stride = section->stride;
        view->count = section->count;
        view->flags = section->flags;
        view->alignment = section->alignment;
        return NVIVO_STATUS_OK;
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (...) {
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}

NVivoStatus nvivo_materialize_gpu_parameters(const std::uint8_t* pack,
                                             std::size_t packSize,
                                             NVivoByteBuffer* output) {
    if (output == nullptr) return NVIVO_STATUS_INVALID_ARGUMENT;
    output->data = nullptr;
    output->size = 0;
    if (pack == nullptr && packSize != 0) return NVIVO_STATUS_INVALID_ARGUMENT;

    try {
        const auto inspection = nvivo::inspectProgramPack(bytesOf(pack, packSize), true);
        if (!inspection.valid) return NVIVO_STATUS_INVALID_PACK;
        const auto* section = findSection(
            inspection,
            static_cast<std::uint32_t>(nvivo::PackSectionType::parameters)
        );
        if (section == nullptr || section->stride != sizeof(nvivo::ParameterRecord)) {
            return NVIVO_STATUS_INVALID_PACK;
        }
        if (section->offset > packSize || section->size > packSize - section->offset) {
            return NVIVO_STATUS_INVALID_PACK;
        }

        const std::size_t count = section->count;
        if (count > std::numeric_limits<std::size_t>::max() / sizeof(NVivoGPUParameterRecord)) {
            return NVIVO_STATUS_RESOURCE_LIMIT;
        }
        const std::size_t outputBytes = count * sizeof(NVivoGPUParameterRecord);
        if (outputBytes == 0) return NVIVO_STATUS_OK;

        auto* destination = static_cast<NVivoGPUParameterRecord*>(std::malloc(outputBytes));
        if (destination == nullptr) return NVIVO_STATUS_OUT_OF_MEMORY;
        const auto* source = pack + static_cast<std::size_t>(section->offset);

        for (std::size_t index = 0; index < count; ++index) {
            nvivo::ParameterRecord record;
            std::memcpy(
                &record,
                source + index * sizeof(nvivo::ParameterRecord),
                sizeof(record)
            );
            destination[index] = {
                narrowFinite(record.value),
                narrowFinite(record.minimum),
                narrowFinite(record.maximum),
                record.flags
            };
        }

        output->data = reinterpret_cast<std::uint8_t*>(destination);
        output->size = outputBytes;
        return NVIVO_STATUS_OK;
    } catch (const std::bad_alloc&) {
        return NVIVO_STATUS_OUT_OF_MEMORY;
    } catch (...) {
        return NVIVO_STATUS_INTERNAL_ERROR;
    }
}

} // extern "C"
