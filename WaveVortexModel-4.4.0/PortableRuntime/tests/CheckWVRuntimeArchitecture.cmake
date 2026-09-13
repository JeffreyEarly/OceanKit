cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED WV_REPOSITORY_ROOT)
    message(FATAL_ERROR "WV_REPOSITORY_ROOT is required.")
endif()

file(GLOB_RECURSE kernel_sources LIST_DIRECTORIES false
    "${WV_REPOSITORY_ROOT}/CompiledKernel/include/WaveVortexKernel/*.hpp"
    "${WV_REPOSITORY_ROOT}/CompiledKernel/src/*.cpp")
set(kernel_forbidden "mex.h" "matrix.h" "matlabdataarray" "netcdf.h"
    "fftw3.h" "accelerate/accelerate.h" "dispatch/dispatch.h")
foreach(source_path IN LISTS kernel_sources)
    file(READ "${source_path}" source)
    string(TOLOWER "${source}" source_lower)
    foreach(token IN LISTS kernel_forbidden)
        string(FIND "${source_lower}" "${token}" position)
        if(NOT position EQUAL -1)
            message(FATAL_ERROR "CompiledKernel source has forbidden dependency '${token}': ${source_path}")
        endif()
    endforeach()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/CMakeLists.txt" runtime_cmake)
string(FIND "${runtime_cmake}" "WaveVortex::Checkpoint" checkpoint_alias)
if(NOT checkpoint_alias EQUAL -1)
    message(FATAL_ERROR "The obsolete WaveVortex::Checkpoint alias was restored.")
endif()

file(GLOB_RECURSE runtime_production_sources LIST_DIRECTORIES false
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/*.hpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.cpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.hpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/app/*.cpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/app/*.hpp")
foreach(source_path IN LISTS runtime_production_sources)
    file(READ "${source_path}" source)
    foreach(token
            WVObserverFactoryRegistry
            WVOutputScheduleFactoryRegistry
            WVForcingFactoryRegistry)
        string(FIND "${source}" "${token}" position)
        if(NOT position EQUAL -1)
            message(FATAL_ERROR
                "Portable runtime restored retired process-global registry ${token}: ${source_path}")
        endif()
    endforeach()
    if(source MATCHES "seal[A-Za-z_]*(Factory)?Registration|sealRegistration")
        message(FATAL_ERROR
            "Portable runtime restored a registration-seal API: ${source_path}")
    endif()
    if(source MATCHES
       "static[ \t\r\n]+std::(map|unordered_map|mutex)[ \t\r\n]*[<A-Za-z_]")
        message(FATAL_ERROR
            "Portable runtime contains process-global mutable registry storage: ${source_path}")
    endif()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/app/WaveVortexRun.cpp" runner)
foreach(token "netcdf.h" "nc_open(" "nc_create(" "nc_def_")
    string(FIND "${runner}" "${token}" position)
    if(NOT position EQUAL -1)
        message(FATAL_ERROR "The standalone CLI bypasses the persistence API with '${token}'.")
    endif()
endforeach()
string(FIND "${runner}" "WVConstantStratificationIntegrationSystem::create" direct_system_create)
if(NOT direct_system_create EQUAL -1)
    message(FATAL_ERROR "The standalone CLI bypasses the WVModel façade.")
endif()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVModel.cpp" model)
foreach(token "addNewEvenlySpacedOutputGroup" "WVModelOutputGroup::evenlySpaced")
    string(FIND "${model}" "${token}" position)
    if(NOT position EQUAL -1)
        message(FATAL_ERROR
            "WVModel translates restored output records through a schedule-specific builder: ${token}")
    endif()
endforeach()
string(FIND "${model}" "WVModelOutputConfiguration::compile" record_compiler)
if(record_compiler EQUAL -1)
    message(FATAL_ERROR
        "WVModel does not route restored records through the authoritative output compiler.")
endif()
string(FIND "${model}" "WVOutputDriver driver(" local_output_driver)
if(NOT local_output_driver EQUAL -1)
    message(FATAL_ERROR
        "WVModel destroys retry state by constructing a call-local output driver.")
endif()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVModel.hpp" model_header)
foreach(token "internalIntegrationSystem" "internalIntegrator")
    string(FIND "${model_header}" "${token}" public_internal_service)
    if(NOT public_internal_service EQUAL -1)
        message(FATAL_ERROR
            "WVModel exposes obsolete provisional internal service '${token}'.")
    endif()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVModelOutputNetCDF.hpp" output_source_header)
foreach(token
        "WVOutputValueType"
        "WVObserverOutputCadence"
        "WVObserverOutputAttribute"
        "WVObserverOutputVariableSpecification"
        "WVObserverOutputValueView"
        "specifications(")
    string(FIND "${output_source_header}" "${token}" obsolete_source_adapter)
    if(NOT obsolete_source_adapter EQUAL -1)
        message(FATAL_ERROR
            "WVObserverSampleSource exposes obsolete fixed-shape adapter '${token}'.")
    endif()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/WVExtensionCatalog.hpp" extension_catalog_header)
foreach(token
        "WVLegacyObserverOperationBinder"
        "WVLegacyObserverOperationResolver"
        "WVLegacyObserverPersistenceMetadata"
        "legacyOperationResolver"
        "legacyPersistence")
    string(FIND "${extension_catalog_header}" "${token}" public_legacy_seam)
    if(NOT public_legacy_seam EQUAL -1)
        message(FATAL_ERROR
            "The extension catalog exposes obsolete provisional legacy seam '${token}'.")
    endif()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVModelOutputConfiguration.cpp" output_compiler)
string(FIND "${output_compiler}" "resolveOutputPlan" schema_preflight)
string(FIND "${output_compiler}" "WVPortableObserverDescriptor::create" observer_construction)
if(schema_preflight EQUAL -1 OR observer_construction EQUAL -1 OR
   schema_preflight GREATER observer_construction)
    message(FATAL_ERROR
        "Canonical schema preflight must precede observer-provider construction.")
endif()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVModelOutputNetCDFReader.cpp" output_reader)
string(FIND "${output_reader}" "maximumValidationChunkElements = 4096" bounded_payload_validation)
string(FIND "${output_reader}" "slab(slabSize)" state_sized_payload_validation)
string(FIND "${output_reader}" "resolveOutputPlan" raw_schema_resolution)
if(bounded_payload_validation EQUAL -1 OR
   NOT state_sized_payload_validation EQUAL -1 OR
   NOT raw_schema_resolution EQUAL -1)
    message(FATAL_ERROR
        "Raw output inspection must retain bounded data-only payload validation.")
endif()

foreach(relative_path
        "PortableRuntime/src/WVModelOutputNetCDFWriter.cpp"
        "PortableRuntime/src/WVObserverOutputEvaluationService.cpp")
    file(READ "${WV_REPOSITORY_ROOT}/${relative_path}" source)
    if(source MATCHES "WVPortableObserverRecord[ \t\r\n]+descriptor(Record)?[ \t\r\n]*;")
        message(FATAL_ERROR
            "${relative_path} retains a duplicate compiled observer/output graph.")
    endif()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/CompiledKernel/adapters/native-fftw/wv_compiled_backend_mex.cpp" mex_gateway)
string(FIND "${mex_gateway}" "std::unique_ptr<WVModel>" model_owner)
string(FIND "${mex_gateway}" "std::unique_ptr<WVTransformConstantStratificationKernel>" kernel_owner)
if(model_owner EQUAL -1 OR NOT kernel_owner EQUAL -1)
    message(FATAL_ERROR "The production MEX handle must own WVModel rather than a raw kernel.")
endif()

set(numerical_sources
    "PortableRuntime/src/WVRungeKutta.cpp"
    "PortableRuntime/src/WVIntegrationState.cpp"
    "PortableRuntime/src/WVForcingEngine.cpp"
    "PortableRuntime/src/WVFieldEvaluationService.cpp")
foreach(relative_path IN LISTS numerical_sources)
    file(READ "${WV_REPOSITORY_ROOT}/${relative_path}" source)
    foreach(token "netcdf.h" "WVCheckpointWriter" "WVModelOutputNetCDF")
        string(FIND "${source}" "${token}" position)
        if(NOT position EQUAL -1)
            message(FATAL_ERROR "Portable numerical source owns persistence behavior: ${relative_path}")
        endif()
    endforeach()
endforeach()

set(catalog_header
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/WaveVortexRuntime/generated/WVPortableVariableCatalog.hpp")
if(NOT EXISTS "${catalog_header}")
    message(FATAL_ERROR "The generated portable variable catalog is missing.")
endif()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVFieldEvaluationService.cpp" field_service)
foreach(token "constexpr const char *fieldNames" "request.fieldName ==" "request.fieldName !=")
    string(FIND "${field_service}" "${token}" position)
    if(NOT position EQUAL -1)
        message(FATAL_ERROR "Field evaluation restored runtime name dispatch or duplicate metadata: ${token}")
    endif()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVObserverOutputEvaluationService.cpp" output_service)
foreach(token "struct FieldMetadata" "linearInitialOnly(const std::string")
    string(FIND "${output_service}" "${token}" position)
    if(NOT position EQUAL -1)
        message(FATAL_ERROR "Observer output restored duplicate variable metadata: ${token}")
    endif()
endforeach()

message(STATUS "Portable runtime architecture policy passed.")
