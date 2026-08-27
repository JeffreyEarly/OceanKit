cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED WV_REPOSITORY_ROOT)
    message(FATAL_ERROR "WV_REPOSITORY_ROOT is required.")
endif()

set(generic_sources
    "PortableRuntime/src/WVRungeKutta.cpp"
    "PortableRuntime/include/WaveVortexRuntime/WVRungeKutta.hpp"
    "PortableRuntime/include/WaveVortexRuntime/WVIntegrationContracts.hpp")
foreach(relative_path IN LISTS generic_sources)
    file(READ "${WV_REPOSITORY_ROOT}/${relative_path}" source)
    foreach(token
            "WVTransformBarotropicQG"
            "WVBarotropicQGIntegrationSystem"
            "coefficients.Ap"
            "coefficients.Am"
            "coefficients.A0")
        string(FIND "${source}" "${token}" position)
        if(NOT position EQUAL -1)
            message(FATAL_ERROR
                "Generic integrator contract contains transform dispatch or legacy coefficient access '${token}': ${relative_path}")
        endif()
    endforeach()
endforeach()

# Transform-specific coefficient names, legacy shapes, and persisted-transform
# dispatch belong only in explicitly named adapters. These generic model,
# orchestration, reader, writer, and CLI sources consume resolved layouts and
# services.
set(generic_boundary_sources
    "PortableRuntime/src/WVModel.cpp"
    "PortableRuntime/src/WVOutputOrchestration.cpp"
    "PortableRuntime/src/WVModelOutputNetCDFReader.cpp"
    "PortableRuntime/src/WVModelOutputNetCDFWriter.cpp"
    "PortableRuntime/app/WaveVortexRun.cpp")
foreach(relative_path IN LISTS generic_boundary_sources)
    file(READ "${WV_REPOSITORY_ROOT}/${relative_path}" source)
    foreach(token
            "\"Ap\""
            "\"Am\""
            "\"A0\""
            "coefficients.Ap"
            "coefficients.Am"
            "coefficients.A0"
            "hasLegacyCoefficientTriple()"
            "coefficientShape()")
        string(FIND "${source}" "${token}" position)
        if(NOT position EQUAL -1)
            message(FATAL_ERROR
                "Generic runtime boundary contains transform-specific coefficient behavior '${token}': ${relative_path}")
        endif()
    endforeach()
endforeach()

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVBarotropicQGIntegrationSystem.cpp"
    system_source)
foreach(required "coefficientFamilyView" "forcingEngine_->evaluateRightHandSide")
    string(FIND "${system_source}" "${required}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR
            "Barotropic QG system does not remain behind the resolved coefficient-family boundary: ${required}")
    endif()
endforeach()

foreach(forbidden "coefficients.Ap" "coefficients.Am" "WVForcingKind")
    string(FIND "${system_source}" "${forbidden}" position)
    if(NOT position EQUAL -1)
        message(FATAL_ERROR
            "Barotropic QG system leaked a legacy family or closed forcing dispatch: ${forbidden}")
    endif()
endforeach()

message(STATUS "Barotropic QG architecture policy passed.")
