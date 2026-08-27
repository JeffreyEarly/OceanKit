file(GLOB_RECURSE forcing_sources
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/*.hpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.cpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.hpp")

set(forcing_class_literals
    WVNonlinearAdvection
    WVAdaptiveDamping
    WVFixedAmplitudeForcing
    WVBottomFrictionQuadratic
    WVBottomFrictionLinear
    WVPseudoTopographicWaveGeneration
    WVBetaPlanePVAdvection)
list(APPEND forcing_class_literals WVNarrowBandGeostrophicForcing)

foreach(source IN LISTS forcing_sources)
    if(source MATCHES "/WVForcingContracts\.cpp$")
        continue()
    endif()
    file(READ "${source}" contents)
    foreach(class_name IN LISTS forcing_class_literals)
        string(FIND "${contents}" "\"${class_name}\"" occurrence)
        if(NOT occurrence EQUAL -1)
            message(FATAL_ERROR
                "${source} duplicates paired forcing class literal ${class_name}; "
                "register it in WVForcingContracts instead.")
        endif()
    endforeach()
endforeach()

foreach(source IN LISTS forcing_sources)
    file(READ "${source}" contents)
    foreach(closed_token "WVForcingKind" "WVForcingPayload")
        string(FIND "${contents}" "${closed_token}" occurrence)
        if(NOT occurrence EQUAL -1)
            message(FATAL_ERROR
                "${source} retains closed forcing dispatch token ${closed_token}.")
        endif()
    endforeach()
endforeach()

file(READ "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVForcingEngine.cpp" engine)
foreach(token "typeIdentifier ==" "typeIdentifier !=" "matlabClassName")
    string(FIND "${engine}" "${token}" occurrence)
    if(NOT occurrence EQUAL -1)
        message(FATAL_ERROR
            "WVForcingEngine performs name-based forcing dispatch: ${token}")
    endif()
endforeach()

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVBarotropicQGForcingEngine.cpp"
    qg_engine)
foreach(token "typeIdentifier ==" "typeIdentifier !=" "matlabClassName")
    string(FIND "${qg_engine}" "${token}" occurrence)
    if(NOT occurrence EQUAL -1)
        message(FATAL_ERROR
            "WVBarotropicQGForcingEngine performs name-based forcing dispatch: ${token}")
    endif()
endforeach()
