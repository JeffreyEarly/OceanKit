file(GLOB_RECURSE observer_sources
    "${WV_REPOSITORY_ROOT}/PortableRuntime/include/*.hpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.cpp"
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/*.hpp")

set(observer_class_literals
    WVCoefficients
    WVEulerianFields
    WVMooring
    WVLagrangianParticles
    WVTracer)

foreach(source IN LISTS observer_sources)
    if(source MATCHES "/WVObserverAdapter\\.(cpp|hpp)$")
        continue()
    endif()
    file(READ "${source}" contents)
    foreach(class_name IN LISTS observer_class_literals)
        string(FIND "${contents}" "\"${class_name}\"" occurrence)
        if(NOT occurrence EQUAL -1)
            message(FATAL_ERROR
                "${source} duplicates built-in observer class literal ${class_name}; "
                "register it in WVObserverAdapter instead.")
        endif()
    endforeach()
endforeach()

# Registration sites may name the implementation they construct, but no
# production source may branch on those serialized identities.
foreach(source IN LISTS observer_sources)
    file(READ "${source}" contents)
    foreach(class_name IN LISTS observer_class_literals)
        if(contents MATCHES
           "typeIdentifier[ \t\r\n]*[!=]=[ \t\r\n]*\"${class_name}\"|\"${class_name}\"[ \t\r\n]*[!=]=[ \t\r\n]*[^;\n]*typeIdentifier")
            message(FATAL_ERROR
                "${source} branches on built-in observer identity ${class_name}; "
                "resolve the operation through WVExtensionCatalog metadata.")
        endif()
    endforeach()
endforeach()

set(generic_observation_consumers
    "PortableRuntime/src/WVModelOutputNetCDFWriter.cpp"
    "PortableRuntime/src/WVModelOutputNetCDFReader.cpp"
    "PortableRuntime/src/WVObserverOutputEvaluationService.cpp"
    "PortableRuntime/src/WVOutputOrchestration.cpp"
    "PortableRuntime/src/WVConstantStratificationIntegrationSystem.cpp")
set(observer_kind_dispatch
    "observerImplementation("
    "recordsCoefficients("
    "recordsEulerianFields("
    "recordsFixedProfiles("
    "recordsFixedPoints("
    "recordsMovingParticles("
    "recordsTracerState(")
foreach(relative_path IN LISTS generic_observation_consumers)
    file(READ "${WV_REPOSITORY_ROOT}/${relative_path}" contents)
    foreach(token IN LISTS observer_kind_dispatch)
        string(FIND "${contents}" "${token}" occurrence)
        if(NOT occurrence EQUAL -1)
            message(FATAL_ERROR
                "${relative_path} branches on observer implementation kind via ${token}; "
                "route observation data through schemas and batches instead.")
        endif()
    endforeach()
endforeach()

file(READ
    "${WV_REPOSITORY_ROOT}/PortableRuntime/src/WVObserverOutputEvaluationService.cpp"
    evaluator_contents)
foreach(token IN ITEMS WVObserverSamplingTopology "executionPlan(")
    string(FIND "${evaluator_contents}" "${token}" occurrence)
    if(NOT occurrence EQUAL -1)
        message(FATAL_ERROR
            "WVObserverOutputEvaluationService.cpp retains observer-topology dispatch via ${token}; "
            "resolved observer providers must declare plans and produce batches.")
    endif()
endforeach()

set(retired_dispatch_symbols
    WVObserverKind
    WVObserverStateContract
    WVObserverOutputRule
    WVObserverSamplingTopology
    WVObserverIntegratedOperation
    legacyExecutionPlan
    recordsCoefficients
    recordsEulerianFields
    recordsFixedProfiles
    recordsFixedPoints
    recordsMovingParticles
    recordsTracerState
    contributesRightHandSide
    ownsParticleState
    ownsTracerState)

foreach(source IN LISTS observer_sources)
    file(READ "${source}" contents)
    foreach(symbol IN LISTS retired_dispatch_symbols)
        string(FIND "${contents}" "${symbol}" occurrence)
        if(NOT occurrence EQUAL -1)
            message(FATAL_ERROR
                "${source} retains retired closed observer dispatch symbol ${symbol}.")
        endif()
    endforeach()
endforeach()
