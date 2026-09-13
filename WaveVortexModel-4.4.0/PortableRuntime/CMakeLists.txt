cmake_minimum_required(VERSION 3.20)
project(WaveVortexPortableRuntime VERSION 1.0.0 LANGUAGES CXX)

include(CTest)
option(WV_RUNTIME_ENABLE_NATIVE_FFTW "Build wave-vortex-run with the pinned native FFTW adapter" OFF)
option(WV_RUNTIME_BUILD_BENCHMARKS "Build author-only portable-runtime benchmark workers" OFF)
set(WV_RUNTIME_FFTW_ROOT "" CACHE PATH "Prefix containing the validated FFTW 3.3.11 include and lib folders")

if(NOT TARGET WaveVortex::Kernel)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../CompiledKernel ${CMAKE_CURRENT_BINARY_DIR}/kernel)
endif()

set(WV_RUNTIME_COMPACT_VARIABLE_POLICY_DEFAULT OFF)
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64" AND
   WV_RUNTIME_ENABLE_NATIVE_FFTW AND WV_ENABLE_ACCELERATE)
    set(WV_RUNTIME_COMPACT_VARIABLE_POLICY_DEFAULT ON)
endif()
option(WV_RUNTIME_ENABLE_COMPACT_VARIABLE_POLICY
    "Select the compact variable-kernel policy in the native runner"
    ${WV_RUNTIME_COMPACT_VARIABLE_POLICY_DEFAULT})
unset(WV_RUNTIME_COMPACT_VARIABLE_POLICY_DEFAULT)

find_package(netCDF CONFIG QUIET)
if(TARGET netCDF::netcdf)
    set(WAVEVORTEX_NETCDF_TARGET netCDF::netcdf)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)
    set(WAVEVORTEX_NETCDF_TARGET PkgConfig::NETCDF)
endif()

add_library(WaveVortexPortableRuntime STATIC
    include/WaveVortexRuntime/WVNativeVariablePolicy.hpp
    src/WVNativeVariablePolicy.cpp
    include/WaveVortexRuntime/WVExtensionCatalog.hpp
    src/WVExtensionCatalog.cpp
    src/WVPortableImplementationContract.cpp
    src/WVPortableTypedRecord.cpp
    src/WVObservation.cpp
    src/WVOutputSchedule.cpp
    src/WVForcingContracts.cpp
    src/WVNetCDF.cpp
    src/WVStratifiedModalRecord.cpp
    src/WVForcingScheduleDecoder.cpp
    src/WVCheckpointTransformAdapters.cpp
    src/WVCheckpointWriter.cpp
    src/WVForcingEngine.cpp
    src/WVIntegrationState.cpp
    src/WVNoMotionProfile.cpp
    src/WVNoMotionProfileRecovery.cpp
    src/WVNoMotionDensityMoments.cpp
    src/WVDensityEventEvaluation.cpp
    src/WVRungeKutta.cpp
    include/WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp
    src/WVBarotropicQGForcingEngine.cpp
    src/WVStratifiedQGForcingEngine.cpp
    src/WVHydrostaticForcingEngine.cpp
    src/WVBoussinesqForcingEngine.cpp
    src/WVBarotropicQGIntegrationSystem.cpp
    src/WVStratifiedQGIntegrationSystem.cpp
    src/WVHydrostaticIntegrationSystem.cpp
    src/WVBoussinesqIntegrationSystem.cpp
    src/WVBarotropicQGFieldEvaluationAdapter.cpp
    src/WVStratifiedFieldEvaluationAdapter.cpp
    src/WVOutputTransformAdapters.cpp
    src/WVOutputOrchestration.cpp
    src/WVObserverAdapter.cpp
    src/WVLegacyObserverCompatibility.hpp
    src/WVObserverContracts.cpp
    src/WVObserverSampleSource.cpp
    src/WVDiagnosticFieldPlan.cpp
    src/WVFieldEvaluationService.cpp
    src/WVConstantStratificationIntegrationSystem.cpp
    src/WVObserverOutputEvaluationService.cpp
    src/WVModelOutputTransformAdapters.cpp
    src/WVModelOutputNetCDFReader.cpp
    src/WVLegacyObservationNetCDFAdapter.cpp
    src/WVModelOutputNetCDFWriter.cpp
    src/WVModelOutputConfiguration.cpp
    src/WVModelInternalAccess.hpp
    src/WVModel.cpp
    src/WVModelTransformAdapters.cpp
)
add_library(WaveVortex::PortableRuntime ALIAS WaveVortexPortableRuntime)

target_compile_features(WaveVortexPortableRuntime PUBLIC cxx_std_17)
target_include_directories(WaveVortexPortableRuntime PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(WaveVortexPortableRuntime PUBLIC WaveVortex::Kernel PRIVATE ${WAVEVORTEX_NETCDF_TARGET})
set_target_properties(WaveVortexPortableRuntime PROPERTIES
    CXX_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON
)

add_library(WaveVortexRunner STATIC
    include/WaveVortexRuntime/WVRunner.hpp
    app/WaveVortexRun.cpp
    app/WVRunRequest.cpp
    app/WVRunnerVariablePolicy.cpp
)
add_library(WaveVortex::Runner ALIAS WaveVortexRunner)
target_compile_features(WaveVortexRunner PUBLIC cxx_std_17)
target_include_directories(WaveVortexRunner PRIVATE
    app
    src
    third_party
)
target_link_libraries(WaveVortexRunner PUBLIC WaveVortex::PortableRuntime PRIVATE WaveVortex::ReferenceFFT)
set_target_properties(WaveVortexRunner PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)

add_executable(wave-vortex-run app/WaveVortexRunMain.cpp)
target_compile_features(wave-vortex-run PRIVATE cxx_std_17)
target_link_libraries(wave-vortex-run PRIVATE WaveVortex::Runner)
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/.. RESULT_VARIABLE WV_RUNTIME_GIT_RESULT OUTPUT_VARIABLE WV_RUNTIME_SOURCE_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT WV_RUNTIME_GIT_RESULT EQUAL 0 OR NOT WV_RUNTIME_SOURCE_COMMIT)
    file(READ "${CMAKE_CURRENT_LIST_DIR}/../resources/mpackage.json" WV_RUNTIME_PACKAGE_MANIFEST)
    string(JSON WV_RUNTIME_PACKAGE_VERSION GET "${WV_RUNTIME_PACKAGE_MANIFEST}" version)
    set(WV_RUNTIME_SOURCE_COMMIT "WaveVortexModel-${WV_RUNTIME_PACKAGE_VERSION}")
endif()
target_compile_definitions(WaveVortexRunner PRIVATE WV_RUNTIME_SOURCE_COMMIT="${WV_RUNTIME_SOURCE_COMMIT}")
target_compile_definitions(WaveVortexRunner PRIVATE WV_RUNTIME_HAS_DENSE_OUTPUT=1)
target_compile_definitions(WaveVortexRunner PRIVATE
    WV_RUNTIME_ENABLE_COMPACT_VARIABLE_POLICY=$<BOOL:${WV_RUNTIME_ENABLE_COMPACT_VARIABLE_POLICY}>)

if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
    if(NOT APPLE OR NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        message(FATAL_ERROR "The optimized native FFTW runtime provider currently supports Apple silicon only.")
    endif()
    find_path(WV_RUNTIME_FFTW_INCLUDE fftw3.h PATHS ${WV_RUNTIME_FFTW_ROOT}/include NO_DEFAULT_PATH REQUIRED)
    find_library(WV_RUNTIME_FFTW_BASE NAMES fftw3 libfftw3.3 PATHS ${WV_RUNTIME_FFTW_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    find_library(WV_RUNTIME_FFTW_THREADS NAMES fftw3_threads libfftw3_threads.3 PATHS ${WV_RUNTIME_FFTW_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    add_library(WaveVortexNativeFFTW STATIC ../CompiledKernel/adapters/native-fftw/WVNativeFFTWEngine.cpp)
    target_compile_features(WaveVortexNativeFFTW PUBLIC cxx_std_17)
    target_include_directories(WaveVortexNativeFFTW PUBLIC ../CompiledKernel/adapters/native-fftw ${WV_RUNTIME_FFTW_INCLUDE})
    target_link_libraries(WaveVortexNativeFFTW PUBLIC WaveVortex::Kernel ${WV_RUNTIME_FFTW_THREADS} ${WV_RUNTIME_FFTW_BASE} ${CMAKE_DL_LIBS})
    set_target_properties(WaveVortexNativeFFTW PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(WaveVortexRunner PRIVATE WaveVortexNativeFFTW)
    target_compile_definitions(WaveVortexRunner PRIVATE WV_RUNTIME_HAS_NATIVE_FFTW=1 WV_RUNTIME_EXPECTED_FFTW_ROOT="${WV_RUNTIME_FFTW_ROOT}")
    set_target_properties(wave-vortex-run PROPERTIES BUILD_RPATH "${WV_RUNTIME_FFTW_ROOT}/lib")
else()
    target_compile_definitions(WaveVortexRunner PRIVATE WV_RUNTIME_HAS_NATIVE_FFTW=0 WV_RUNTIME_EXPECTED_FFTW_ROOT="")
endif()

if(WV_RUNTIME_ENABLE_COMPACT_VARIABLE_POLICY)
    if(NOT WV_RUNTIME_ENABLE_NATIVE_FFTW)
        message(FATAL_ERROR "The compact variable runner policy requires the native FFTW provider.")
    endif()
    if(NOT APPLE OR NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        message(FATAL_ERROR "The compact variable runner policy currently supports Apple silicon only.")
    endif()
    if(NOT WV_ENABLE_ACCELERATE)
        message(FATAL_ERROR "The compact variable runner policy requires the Accelerate matrix adapter.")
    endif()
    target_include_directories(WaveVortexRunner PRIVATE
        ../CompiledKernel/adapters/accelerate)
    target_link_libraries(WaveVortexRunner PRIVATE WaveVortex::AccelerateMatrix)
endif()

if(WV_RUNTIME_BUILD_BENCHMARKS)
    if(NOT WV_RUNTIME_ENABLE_NATIVE_FFTW)
        message(FATAL_ERROR "Portable-runtime benchmarks require the validated native FFTW provider.")
    endif()
    add_executable(wv-standalone-nonlinear-flux-benchmark
        ../Benchmarks/standalone/WVStandaloneNonlinearFluxBenchmark.cpp)
    target_compile_features(wv-standalone-nonlinear-flux-benchmark PRIVATE cxx_std_17)
    target_include_directories(wv-standalone-nonlinear-flux-benchmark PRIVATE
        src
        ../CompiledKernel/adapters/native-fftw
    )
    target_link_libraries(wv-standalone-nonlinear-flux-benchmark PRIVATE
        WaveVortex::PortableRuntime
        WaveVortexNativeFFTW
        ${WAVEVORTEX_NETCDF_TARGET}
    )
    target_compile_definitions(wv-standalone-nonlinear-flux-benchmark PRIVATE
        WV_RUNTIME_EXPECTED_FFTW_ROOT="${WV_RUNTIME_FFTW_ROOT}"
        WV_RUNTIME_SOURCE_COMMIT="${WV_RUNTIME_SOURCE_COMMIT}"
    )
    set_target_properties(wv-standalone-nonlinear-flux-benchmark PROPERTIES
        BUILD_RPATH "${WV_RUNTIME_FFTW_ROOT}/lib")
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(wv-standalone-nonlinear-flux-benchmark PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()

    add_executable(wv-standalone-barotropic-qg-benchmark
        ../Benchmarks/standalone/WVStandaloneBarotropicQGBenchmark.cpp)
    target_compile_features(wv-standalone-barotropic-qg-benchmark PRIVATE
        cxx_std_17)
    target_include_directories(wv-standalone-barotropic-qg-benchmark PRIVATE
        ../CompiledKernel/adapters/native-fftw
        ../CompiledKernel/adapters/reference)
    target_link_libraries(wv-standalone-barotropic-qg-benchmark PRIVATE
        WaveVortex::Kernel
        WaveVortex::ReferenceFFT
        WaveVortexNativeFFTW)
    set_target_properties(wv-standalone-barotropic-qg-benchmark PROPERTIES
        BUILD_RPATH "${WV_RUNTIME_FFTW_ROOT}/lib")
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(wv-standalone-barotropic-qg-benchmark PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
    target_compile_options(WaveVortexPortableRuntime PRIVATE -Wall -Wextra -Wpedantic -Werror)
    target_compile_options(WaveVortexRunner PRIVATE -Wall -Wextra -Wpedantic -Werror)
    target_compile_options(wave-vortex-run PRIVATE -Wall -Wextra -Wpedantic -Werror)
    if(TARGET WaveVortexNativeFFTW)
        target_compile_options(WaveVortexNativeFFTW PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(WVBarotropicQGToleranceDump
        tests/WVBarotropicQGToleranceDump.cpp)
    target_compile_features(WVBarotropicQGToleranceDump PRIVATE cxx_std_17)
    target_include_directories(WVBarotropicQGToleranceDump PRIVATE
        ../CompiledKernel/adapters/reference)
    target_link_libraries(WVBarotropicQGToleranceDump PRIVATE
        WaveVortex::PortableRuntime
        WaveVortex::ReferenceFFT
        ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVBarotropicQGToleranceDump PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()

    add_executable(WVBarotropicQGForcingDump
        tests/WVBarotropicQGForcingDump.cpp)
    target_compile_features(WVBarotropicQGForcingDump PRIVATE cxx_std_17)
    target_include_directories(WVBarotropicQGForcingDump PRIVATE
        ../CompiledKernel/adapters/reference)
    target_link_libraries(WVBarotropicQGForcingDump PRIVATE
        WaveVortex::PortableRuntime
        WaveVortex::ReferenceFFT)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVBarotropicQGForcingDump PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()

    function(wv_add_runtime_test target source test_name)
        add_executable(${target} tests/${source})
        target_compile_features(${target} PRIVATE cxx_std_17)
        target_include_directories(${target} PRIVATE
            src
            ../CompiledKernel/adapters/reference
        )
        target_link_libraries(${target} PRIVATE
            WaveVortex::PortableRuntime
            WaveVortex::ReferenceFFT
            ${WAVEVORTEX_NETCDF_TARGET}
        )
        target_compile_definitions(${target} PRIVATE
            WV_CHECKPOINT_FIXTURE_DIR="${CMAKE_CURRENT_LIST_DIR}/tests/fixtures"
            WV_RUNTIME_FIXTURE_DIR="${CMAKE_CURRENT_LIST_DIR}/tests/fixtures"
            WV_RUNTIME_RUNNER="$<TARGET_FILE:wave-vortex-run>"
            WV_TEST_NATIVE_FFTW=$<BOOL:${WV_RUNTIME_ENABLE_NATIVE_FFTW}>
        )
        if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
            target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
        endif()
        add_dependencies(${target} wave-vortex-run)
        add_test(NAME ${test_name} COMMAND ${target})
    endfunction()

    wv_add_runtime_test(TestWVPortableImplementationContract TestWVPortableImplementationContract.cpp portable-implementation-contract)
    wv_add_runtime_test(TestWVPortableTypedRecord TestWVPortableTypedRecord.cpp portable-typed-record)
    wv_add_runtime_test(TestWVObservation TestWVObservation.cpp observation-schema-batch)
    wv_add_runtime_test(TestWVPortableVariableCatalog TestWVPortableVariableCatalog.cpp portable-variable-catalog)
    target_sources(TestWVPortableVariableCatalog PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVPortableVariableCatalog PRIVATE ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVNoMotionProfile TestWVNoMotionProfile.cpp no-motion-profile)
    target_sources(TestWVNoMotionProfile PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVNoMotionProfile PRIVATE ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVNoMotionProfileRecovery TestWVNoMotionProfileRecovery.cpp no-motion-profile-recovery)
    target_sources(TestWVNoMotionProfileRecovery PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVNoMotionProfileRecovery PRIVATE ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVDensityEventEvaluation TestWVDensityEventEvaluation.cpp density-event-evaluation)
    target_sources(TestWVDensityEventEvaluation PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVDensityEventEvaluation PRIVATE src ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVDensityFieldEvaluation TestWVDensityFieldEvaluation.cpp density-field-evaluation)
    target_sources(TestWVDensityFieldEvaluation PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVDensityFieldEvaluation PRIVATE src ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVObserverAdapter TestWVObserverAdapter.cpp observer-adapter)
    wv_add_runtime_test(TestWVExtensionCatalog TestWVExtensionCatalog.cpp extension-catalog)
    wv_add_runtime_test(TestWVForcingCompatibility TestWVForcingCompatibility.cpp forcing-compatibility-catalog)
    add_executable(TestWVStratifiedModalRecord tests/TestWVStratifiedModalRecord.cpp
        ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVStratifiedModalRecord PRIVATE ../tools/compiled-kernel/tests)
    target_link_libraries(TestWVStratifiedModalRecord PRIVATE WaveVortex::PortableRuntime ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVStratifiedModalRecord PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME stratified-modal-record COMMAND TestWVStratifiedModalRecord)
    add_executable(TestWVOwnedStratifiedModalSource tests/TestWVOwnedStratifiedModalSource.cpp)
    target_link_libraries(TestWVOwnedStratifiedModalSource PRIVATE WaveVortex::PortableRuntime ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVOwnedStratifiedModalSource PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME owned-stratified-modal-source COMMAND TestWVOwnedStratifiedModalSource)
    add_executable(TestWVMatlabRawPrimitives tests/TestWVMatlabRawPrimitives.cpp)
    target_link_libraries(TestWVMatlabRawPrimitives PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVMatlabRawPrimitives PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME matlab-raw-primitives COMMAND TestWVMatlabRawPrimitives)
    add_executable(TestWVBarotropicMatlabPrimitives tests/TestWVBarotropicMatlabPrimitives.cpp)
    target_link_libraries(TestWVBarotropicMatlabPrimitives PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVBarotropicMatlabPrimitives PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME barotropic-matlab-primitives COMMAND TestWVBarotropicMatlabPrimitives)
    add_executable(TestWVConstantMatlabPrimitives tests/TestWVConstantMatlabPrimitives.cpp)
    target_link_libraries(TestWVConstantMatlabPrimitives PRIVATE WaveVortex::Kernel WaveVortex::ReferenceFFT)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVConstantMatlabPrimitives PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME constant-matlab-primitives COMMAND TestWVConstantMatlabPrimitives)
    if(NOT TARGET TestWVAllocationProbe)
        add_executable(TestWVAllocationProbe
            ../tools/compiled-kernel/tests/TestWVAllocationProbe.cpp
            ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
        target_compile_features(TestWVAllocationProbe PRIVATE cxx_std_17)
        if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
            target_compile_options(TestWVAllocationProbe PRIVATE -Wall -Wextra -Wpedantic -Werror)
        endif()
    endif()
    add_test(NAME allocation-probe COMMAND TestWVAllocationProbe)
    wv_add_runtime_test(TestWVCheckpointReader TestWVCheckpointReader.cpp checkpoint-reader)
    wv_add_runtime_test(TestWVCheckpointWriter TestWVCheckpointWriter.cpp checkpoint-writer)
    wv_add_runtime_test(TestWVVariableEvaluation TestWVVariableEvaluation.cpp variable-evaluation)
    wv_add_runtime_test(TestWVPreparedFieldCache TestWVPreparedFieldCache.cpp prepared-field-cache)
    target_sources(TestWVPreparedFieldCache PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    wv_add_runtime_test(TestWVForcingDiagnosticSession TestWVForcingDiagnosticSession.cpp forcing-diagnostic-session)
    wv_add_runtime_test(TestWVForcingAndRK4 TestWVForcingAndRK4.cpp forcing-and-rk4)
    wv_add_runtime_test(TestWVForcingTendencies TestWVForcingTendencies.cpp forcing-tendencies)
    wv_add_runtime_test(TestWVUnifiedIntegration TestWVUnifiedIntegration.cpp unified-integration)
    wv_add_runtime_test(TestWVBarotropicQGIntegration TestWVBarotropicQGIntegration.cpp barotropic-qg-integration)
    wv_add_runtime_test(TestWVBarotropicQGForcing TestWVBarotropicQGForcing.cpp barotropic-qg-forcing)
    wv_add_runtime_test(TestWVBarotropicQGModelOutput TestWVBarotropicQGModelOutput.cpp barotropic-qg-model-output)
    wv_add_runtime_test(TestWVOutputOrchestration TestWVOutputOrchestration.cpp output-orchestration)
    target_sources(TestWVOutputOrchestration PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVOutputOrchestration PRIVATE ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVOutputOccurrence TestWVOutputOccurrence.cpp output-occurrence)
    wv_add_runtime_test(TestWVModelOutputNetCDF TestWVModelOutputNetCDF.cpp model-output-netcdf)
    wv_add_runtime_test(TestWVModelOutputConfiguration TestWVModelOutputConfiguration.cpp model-output-configuration)
    wv_add_runtime_test(TestWVModel TestWVModel.cpp model-facade)
    target_sources(TestWVModel PRIVATE ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(TestWVModel PRIVATE ../tools/compiled-kernel/tests)
    wv_add_runtime_test(TestWVFieldEvaluation TestWVFieldEvaluation.cpp field-evaluation)
    wv_add_runtime_test(TestWVObserverOutputEvaluation TestWVObserverOutputEvaluation.cpp observer-output-evaluation)
    wv_add_runtime_test(TestWVObservationOccurrence TestWVObservationOccurrence.cpp observation-occurrence)
    wv_add_runtime_test(TestWVLagrangianParticles TestWVLagrangianParticles.cpp lagrangian-particles)
    wv_add_runtime_test(TestWVRunnerVariablePolicy TestWVRunnerVariablePolicy.cpp runner-variable-policy)
    target_include_directories(TestWVRunnerVariablePolicy PRIVATE app)
    target_link_libraries(TestWVRunnerVariablePolicy PRIVATE WaveVortex::Runner)
    wv_add_runtime_test(TestWaveVortexRun TestWaveVortexRun.cpp standalone-runner)
    target_link_libraries(TestWaveVortexRun PRIVATE WaveVortex::Runner)
    target_compile_definitions(TestWaveVortexRun PRIVATE
        WV_TEST_COMPACT_VARIABLE_POLICY=$<BOOL:${WV_RUNTIME_ENABLE_COMPACT_VARIABLE_POLICY}>)

    add_executable(WVForwardIntegrationProbe tests/WVForwardIntegrationProbe.cpp)
    target_link_libraries(WVForwardIntegrationProbe PRIVATE WaveVortex::Runner)
    add_test(NAME forward-integration-probe COMMAND ${CMAKE_COMMAND}
        -DPROBE=$<TARGET_FILE:WVForwardIntegrationProbe>
        -DFIXTURE=${CMAKE_CURRENT_LIST_DIR}/tests/fixtures/forcing-mixed-hydrostatic.nc
        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/forward-integration-probe-contract
        -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVForwardIntegrationProbe.cmake)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVForwardIntegrationProbe PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()

    if(TARGET WaveVortexNativeFFTW)
        add_executable(TestWVBarotropicQGNativeFFTW
            tests/TestWVBarotropicQGNativeFFTW.cpp)
        target_compile_features(TestWVBarotropicQGNativeFFTW PRIVATE cxx_std_17)
        target_include_directories(TestWVBarotropicQGNativeFFTW PRIVATE
            ../CompiledKernel/adapters/native-fftw
            ../CompiledKernel/adapters/reference)
        target_link_libraries(TestWVBarotropicQGNativeFFTW PRIVATE
            WaveVortex::Kernel
            WaveVortex::ReferenceFFT
            WaveVortexNativeFFTW)
        target_compile_definitions(TestWVBarotropicQGNativeFFTW PRIVATE
            WV_RUNTIME_EXPECTED_FFTW_ROOT="${WV_RUNTIME_FFTW_ROOT}")
        set_target_properties(TestWVBarotropicQGNativeFFTW PROPERTIES
            BUILD_RPATH "${WV_RUNTIME_FFTW_ROOT}/lib")
        if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
            target_compile_options(TestWVBarotropicQGNativeFFTW PRIVATE
                -Wall -Wextra -Wpedantic -Werror)
        endif()
        add_test(NAME barotropic-qg-native-fftw
            COMMAND TestWVBarotropicQGNativeFFTW)
    endif()

    add_executable(TestWVExtendedRunner tests/TestWVExtendedRunner.cpp)
    target_compile_features(TestWVExtendedRunner PRIVATE cxx_std_17)
    target_include_directories(TestWVExtendedRunner PRIVATE tests)
    target_link_libraries(TestWVExtendedRunner PRIVATE WaveVortex::Runner)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVExtendedRunner PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME extended-runner-composition COMMAND TestWVExtendedRunner)
    target_compile_definitions(TestWVModelOutputNetCDF PRIVATE
        WV_RUNTIME_EXTENDED_RUNNER="$<TARGET_FILE:TestWVExtendedRunner>"
    )
    add_dependencies(TestWVModelOutputNetCDF TestWVExtendedRunner)

    add_executable(TestWVRunRequest
        tests/TestWVRunRequest.cpp
        app/WVRunRequest.cpp
    )
    target_compile_features(TestWVRunRequest PRIVATE cxx_std_17)
    target_include_directories(TestWVRunRequest PRIVATE app third_party include)
    target_compile_definitions(TestWVRunRequest PRIVATE
        WV_RUNTIME_FIXTURE_DIR="${CMAKE_CURRENT_LIST_DIR}/tests/fixtures"
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVRunRequest PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME run-request COMMAND TestWVRunRequest)

    add_test(NAME observer-dispatch-source-policy
        COMMAND ${CMAKE_COMMAND}
            -DWV_REPOSITORY_ROOT=${CMAKE_CURRENT_LIST_DIR}/..
            -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVObserverDispatch.cmake
    )
    add_test(NAME forcing-dispatch-source-policy
        COMMAND ${CMAKE_COMMAND}
            -DWV_REPOSITORY_ROOT=${CMAKE_CURRENT_LIST_DIR}/..
            -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVForcingDispatch.cmake
    )
    add_test(NAME runtime-architecture-source-policy
        COMMAND ${CMAKE_COMMAND}
            -DWV_REPOSITORY_ROOT=${CMAKE_CURRENT_LIST_DIR}/..
            -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVRuntimeArchitecture.cmake
    )
    add_test(NAME barotropic-qg-architecture-source-policy
        COMMAND ${CMAKE_COMMAND}
            -DWV_REPOSITORY_ROOT=${CMAKE_CURRENT_LIST_DIR}/..
            -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVBarotropicQGArchitecture.cmake
    )
    add_test(NAME output-schedule-dispatch-source-policy
        COMMAND ${CMAKE_COMMAND}
            -DWV_REPOSITORY_ROOT=${CMAKE_CURRENT_LIST_DIR}/..
            -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVOutputScheduleDispatch.cmake
    )
    add_test(NAME observation-occurrence-source-policy
        COMMAND ${CMAKE_COMMAND}
            -DWV_REPOSITORY_ROOT=${CMAKE_CURRENT_LIST_DIR}/..
            -P ${CMAKE_CURRENT_LIST_DIR}/tests/CheckWVObservationOccurrence.cmake
    )
endif()

# Scientific modal reader probe for MATLAB-authored parity fixtures.
if(BUILD_TESTING)
    add_executable(WVStratifiedModalDump tests/WVStratifiedModalDump.cpp)
    target_include_directories(WVStratifiedModalDump PRIVATE third_party)
    target_link_libraries(WVStratifiedModalDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVStratifiedModalDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(WVStratifiedQGKernelDump tests/WVStratifiedQGKernelDump.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(WVStratifiedQGKernelDump PRIVATE third_party)
    target_link_libraries(WVStratifiedQGKernelDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT WaveVortex::AccelerateMatrix)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVStratifiedQGKernelDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVStratifiedQGKernelDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVStratifiedQGKernelDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVStratifiedQGKernelDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(WVStratifiedQGFieldDump tests/WVStratifiedQGFieldDump.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(WVStratifiedQGFieldDump PRIVATE third_party)
    target_link_libraries(WVStratifiedQGFieldDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVStratifiedQGFieldDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVStratifiedQGFieldDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVStratifiedQGFieldDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVStratifiedQGFieldDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_executable(WVDiagnosticFieldDump tests/WVDiagnosticFieldDump.cpp)
    target_include_directories(WVDiagnosticFieldDump PRIVATE third_party)
    target_link_libraries(WVDiagnosticFieldDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVDiagnosticFieldDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVDiagnosticFieldDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVDiagnosticFieldDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVDiagnosticFieldDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_executable(WVNoMotionProfileDump tests/WVNoMotionProfileDump.cpp)
    target_include_directories(WVNoMotionProfileDump PRIVATE third_party)
    target_link_libraries(WVNoMotionProfileDump PRIVATE WaveVortex::PortableRuntime)
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVNoMotionProfileDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_executable(WVNoMotionRecoveryDump tests/WVNoMotionRecoveryDump.cpp)
    target_include_directories(WVNoMotionRecoveryDump PRIVATE third_party)
    target_link_libraries(WVNoMotionRecoveryDump PRIVATE WaveVortex::PortableRuntime)
    add_executable(WVDensityEventDump tests/WVDensityEventDump.cpp)
    target_include_directories(WVDensityEventDump PRIVATE src third_party)
    target_link_libraries(WVDensityEventDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVDensityEventDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVDensityEventDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVDensityEventDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(WVDensityEventDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVNoMotionRecoveryDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(TestWVStratifiedQGKernel tests/TestWVStratifiedQGKernel.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_link_libraries(TestWVStratifiedQGKernel PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVStratifiedQGKernel PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME stratified-qg-kernel COMMAND TestWVStratifiedQGKernel)
endif()

if(BUILD_TESTING)
    add_executable(WVStableForcingDump tests/WVStableForcingDump.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(WVStableForcingDump PRIVATE third_party)
    target_link_libraries(WVStableForcingDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVStableForcingDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVStableForcingDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVStableForcingDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVStableForcingDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    foreach(probe WVStratifiedQGLifecycleProbe WVHydrostaticLifecycleProbe WVBoussinesqLifecycleProbe)
        add_executable(${probe} tests/WVStratifiedLifecycleProbe.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
        target_include_directories(${probe} PRIVATE third_party)
        target_link_libraries(${probe} PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT)
        if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
            target_link_libraries(${probe} PRIVATE WaveVortexNativeFFTW)
            target_compile_definitions(${probe} PRIVATE WV_TEST_NATIVE_FFTW=1)
        else()
            target_compile_definitions(${probe} PRIVATE WV_TEST_NATIVE_FFTW=0)
        endif()
        if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
            target_compile_options(${probe} PRIVATE -Wall -Wextra -Wpedantic -Werror)
        endif()
    endforeach()
endif()

if(BUILD_TESTING)
    add_executable(WVHydrostaticKernelDump tests/WVHydrostaticKernelDump.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(WVHydrostaticKernelDump PRIVATE third_party)
    target_link_libraries(WVHydrostaticKernelDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT WaveVortex::AccelerateMatrix)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVHydrostaticKernelDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVHydrostaticKernelDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVHydrostaticKernelDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVHydrostaticKernelDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(TestWVHydrostaticKernel tests/TestWVHydrostaticKernel.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_link_libraries(TestWVHydrostaticKernel PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVHydrostaticKernel PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME hydrostatic-kernel COMMAND TestWVHydrostaticKernel)

    add_executable(TestWVHydrostaticRuntime tests/TestWVHydrostaticRuntime.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_link_libraries(TestWVHydrostaticRuntime PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVHydrostaticRuntime PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME hydrostatic-runtime COMMAND TestWVHydrostaticRuntime)
endif()

if(BUILD_TESTING)
    add_executable(WVBoussinesqKernelDump tests/WVBoussinesqKernelDump.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_include_directories(WVBoussinesqKernelDump PRIVATE third_party)
    target_link_libraries(WVBoussinesqKernelDump PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT WaveVortex::AccelerateMatrix)
    if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
        target_link_libraries(WVBoussinesqKernelDump PRIVATE WaveVortexNativeFFTW)
        target_compile_definitions(WVBoussinesqKernelDump PRIVATE WV_TEST_NATIVE_FFTW=1)
    else()
        target_compile_definitions(WVBoussinesqKernelDump PRIVATE WV_TEST_NATIVE_FFTW=0)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(WVBoussinesqKernelDump PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(TestWVBoussinesqKernel tests/TestWVBoussinesqKernel.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_link_libraries(TestWVBoussinesqKernel PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVBoussinesqKernel PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME boussinesq-kernel COMMAND TestWVBoussinesqKernel)
endif()

if(BUILD_TESTING)
    add_executable(TestWVBoussinesqRuntime tests/TestWVBoussinesqRuntime.cpp ../tools/compiled-kernel/tests/WVAllocationProbe.cpp)
    target_link_libraries(TestWVBoussinesqRuntime PRIVATE WaveVortex::PortableRuntime WaveVortex::ReferenceFFT ${WAVEVORTEX_NETCDF_TARGET})
    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(TestWVBoussinesqRuntime PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()
    add_test(NAME boussinesq-runtime COMMAND TestWVBoussinesqRuntime)
    foreach(test_target TestWVHydrostaticKernel TestWVBoussinesqKernel TestWVHydrostaticRuntime TestWVBoussinesqRuntime)
        if(WV_RUNTIME_ENABLE_NATIVE_FFTW)
            target_link_libraries(${test_target} PRIVATE WaveVortexNativeFFTW)
            target_compile_definitions(${test_target} PRIVATE WV_TEST_NATIVE_FFTW=1)
        else()
            target_compile_definitions(${test_target} PRIVATE WV_TEST_NATIVE_FFTW=0)
        endif()
    endforeach()

endif()
