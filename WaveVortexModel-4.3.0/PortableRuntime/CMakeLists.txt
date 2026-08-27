cmake_minimum_required(VERSION 3.20)
project(WaveVortexPortableRuntime VERSION 1.0.0 LANGUAGES CXX)

include(CTest)
option(WV_RUNTIME_ENABLE_NATIVE_FFTW "Build wave-vortex-run with the pinned native FFTW adapter" OFF)
option(WV_RUNTIME_BUILD_BENCHMARKS "Build author-only portable-runtime benchmark workers" OFF)
set(WV_RUNTIME_FFTW_ROOT "" CACHE PATH "Prefix containing the validated FFTW 3.3.11 include and lib folders")

if(NOT TARGET WaveVortex::Kernel)
    add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../CompiledKernel ${CMAKE_CURRENT_BINARY_DIR}/kernel)
endif()

find_package(netCDF CONFIG QUIET)
if(TARGET netCDF::netcdf)
    set(WAVEVORTEX_NETCDF_TARGET netCDF::netcdf)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)
    set(WAVEVORTEX_NETCDF_TARGET PkgConfig::NETCDF)
endif()

add_library(WaveVortexPortableRuntime STATIC
    include/WaveVortexRuntime/WVExtensionCatalog.hpp
    src/WVExtensionCatalog.cpp
    src/WVPortableImplementationContract.cpp
    src/WVPortableTypedRecord.cpp
    src/WVObservation.cpp
    src/WVOutputSchedule.cpp
    src/WVForcingContracts.cpp
    src/WVNetCDF.cpp
    src/WVForcingScheduleDecoder.cpp
    src/WVCheckpointTransformAdapters.cpp
    src/WVCheckpointWriter.cpp
    src/WVForcingEngine.cpp
    src/WVIntegrationState.cpp
    src/WVRungeKutta.cpp
    include/WaveVortexRuntime/WVBarotropicQGForcingEngine.hpp
    src/WVBarotropicQGForcingEngine.cpp
    src/WVBarotropicQGIntegrationSystem.cpp
    src/WVBarotropicQGFieldEvaluationAdapter.cpp
    src/WVOutputTransformAdapters.cpp
    src/WVOutputOrchestration.cpp
    src/WVObserverAdapter.cpp
    src/WVLegacyObserverCompatibility.hpp
    src/WVObserverContracts.cpp
    src/WVObserverSampleSource.cpp
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
    wv_add_runtime_test(TestWVObserverAdapter TestWVObserverAdapter.cpp observer-adapter)
    wv_add_runtime_test(TestWVExtensionCatalog TestWVExtensionCatalog.cpp extension-catalog)
    wv_add_runtime_test(TestWVCheckpointReader TestWVCheckpointReader.cpp checkpoint-reader)
    wv_add_runtime_test(TestWVCheckpointWriter TestWVCheckpointWriter.cpp checkpoint-writer)
    wv_add_runtime_test(TestWVForcingAndRK4 TestWVForcingAndRK4.cpp forcing-and-rk4)
    wv_add_runtime_test(TestWVUnifiedIntegration TestWVUnifiedIntegration.cpp unified-integration)
    wv_add_runtime_test(TestWVBarotropicQGIntegration TestWVBarotropicQGIntegration.cpp barotropic-qg-integration)
    wv_add_runtime_test(TestWVBarotropicQGForcing TestWVBarotropicQGForcing.cpp barotropic-qg-forcing)
    wv_add_runtime_test(TestWVBarotropicQGModelOutput TestWVBarotropicQGModelOutput.cpp barotropic-qg-model-output)
    wv_add_runtime_test(TestWVOutputOrchestration TestWVOutputOrchestration.cpp output-orchestration)
    wv_add_runtime_test(TestWVOutputOccurrence TestWVOutputOccurrence.cpp output-occurrence)
    wv_add_runtime_test(TestWVModelOutputNetCDF TestWVModelOutputNetCDF.cpp model-output-netcdf)
    wv_add_runtime_test(TestWVModelOutputConfiguration TestWVModelOutputConfiguration.cpp model-output-configuration)
    wv_add_runtime_test(TestWVModel TestWVModel.cpp model-facade)
    wv_add_runtime_test(TestWVFieldEvaluation TestWVFieldEvaluation.cpp field-evaluation)
    wv_add_runtime_test(TestWVObserverOutputEvaluation TestWVObserverOutputEvaluation.cpp observer-output-evaluation)
    wv_add_runtime_test(TestWVObservationOccurrence TestWVObservationOccurrence.cpp observation-occurrence)
    wv_add_runtime_test(TestWVLagrangianParticles TestWVLagrangianParticles.cpp lagrangian-particles)
    wv_add_runtime_test(TestWaveVortexRun TestWaveVortexRun.cpp standalone-runner)

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
    target_include_directories(TestWVRunRequest PRIVATE app third_party)
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
