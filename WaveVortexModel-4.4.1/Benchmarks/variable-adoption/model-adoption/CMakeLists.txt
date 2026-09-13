cmake_minimum_required(VERSION 3.20)
project(WVVariableModelAdoption LANGUAGES CXX)

set(WVM_SOURCE "" CACHE PATH "Frozen control or candidate WVM authoring checkout")
if(NOT EXISTS "${WVM_SOURCE}/PortableRuntime/CMakeLists.txt")
    message(FATAL_ERROR "WVM_SOURCE must identify a WVM authoring checkout")
endif()

set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(WV_RUNTIME_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(WV_RUNTIME_ENABLE_NATIVE_FFTW ON CACHE BOOL "" FORCE)
add_subdirectory("${WVM_SOURCE}/PortableRuntime" runtime)

add_executable(wv-variable-model-adoption WVVariableModelAdoptionWorker.cpp)
target_compile_features(wv-variable-model-adoption PRIVATE cxx_std_17)
target_include_directories(wv-variable-model-adoption PRIVATE
    "${WVM_SOURCE}/CompiledKernel/adapters/native-fftw")
target_link_libraries(wv-variable-model-adoption PRIVATE
    WaveVortex::PortableRuntime WaveVortexNativeFFTW)

if(EXISTS "${WVM_SOURCE}/PortableRuntime/include/WaveVortexRuntime/WVVariableKernelServices.hpp")
    target_include_directories(wv-variable-model-adoption PRIVATE
        "${WVM_SOURCE}/CompiledKernel/adapters/accelerate")
    target_link_libraries(wv-variable-model-adoption PRIVATE
        WaveVortex::AccelerateMatrix)
    target_compile_definitions(wv-variable-model-adoption PRIVATE
        WV_MODEL_VARIABLE_SERVICES=1)
else()
    target_compile_definitions(wv-variable-model-adoption PRIVATE
        WV_MODEL_VARIABLE_SERVICES=0)
endif()

execute_process(
    COMMAND git -C "${WVM_SOURCE}" rev-parse HEAD
    OUTPUT_VARIABLE model_source_commit
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY)
execute_process(
    COMMAND git -C "${WVM_SOURCE}" diff --quiet HEAD -- CompiledKernel PortableRuntime
    RESULT_VARIABLE model_source_worktree_dirty)
execute_process(
    COMMAND git -C "${WVM_SOURCE}" diff --cached --quiet HEAD -- CompiledKernel PortableRuntime
    RESULT_VARIABLE model_source_index_dirty)
if(model_source_worktree_dirty EQUAL 0 AND model_source_index_dirty EQUAL 0)
    set(model_compiled_source_dirty 0)
else()
    set(model_compiled_source_dirty 1)
endif()
file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/WVVariableModelAdoptionWorker.cpp"
    model_worker_source_sha256)
target_compile_definitions(wv-variable-model-adoption PRIVATE
    WV_RUNTIME_EXPECTED_FFTW_ROOT="${WV_RUNTIME_FFTW_ROOT}"
    WV_MODEL_SOURCE_COMMIT="${model_source_commit}"
    WV_MODEL_SOURCE_ROOT="${WVM_SOURCE}"
    WV_MODEL_COMPILED_SOURCE_DIRTY=${model_compiled_source_dirty}
    WV_MODEL_WORKER_SOURCE_SHA256="${model_worker_source_sha256}"
    WV_MODEL_BUILD_TYPE="${CMAKE_BUILD_TYPE}")
set_target_properties(wv-variable-model-adoption PROPERTIES
    CXX_EXTENSIONS OFF
    BUILD_RPATH "${WV_RUNTIME_FFTW_ROOT}/lib")

if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
    target_compile_options(wv-variable-model-adoption PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
