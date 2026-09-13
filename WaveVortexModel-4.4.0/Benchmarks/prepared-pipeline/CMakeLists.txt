cmake_minimum_required(VERSION 3.20)
project(WVPreparedPipeline LANGUAGES CXX)

set(WVM_SOURCE "" CACHE PATH "WVM authoring checkout")
if(NOT EXISTS "${WVM_SOURCE}/PortableRuntime/CMakeLists.txt")
  message(FATAL_ERROR "WVM_SOURCE must identify an authoring checkout")
endif()

option(BUILD_TESTING "Build the WaveVortex test targets" OFF)
set(WV_RUNTIME_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(WV_RUNTIME_ENABLE_NATIVE_FFTW ON CACHE BOOL "" FORCE)
add_subdirectory("${WVM_SOURCE}/PortableRuntime" runtime)

execute_process(
  COMMAND git -C "${WVM_SOURCE}" rev-parse HEAD
  OUTPUT_VARIABLE WVM_SOURCE_REVISION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE WVM_GIT_RESULT)
if(NOT WVM_GIT_RESULT EQUAL 0)
  set(WVM_SOURCE_REVISION "unknown")
endif()
execute_process(
  COMMAND git -C "${WVM_SOURCE}" status --porcelain
  OUTPUT_VARIABLE WVM_SOURCE_STATUS
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(WVM_SOURCE_STATUS STREQUAL "")
  set(WVM_SOURCE_DIRTY 0)
else()
  set(WVM_SOURCE_DIRTY 1)
endif()
execute_process(
  COMMAND git -C "${WVM_SOURCE}" diff --no-ext-diff HEAD
  OUTPUT_VARIABLE WVM_SOURCE_DIFF)
string(SHA256 WVM_SOURCE_DIFF_SHA256 "${WVM_SOURCE_DIFF}")
string(SHA256 WVM_SOURCE_STATUS_SHA256 "${WVM_SOURCE_STATUS}")

add_executable(wv-prepared-pipeline WVPreparedPipelineWorker.cpp)
file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/WVPreparedPipelineWorker.cpp"
  WV_BENCHMARK_WORKER_SHA256)
file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt"
  WV_BENCHMARK_CMAKE_SHA256)
target_compile_features(wv-prepared-pipeline PRIVATE cxx_std_17)
target_include_directories(wv-prepared-pipeline PRIVATE
  "${WVM_SOURCE}/PortableRuntime/third_party")
target_link_libraries(wv-prepared-pipeline PRIVATE
  WaveVortex::PortableRuntime WaveVortexNativeFFTW WaveVortex::AccelerateMatrix)
target_compile_definitions(wv-prepared-pipeline PRIVATE
  WV_BENCHMARK_SOURCE_REVISION="${WVM_SOURCE_REVISION}"
  WV_BENCHMARK_SOURCE_DIRTY=${WVM_SOURCE_DIRTY}
  WV_BENCHMARK_SOURCE_DIFF_SHA256="${WVM_SOURCE_DIFF_SHA256}"
  WV_BENCHMARK_SOURCE_STATUS_SHA256="${WVM_SOURCE_STATUS_SHA256}"
  WV_BENCHMARK_WORKER_SHA256="${WV_BENCHMARK_WORKER_SHA256}"
  WV_BENCHMARK_CMAKE_SHA256="${WV_BENCHMARK_CMAKE_SHA256}"
  WV_BENCHMARK_COMPILER="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")
set_target_properties(wv-prepared-pipeline PROPERTIES CXX_EXTENSIONS OFF
  BUILD_RPATH "${WV_RUNTIME_FFTW_ROOT}/lib")
if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
  target_compile_options(wv-prepared-pipeline PRIVATE
    -Wall -Wextra -Wpedantic -Werror)
endif()
