cmake_minimum_required(VERSION 3.20)
project(WaveVortexKernel VERSION 0.1.0 LANGUAGES CXX)
find_package(Threads REQUIRED)

add_library(WaveVortexKernel STATIC
    src/WVKernelTypes.cpp
    src/WVTransformBarotropicQGKernel.cpp
    src/WVTransformConstantStratificationKernel.cpp
)
add_library(WaveVortex::Kernel ALIAS WaveVortexKernel)

add_library(WaveVortexReferenceFFT STATIC
    adapters/reference/WVReferenceFFTEngine.cpp
)
add_library(WaveVortex::ReferenceFFT ALIAS WaveVortexReferenceFFT)
target_compile_features(WaveVortexReferenceFFT PUBLIC cxx_std_17)
target_include_directories(WaveVortexReferenceFFT PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters/reference>
)
target_link_libraries(WaveVortexReferenceFFT PUBLIC WaveVortexKernel)
set_target_properties(WaveVortexReferenceFFT PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)

target_compile_features(WaveVortexKernel PUBLIC cxx_std_17)
target_include_directories(WaveVortexKernel PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(WaveVortexKernel PUBLIC Threads::Threads)
set_target_properties(WaveVortexKernel PROPERTIES
    CXX_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON
)

if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
    target_compile_options(WaveVortexKernel PRIVATE -Wall -Wextra -Wpedantic -Werror)
    target_compile_options(WaveVortexReferenceFFT PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
