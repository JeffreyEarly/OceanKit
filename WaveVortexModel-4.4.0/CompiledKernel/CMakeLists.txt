cmake_minimum_required(VERSION 3.20)
project(WaveVortexKernel VERSION 0.1.0 LANGUAGES CXX)
find_package(Threads REQUIRED)

# The qualified retained compact schedule is the default. Explicit OFF keeps
# the frozen control available for comparison and reproducibility.
option(WV_KERNEL_COMPACT_CONSTANT_CANDIDATE "Select the retained compact constant schedule; OFF selects the frozen control" ON)
set(WV_KERNEL_COMPACT_HORIZONTAL_WORKERS "12" CACHE STRING "Prepared outer workers for the retained compact horizontal schedule")
set(WV_KERNEL_COMPACT_POINTWISE_WORKERS "12" CACHE STRING "Prepared pointwise workers for the retained compact constant schedule")
foreach(worker_setting WV_KERNEL_COMPACT_HORIZONTAL_WORKERS WV_KERNEL_COMPACT_POINTWISE_WORKERS)
    if(NOT "${${worker_setting}}" MATCHES "^[1-9][0-9]*$")
        message(FATAL_ERROR "${worker_setting} must be a positive integer")
    endif()
endforeach()

add_library(WaveVortexKernel STATIC
    src/WVRetainedHorizontalOperator.cpp
    src/WVPreparedVerticalOperator.cpp
    src/WVOwnedStratifiedModalSource.cpp
    src/WVTransformStratifiedQGKernel.cpp
    src/WVTransformHydrostaticKernel.cpp
    src/WVTransformBoussinesqKernel.cpp
    src/WVKernelTypes.cpp
    src/WVTransformBarotropicQGKernel.cpp
    src/WVTransformConstantStratificationKernel.cpp
)
add_library(WaveVortex::Kernel ALIAS WaveVortexKernel)
target_compile_definitions(WaveVortexKernel PUBLIC
    WV_KERNEL_COMPACT_CONSTANT_CANDIDATE=$<BOOL:${WV_KERNEL_COMPACT_CONSTANT_CANDIDATE}>
    WV_KERNEL_COMPACT_HORIZONTAL_WORKERS=${WV_KERNEL_COMPACT_HORIZONTAL_WORKERS}
    WV_KERNEL_COMPACT_POINTWISE_WORKERS=${WV_KERNEL_COMPACT_POINTWISE_WORKERS})

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

# Optional native matrix service; the portable core never includes vendor APIs.
option(WV_ENABLE_ACCELERATE "Enable the Apple Accelerate matrix adapter" ${APPLE})
if(WV_ENABLE_ACCELERATE AND NOT APPLE)
    message(FATAL_ERROR "Accelerate is available only on Apple platforms.")
endif()
add_library(WaveVortexAccelerateMatrix STATIC adapters/accelerate/WVAccelerateMatrixBackend.cpp)
add_library(WaveVortex::AccelerateMatrix ALIAS WaveVortexAccelerateMatrix)
target_compile_features(WaveVortexAccelerateMatrix PUBLIC cxx_std_17)
target_include_directories(WaveVortexAccelerateMatrix PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/adapters/accelerate>)
target_link_libraries(WaveVortexAccelerateMatrix PUBLIC WaveVortex::Kernel)
target_compile_definitions(WaveVortexAccelerateMatrix PRIVATE WV_HAVE_ACCELERATE=$<BOOL:${WV_ENABLE_ACCELERATE}>)
if(WV_ENABLE_ACCELERATE)
    target_link_libraries(WaveVortexAccelerateMatrix PRIVATE "-framework Accelerate")
endif()
set_target_properties(WaveVortexAccelerateMatrix PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
    target_compile_options(WaveVortexAccelerateMatrix PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
