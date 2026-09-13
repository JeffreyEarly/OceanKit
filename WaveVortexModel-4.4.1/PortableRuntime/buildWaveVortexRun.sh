#!/bin/sh
set -eu

script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repository_root=$(CDPATH= cd -- "$script_directory/.." && pwd)
manifest="$repository_root/CompiledKernel/native-fftw-provider.env"
. "$manifest"

build_directory=${1:-"$repository_root/.compiled-backend-cache/runtime-build"}
cache_root=${WV_RUNTIME_CACHE_ROOT:-"$repository_root/.compiled-backend-cache"}
archive="$cache_root/downloads/fftw-$WV_FFTW_VERSION.tar.gz"
source_directory="$cache_root/source/fftw-$WV_FFTW_VERSION"
provider_root="$cache_root/provider/$WV_FFTW_PROVIDER_ID"
provider_build="$cache_root/build/$WV_FFTW_PROVIDER_ID"

mkdir -p "$cache_root/downloads" "$cache_root/source" "$provider_build"
if [ ! -f "$provider_root/lib/libfftw3.3.dylib" ] || [ ! -f "$provider_root/lib/libfftw3_threads.3.dylib" ]; then
    if [ ! -f "$archive" ]; then
        /usr/bin/curl --fail --location --output "$archive" "$WV_FFTW_SOURCE_URL"
    fi
    actual=$(/usr/bin/shasum -a 256 "$archive" | /usr/bin/awk '{print $1}')
    if [ "$actual" != "$WV_FFTW_SOURCE_SHA256" ]; then
        echo "FFTW archive checksum mismatch." >&2
        exit 2
    fi
    if [ ! -x "$source_directory/configure" ]; then
        rm -rf "$source_directory"
        /usr/bin/tar -xzf "$archive" -C "$cache_root/source"
    fi
    sdk=$(/usr/bin/xcrun --sdk macosx --show-sdk-path)
    flags="$WV_FFTW_COMPILER_FLAGS -mmacosx-version-min=$WV_FFTW_DEPLOYMENT_TARGET -isysroot $sdk"
    cd "$provider_build"
    env SDKROOT="$sdk" MACOSX_DEPLOYMENT_TARGET="$WV_FFTW_DEPLOYMENT_TARGET" CC=/usr/bin/clang CFLAGS="$flags" "$source_directory/configure" --prefix="$provider_root" $WV_FFTW_CONFIGURE_FLAGS
    /usr/bin/make -j"$(/usr/sbin/sysctl -n hw.logicalcpu)"
    /usr/bin/make install
fi

cmake -S "$repository_root/PortableRuntime" -B "$build_directory" -DCMAKE_BUILD_TYPE=Release -DWV_RUNTIME_ENABLE_NATIVE_FFTW=ON -DWV_RUNTIME_FFTW_ROOT="$provider_root"
cmake --build "$build_directory" --parallel --target wave-vortex-run
printf '%s\n' "$build_directory/wave-vortex-run"
