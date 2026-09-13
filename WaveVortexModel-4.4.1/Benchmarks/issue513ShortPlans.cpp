// Standalone diagnostic for the WVM constant-stratification vertical plans.
//
// This intentionally mirrors verticalSpecification(c,1,1,1,sine): Nz=129,
// scalar transform stride two, and guru64 batch dimensions
// {2,1,1}, {1,258,258}, {1,258,258}.  It is a diagnostic only; it does not
// enter the benchmark build or the production source tree.

#include <fftw3.h>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using Clock = std::chrono::steady_clock;
constexpr std::size_t nz = 129;
constexpr std::size_t fullRows = 256 * 256;
constexpr unsigned flags = FFTW_MEASURE | FFTW_UNALIGNED;

struct Plan {
    fftw_plan value = nullptr;
    // The production planner creates the plan from an aligned origin and
    // executes it with the sine data pointer offset by two scalar doubles.
    std::size_t offset = 0;
    std::size_t rows = 0;
    std::size_t length = 0;
    bool sine = false;

    Plan(std::size_t rowCount, bool sineValue, std::size_t workers)
        : rows(rowCount), length(sineValue ? nz - 2 : nz), sine(sineValue) {
        if (fftw_init_threads() == 0)
            throw std::runtime_error("fftw_init_threads failed");
        fftw_plan_with_nthreads(static_cast<int>(workers));

        // Two interleaved complex channels are represented as two real
        // batches.  For DST-I the first physical endpoint is skipped.
        const std::size_t doubles = 2 * nz * rows;
        storage = static_cast<double*>(fftw_malloc(doubles * sizeof(double)));
        if (storage == nullptr) throw std::bad_alloc();
        std::fill(storage, storage + doubles, 0.0);

        fftw_iodim64 transform{static_cast<ptrdiff_t>(length), 2, 2};
        fftw_iodim64 batches[] = {
            {2, 1, 1},
            {static_cast<ptrdiff_t>(rows),
             static_cast<ptrdiff_t>(2 * nz),
             static_cast<ptrdiff_t>(2 * nz)},
            {1, static_cast<ptrdiff_t>(2 * nz),
             static_cast<ptrdiff_t>(2 * nz)}};
        offset = sine ? 2 : 0;
        const fftw_r2r_kind kind = sine ? FFTW_RODFT00 : FFTW_REDFT00;
        value = fftw_plan_guru64_r2r(
            1, &transform, 3, batches, storage, storage,
            &kind, flags);
        if (value == nullptr) throw std::runtime_error("plan creation failed");
    }

    ~Plan() {
        if (value != nullptr) fftw_destroy_plan(value);
        fftw_free(storage);
    }
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;

    void execute() const { fftw_execute_r2r(value, storage + offset, storage + offset); }

    void reset() {
        const std::size_t doubles = 2 * nz * rows;
        for (std::size_t i = 0; i < doubles; ++i)
            storage[i] = 1.0e-12 * static_cast<double>((i % 17) + 1);
        if (sine) {
            for (std::size_t row = 0; row < rows; ++row) {
                storage[2 * nz * row] = 0.0;
                storage[2 * nz * row + 1] = 0.0;
                storage[2 * nz * row + 2 * (nz - 1)] = 0.0;
                storage[2 * nz * row + 2 * (nz - 1) + 1] = 0.0;
            }
        }
    }

    double* storage = nullptr;
};

struct RunResult {
    double seconds = 0.0;
    bool finite = true;
    bool parityBoundary = true;
};

RunResult roundtripSeconds(Plan& plan, std::size_t repetitions) {
    // REDFT00 and RODFT00 are self-inverse with factor 2*(Nz-1).  Scaling
    // between the two calls prevents growth while keeping the plan hot.
    const double inverseScale = 1.0 / (2.0 * static_cast<double>(nz - 1));
    plan.reset();
    const auto start = Clock::now();
    for (std::size_t i = 0; i < repetitions; ++i) {
        plan.execute();
        plan.execute();
        // A complete type-I round trip is scaled once, after both unnormalized
        // transforms, so repeated iterations retain the same amplitude.
        const std::size_t doubles = 2 * nz * plan.rows;
        for (std::size_t j = 0; j < doubles; ++j) plan.storage[j] *= inverseScale;
    }
    RunResult result;
    result.seconds = std::chrono::duration<double>(Clock::now() - start).count();
    const std::size_t doubles = 2 * nz * plan.rows;
    for (std::size_t j = 0; j < doubles; ++j)
        result.finite = result.finite && std::isfinite(plan.storage[j]);
    if (plan.sine) {
        for (std::size_t row = 0; row < plan.rows; ++row) {
            const std::size_t base = 2 * nz * row;
            result.parityBoundary = result.parityBoundary &&
                plan.storage[base] == 0.0 && plan.storage[base + 1] == 0.0 &&
                plan.storage[base + 2 * (nz - 1)] == 0.0 &&
                plan.storage[base + 2 * (nz - 1) + 1] == 0.0;
        }
    }
    return result;
}

std::size_t parsePositive(const char* text, const char* name) {
    char* end = nullptr;
    const auto value = std::strtoull(text, &end, 10);
    if (end == text || *end != '\0' || value == 0 ||
        value > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument(std::string("invalid ") + name);
    return static_cast<std::size_t>(value);
}

void printPlan(const char* name, const Plan& plan) {
    char* text = fftw_sprint_plan(plan.value);
    std::cout << "plan_begin=" << name << "\n";
    if (text != nullptr) {
        std::cout << text;
        fftw_free(text);
    }
    std::cout << "plan_end=" << name << "\n";
}

void run(std::size_t rows, std::size_t workers, std::size_t baseRepetitions) {
    Plan dct(rows, false, workers);
    Plan dst(rows, true, workers);
    const std::size_t repetitions =
        baseRepetitions * (fullRows / rows);
    std::cout << std::setprecision(12)
              << "workers=" << workers << " rows=" << rows
              << " nz=" << nz << " full_rows=" << fullRows
              << " repetitions=" << repetitions << " flags=MEASURE|UNALIGNED\n";
    printPlan(rows == fullRows ? "dct_full" : "dct_tile", dct);
    printPlan(rows == fullRows ? "dst_full" : "dst_tile", dst);
    const RunResult dctResult = roundtripSeconds(dct, repetitions);
    const RunResult dstResult = roundtripSeconds(dst, repetitions);
    std::cout << "dct_roundtrip_seconds=" << dctResult.seconds << "\n"
              << "dst_roundtrip_seconds=" << dstResult.seconds << "\n"
              << "dct_roundtrip_seconds_per_call=" << dctResult.seconds / repetitions << "\n"
              << "dst_roundtrip_seconds_per_call=" << dstResult.seconds / repetitions << "\n"
              << "dct_finite=" << (dctResult.finite ? 1 : 0) << "\n"
              << "dst_finite=" << (dstResult.finite ? 1 : 0) << "\n"
              << "dst_parity_boundary_ok=" << (dstResult.parityBoundary ? 1 : 0) << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: issue513ShortPlans THREADS [REPETITIONS]\n";
        return 2;
    }
    try {
        const std::size_t workers = parsePositive(argv[1], "thread count");
        const std::size_t repetitions = argc == 3 ? parsePositive(argv[2], "repetitions") : 1;
        run(1, workers, repetitions);
        run(64, workers, repetitions);
        run(256, workers, repetitions);
        fftw_cleanup_threads();
    } catch (const std::exception& error) {
        std::cerr << "status=failed message=" << error.what() << "\n";
        return 1;
    }
    return 0;
}
