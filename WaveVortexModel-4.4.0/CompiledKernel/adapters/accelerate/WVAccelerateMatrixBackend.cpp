#include "WVAccelerateMatrixBackend.hpp"
#include <climits>
#include <new>
#if WV_HAVE_ACCELERATE
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
#endif

namespace wavevortex {
#if WV_HAVE_ACCELERATE
namespace {
class AccelerateBackend final : public WVVerticalMatrixBackend {
public:
    const char* identifier() const noexcept override { return "accelerate-gemm-v1"; }
    std::size_t maximumDimension() const noexcept override { return INT_MAX; }
    std::size_t persistentBytes() const noexcept override { return sizeof(*this); }
    bool supportsConcurrentCalls() const noexcept override { return true; }
    void split(std::size_t m, std::size_t k, std::size_t n, const double* a, const double* br, const double* bi,
        std::size_t ldb, double* cr, double* ci, std::size_t ldc, double beta) const noexcept override {
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,static_cast<int>(m),static_cast<int>(n),static_cast<int>(k),1,a,static_cast<int>(m),br,static_cast<int>(ldb),beta,cr,static_cast<int>(ldc));
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,static_cast<int>(m),static_cast<int>(n),static_cast<int>(k),1,a,static_cast<int>(m),bi,static_cast<int>(ldb),beta,ci,static_cast<int>(ldc));
    }
    void interleaved(std::size_t m, std::size_t k, std::size_t n, const WVComplex64* a, const WVComplex64* b,
        std::size_t ldb, WVComplex64* c, std::size_t ldc, double beta) const noexcept override {
        static_assert(sizeof(WVComplex64) == sizeof(__LAPACK_double_complex));
        static_assert(alignof(WVComplex64) >= alignof(__LAPACK_double_complex));
        const __LAPACK_double_complex alpha{1,0}, betaComplex{beta,0};
        cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,static_cast<int>(m),static_cast<int>(n),static_cast<int>(k),&alpha,reinterpret_cast<const __LAPACK_double_complex*>(a),static_cast<int>(m),reinterpret_cast<const __LAPACK_double_complex*>(b),static_cast<int>(ldb),&betaComplex,reinterpret_cast<__LAPACK_double_complex*>(c),static_cast<int>(ldc));
    }
};
}
#endif
WVKernelStatus WVCreateAccelerateMatrixBackend(std::unique_ptr<WVVerticalMatrixBackend>& result) {
#if WV_HAVE_ACCELERATE
    try { result = std::make_unique<AccelerateBackend>(); return WVKernelStatus::ok(); }
    catch (const std::bad_alloc&) { return {WVKernelStatusCode::allocationFailure,"Accelerate backend allocation failed."}; }
#else
    (void)result;
    return {WVKernelStatusCode::unsupportedOperation,"Accelerate is unavailable in this build."};
#endif
}
} // namespace wavevortex
