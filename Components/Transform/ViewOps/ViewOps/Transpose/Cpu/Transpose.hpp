/**
 * @file Transpose.hpp
 * @brief Transpose operations on Views
 */
#pragma once

// External includes
//

// Project includes
//
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "ViewOps/Transpose/Tags.hpp"
#include "View/View.hpp"

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for cpu backends
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Transpose operator
/// @tparam Tout
/// @tparam Tin
template <class Tout, class Tin, class Perm>
class Op : public UnaryBaseOp<Op<Tout, Tin, Perm>, Tout, Tin>
{
public:
   /// @brief default constructor
   Op() = default;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param in input View
   void applyImpl(Tout& out, const Tin& in);
   /// @brief give access to base class
   friend UnaryBaseOp<Op<Tout, Tin, Perm>, Tout, Tin>;
};

/// @brief inclusive scan
/// @tparam T
/// @param K
template <class T>
void pSum(std::vector<T>& vec){
    static_assert(std::is_integral_v<T>, "T must be of integral type");
    assert(vec.size() > 0);
    for (std::size_t i = 1; i < vec.size(); ++i) {
        vec[i] = vec[i-1] + vec[i];
    }
}

template <class Tout, class Tin, class Perm>
void Op<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("Transpose::Cpu::applyImpl");
    if (std::is_same_v<Perm, p201_t> &&
        std::is_same_v<typename Tin::AttributesType, View::DCCSC3D> &&
        std::is_same_v<typename Tout::AttributesType, View::DCCSC3D>) {
        // dense transpose
        assert(out.size() == in.size());
        assert(out.size() == out.dims()[0]*out.dims()[1]*out.dims()[2]);
        auto I = in.dims()[0];
        auto J = in.dims()[1];
        auto K = in.dims()[2];
        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t j = 0; j < J; ++j) {
                for (std::size_t i = 0; i < I; ++i) {
                    std::size_t ijk = i + j*I + k*I*J;
                    std::size_t jki = j + k*J + i*J*K;
                    assert(ijk < in.size());
                    assert(jki < out.size());
                    out[jki] = in[ijk];
                }
            }
        }
    }
    else if (std::is_same_v<Perm, p201_t> &&
        std::is_same_v<typename Tin::AttributesType, View::S1CLCSC3D> &&
        std::is_same_v<typename Tout::AttributesType, View::DCCSC3D>) {
        // dense transpose
        assert(out.size() == in.size());
        auto I = in.dims()[0];
        auto J = in.dims()[1];
        auto K = in.dims()[2];
        // access S1CLCSC3D
        // cumulative column height is (with ijk) I*k - sum(i)_0^k
        // iSum shifted by 1
        std::vector<std::uint32_t> iSum(K, 0);
        for (std::size_t i = 2; i < K; ++i) {
            iSum[i] = iSum[i-1] + 1;
        }
        pSum(iSum);
        // cumulative row width (with jki)
        // kSum shifted by 1
        std::vector<std::uint32_t> kSum(I);
        kSum[I-1] = K-2;
        for (std::size_t i = I-1; i > 0; --i) {
            if (kSum[i] > 0) {
                kSum[i-1] = kSum[i] - 1;
            }
            else {
                kSum[i-1] = 0;
            }
        }
        pSum(kSum);

        for (std::size_t k = 0; k < K; ++k) {
            std::size_t Iloc = I - k;
            for (std::size_t j = 0; j < J; ++j) {
                for (std::size_t i = 0; i < Iloc; ++i) {
                    std::size_t ijk = i + j*Iloc + (k*I-iSum[k])*J;
                    std::size_t jki = j + k*J + (i*K-kSum[i])*J;
                    assert(ijk < in.size());
                    assert(jki < out.size());
                    out[jki] = in[ijk];
                }
            }
        }
    }
    else {
        throw std::logic_error("transpose not implemented");
    }
}

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
