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



    }
    else {
        // throw std::logic_error("transpose not implemented");
        std::cerr << "transpose op not implemented\n";
        assert(false);
    }
}

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
