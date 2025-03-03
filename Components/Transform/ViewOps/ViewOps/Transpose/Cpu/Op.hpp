/**
 * @file Op.hpp
 * @brief Transpose operations on Views
 */
#pragma once

// External includes
//

// Project includes
//
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "View/View.hpp"
#include "ViewOps/Transpose/Tags.hpp"
#include "ViewOps/Transpose/Cpu/Impl.hpp"

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
   if constexpr (std::is_same_v<Perm, p201_t> &&
                 std::is_same_v<typename Tin::AttributesType, View::DCCSC3D> &&
                 std::is_same_v<typename Tout::AttributesType, View::DCCSC3D>)
   {
      details::implPerm201(out, in);
   }
   else if constexpr (std::is_same_v<Perm, p120_t> &&
                      std::is_same_v<typename Tin::AttributesType,
                         View::DCCSC3D> &&
                      std::is_same_v<typename Tout::AttributesType,
                         View::DCCSC3D>)
   {
      details::implPerm120(out, in);
   }
   else if constexpr (std::is_same_v<Perm, p201_t> &&
                      std::is_same_v<typename Tin::AttributesType,
                         View::S1CLCSC3D> &&
                      std::is_same_v<typename Tout::AttributesType,
                         View::DCCSC3D>)
   {
      details::implPerm201(out, in);
   }
   else if constexpr (std::is_same_v<Perm, p120_t> &&
                      std::is_same_v<typename Tout::AttributesType,
                         View::S1CLCSC3D> &&
                      std::is_same_v<typename Tin::AttributesType,
                         View::DCCSC3D>)
   {
      details::implPerm120(out, in);
   }
   else
   {
      throw std::logic_error("transpose not implemented");
   }
}

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
