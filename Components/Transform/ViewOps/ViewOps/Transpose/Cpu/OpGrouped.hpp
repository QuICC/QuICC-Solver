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
#include "ViewOps/Transpose/Cpu/Impl.hpp"
#include "ViewOps/Transpose/Tags.hpp"

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
class OpGrouped : public UnaryBaseOp<OpGrouped<Tout, Tin, Perm>, Tout, Tin>
{
public:
   /// @brief default constructor
   OpGrouped() = default;
   /// @brief dtor
   ~OpGrouped() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param in input View
   void applyImpl(Tout& out, const Tin& in);
   /// @brief give access to base class
   friend UnaryBaseOp<OpGrouped<Tout, Tin, Perm>, Tout, Tin>;
};

template <class Tout, class Tin, class Perm>
void OpGrouped<Tout, Tin, Perm>::applyImpl(Tout& out, const Tin& in)
{
   assert(out.size() == in.size());
   Profiler::RegionFixture<4> fix("Transpose::Cpu::applyImpl");
   if constexpr (std::is_same_v<Perm, p201_t>)
   {
      for (std::size_t i = 0; i < in.size(); ++i)
      {
         details::implPerm201(out[i], in[i]);
      }
   }
   else if constexpr (std::is_same_v<Perm, p120_t>)
   {
      for (std::size_t i = 0; i < in.size(); ++i)
      {
         details::implPerm120(out[i], in[i]);
      }
   }
   else
   {
      throw std::logic_error("transpose not implemented");
   }
}

} // namespace Cpu
} // namespace Transpose
} // namespace QuICC
