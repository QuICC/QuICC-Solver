/**
 * @file Reduction.hpp
 * @brief Reduction operations on Views
 * Allows for a reduction operation.
 */
#pragma once

// External includes
//
#include <cstdint>

// Project includes
//
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"
#include "View/Attributes.hpp"

namespace QuICC {
/// @brief namespace for Reduction type operations
namespace Reduction {
/// @brief namespace for cpu backends
namespace Cpu {

/// @brief Reduction operator
/// @tparam Tout n dimensional view
/// @tparam Tin  n-1 dimensional view
/// @tparam Dir axis to perform reduction on
template <class Tout, class Tin, std::uint32_t Dir>
class Op : public Operator::UnaryBaseOp<Op<Tout, Tin, Dir>, Tout, Tin>
{
public:
   /// @brief default constructor
   Op() = default;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief action implementation that might modify the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   // void applyImpl(Tout& out, Tin& in);
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<Op<Tout, Tin, Dir>, Tout, Tin>;
};

template <class Tout, class Tin, std::uint32_t Dir>
void Op<Tout, Tin, Dir>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Reduction::Cpu::applyImpl");

   using namespace QuICC::View;

   // check types consistency
   static_assert(out.rank() == in.rank() - 1, "input/output rank mismatch");
   static_assert(
      std::is_same_v<typename Tin::ScalarType, typename Tout::ScalarType>,
      "input/output scalar type mismatch");

   if constexpr (Dir == 0u &&
                 std::is_same_v<typename Tin::AttributesType, DCCSC3D> &&
                 std::is_same_v<typename Tout::AttributesType, CSC>)
   {
      // check minimal meta data consistency
      assert(out.pointers()[0].size() == in.pointers()[1].size());
      assert(out.indices()[0].size() == in.indices()[1].size());
      assert(out.size() == in.size() / in.lds());

      std::size_t offSet3D = 0;
      const auto M = in.dims()[0];
      for (std::size_t i = 0; i < in.indices()[1].size(); ++i)
      {
         typename Tin::ScalarType acc = 0;
         for (std::size_t i = 0; i < M; ++i)
         {
            acc += in[offSet3D + i];
         }
         out[i] = acc;
         offSet3D += in.lds();
      }
   }
   else
   {
      // reduction not implemented
      static_assert(std::is_same_v<typename Tout::AttributesType, void>,
         "Reduction not implemented for this type");
   }

   // apply Reduction
}

} // namespace Cpu
} // namespace Reduction
} // namespace QuICC
