/**
 * @file Pointwise.hpp
 * @brief Pointwise operations on Views
 * Allows for any user defined pointwise operation.
 * The operation is defined via a functor object.
 * Value semantic lets a (good) compiler easily inline and
 * remove the indirect call.
 */
#pragma once

// External includes
//

// Project includes
//
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
/// @brief namespace for Pointwise type operations
namespace Pointwise {
/// @brief namespace for cpu backends
namespace Cpu {

using namespace QuICC::Operator;

/// @brief Pointwise operator
/// @tparam Tout View
/// @tparam Tin View
/// @tparam Functor unary scalar functor
template <class Tout, class Tin, class Functor>
class Op : public UnaryBaseOp<Op<Tout, Tin, Functor>, Tout, Tin>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Tin::ScalarType var)
   Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Tin::ScalarType var)
   Op(Functor f) : _f(f){};
   /// @brief default constructor
   Op() = delete;
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
   friend UnaryBaseOp<Op<Tout, Tin, Functor>, Tout, Tin>;
};

template <class Tout, class Tin, class Functor>
void Op<Tout, Tin, Functor>::applyImpl(Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix("Pointwise::Cpu::applyImpl");

   // check Tout and Tin match Functor op
   using res_t = std::invoke_result_t<Functor, typename Tin::ScalarType>;
   static_assert(std::is_same_v<typename Tout::ScalarType, res_t>,
      "Mismatch in functor or arguments");
   // check same size
   assert(out.size() == in.size());

   // apply pointwise functor
   for (std::size_t i = 0; i < out.size(); ++i)
   {
      out[i] = _f(in[i]);
   }
}

} // namespace Cpu
} // namespace Pointwise
} // namespace QuICC
