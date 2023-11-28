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
/// @brief namespace for Cuda backends
namespace Cuda {

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

} // namespace Cuda
} // namespace Pointwise
} // namespace QuICC
