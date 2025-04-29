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
#include "Operator/Nary.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
/// @brief namespace for Pointwise type operations
namespace Pointwise {
/// @brief namespace for Cuda backends
namespace Cuda {

using namespace QuICC::Operator;

/// @brief Pointwise operator
/// @tparam Functor Nary scalar functor
/// @tparam Tout output View
/// @tparam ...Targs input Views
template <class Functor, class Tout, class... Targs>
class Op : public NaryBaseOp<Op<Functor, Tout, Targs...>, Tout, Targs...>
{
private:
   /// @brief stored functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Functor _f;

public:
   /// @brief capture functor by value
   /// @param f functor, i.e. struct with method
   /// Tout::ScalarType operator()(Targs::ScalarType var, ...)
   Op(Functor f) : _f(f){};
   /// @brief default constructor
   Op() = delete;
   /// @brief dtor
   ~Op() = default;

private:
   /// @brief action implementation
   /// @param out output View
   /// @param ...args input Views
   void applyImpl(Tout& out, const Targs&... args);
   /// @brief give access to base class
   friend NaryBaseOp<Op<Functor, Tout, Targs...>, Tout, Targs...>;
};

} // namespace Cuda
} // namespace Pointwise
} // namespace QuICC
