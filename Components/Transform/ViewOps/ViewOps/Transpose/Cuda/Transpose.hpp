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

namespace QuICC {
/// @brief namespace for Transpose type operations
namespace Transpose {
/// @brief namespace for Cuda backends
namespace Cuda {

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

} // namespace Cuda
} // namespace Transpose
} // namespace QuICC
