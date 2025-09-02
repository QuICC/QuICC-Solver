/**
 * @file IDenseOpOperator.hpp
 * @brief Implementation of the generic interface dense spectral operator
 */

#ifndef QUICC_DENSEOP_IDENSEOPOPERATOR_HPP
#define QUICC_DENSEOP_IDENSEOPOPERATOR_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {
/// @brief namespace for generic dense spectral operator builders
namespace DenseOp {

/// @brief base class for dense operator
class IDenseOpOperator
{
public:
   /// @brief ctor
   IDenseOpOperator() = default;

   /// @brief dtor
   virtual ~IDenseOpOperator() = default;
};

} // namespace DenseOp
} // namespace QuICC

#endif // define QUICC_DENSEOP_IDENSEOPOPERATOR_HPP
