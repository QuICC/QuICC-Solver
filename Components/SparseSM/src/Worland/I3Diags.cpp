/**
 * @file I3Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I3 sparse
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I3Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I3Diags::I3Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
   const int q) :
    IDiags(alpha, dBeta, l, q)
{}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
