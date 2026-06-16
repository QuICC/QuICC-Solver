/**
 * @file I4laplDmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4laplDm
 * sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4laplDmDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I4laplDmDiags::I4laplDmDiags(const Scalar_t alpha, const Scalar_t dBeta,
   const int l, const int q) :
    IDiags(alpha, dBeta, l, q)
{}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
