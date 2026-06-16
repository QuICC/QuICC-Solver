/**
 * @file I4D4mDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4D4m sparse
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4D4mDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I4D4mDiags::I4D4mDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
   const int q) :
    IDiags(alpha, dBeta, l, q)
{}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
