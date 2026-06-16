/**
 * @file I4D3pDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4D3p sparse
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4D3pDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I4D3pDiags::I4D3pDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
   const int q) :
    IDiags(alpha, dBeta, l, q)
{}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
