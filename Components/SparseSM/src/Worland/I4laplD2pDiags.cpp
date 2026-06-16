/**
 * @file I4laplD2pDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4laplD2p
 * sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4laplD2pDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I4laplD2pDiags::I4laplD2pDiags(const Scalar_t alpha, const Scalar_t dBeta,
   const int l, const int q) :
    IDiags(alpha, dBeta, l, q)
{}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
