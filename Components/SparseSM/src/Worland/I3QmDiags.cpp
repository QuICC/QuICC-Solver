/**
 * @file I3QmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3Qm sparse
 * operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I3QmDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

I3QmDiags::I3QmDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l,
   const int q) :
    IDiags(alpha, dBeta, l, q)
{}

} // namespace Worland
} // namespace SparseSM
} // namespace QuICC
