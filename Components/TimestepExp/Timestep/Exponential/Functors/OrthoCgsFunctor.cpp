/**
 * @file OrthoCgsFunctor.cpp
 * @brief Classical Gram-Schmidt orthogonalization
 */

// System includes
//

// Project includes
//
#include "Timestep/Exponential/Functors/OrthoCgsFunctor.hpp"
#include "Types/Typedefs.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

namespace Functors {

OrthoCgsFunctor::OrthoCgsFunctor(const int p)
   : mcP(p)
{}

double OrthoCgsFunctor::operator()(Matrix& matV, Matrix& matH, const int j, const int n)
{
   // CGS Orthogonalization
   int i0 = std::max(0, j + 1 - this->mcP);
   Matrix colH = details::computeAugmentedDot(matV, i0, j, matV, j+1, n);
   for(int i = i0; i <= j; i++)
   {
      matH(i, j) = colH(i-i0, 0);

      matV.col(j+1) -= matH(i,j)*matV.col(i);
   }

   // Norm
   double normV = details::computeAugmented2Norm(matV, j+1, n);
   return normV;
}

int OrthoCgsFunctor::p() const
{
   return this->mcP;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
