/**
 * @file OrthoMgs2Functor.hpp
 * @brief Modified Gram-Schmidt orthogonalization with re-orthogonalization
 */

// System includes
//

// Project includes
//
#include "Timestep/Exponential/Functors/OrthoMgs2Functor.hpp"
#include "Types/Typedefs.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

namespace Functors {

OrthoMgs2Functor::OrthoMgs2Functor(const int p)
   : mcP(p)
{}

double OrthoMgs2Functor::operator()(Matrix& matV, Matrix& matH, const int j, const int n)
{
   // MGS Orthogonalization with re-orthogonalization
   int i0 = std::max(0, j + 1 - this->mcP);
   for(int i = i0; i <= j; i++)
   {
      Matrix colH = details::computeAugmentedDot(matV, i, i, matV, j+1, n);

      matH(i, j) = colH(0, 0);

      matV.col(j+1) -= matH(i,j)*matV.col(i);
   }
   for(int i = i0; i <= j; i++)
   {
      Matrix colH = details::computeAugmentedDot(matV, i, i, matV, j+1, n);

      matH(i, j) += colH(0, 0);

      matV.col(j+1) -= colH(0,0)*matV.col(i);
   }

   // Norm
   double normV = details::computeAugmented2Norm(matV, j+1, n);
   return normV;
}

int OrthoMgs2Functor::p() const
{
   return this->mcP;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
