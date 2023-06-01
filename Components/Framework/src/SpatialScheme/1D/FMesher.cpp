/** 
 * @file FMesher.cpp
 * @brief Source of the F spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/1D/FMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"

namespace QuICC {

namespace SpatialScheme {

   FMesher::FMesher(const GridPurpose::Id purpose)
      : I1DMesher(purpose), mNx(-1)
   {
   }

   void FMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& N = this->mDims.at(0);

      // Get standard dealiased FFT size
      this->mNx = Transform::Fft::Tools::dealiasMixedFft(N+1);
      // Check for optimised FFT sizes
      this->mNx = Transform::Fft::Tools::optimizeFft(this->mNx);
   }

   int FMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int FMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int FMesher::nDealias1D() const
   {
      return this->mNx/2 + 1;
   }

} // SpatialScheme
} // QuICC
