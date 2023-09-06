/** 
 * @file TMesher.cpp
 * @brief Source of the T spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/1D/TMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"

namespace QuICC {

namespace SpatialScheme {

   TMesher::TMesher(const GridPurpose::Id purpose)
      : I1DMesher(purpose), mNx(-1)
   {
   }

   void TMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& N = this->mDims.at(0);

      // Get standard dealiased FFT size
      this->mNx = Transform::Fft::Tools::dealiasCosFft(N+1);
      // Check for optimised FFT sizes
      this->mNx = Transform::Fft::Tools::optimizeFft(this->mNx);
   }

   int TMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int TMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int TMesher::nDealias1D() const
   {
      return this->mNx;
   }

} // SpatialScheme
} // QuICC
