/** 
 * @file AFTMesher.cpp
 * @brief Source of the AFT spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/AFTMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   AFTMesher::AFTMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNr(-1), mNt(-1), mNz(-1)
   {
   }

   void AFTMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& N = this->mDims.at(0);
      int& M = this->mDims.at(1);
      int& Z = this->mDims.at(2);

      // Get standard dealiased FFT size
      this->mNr = Transform::Fft::Tools::dealiasCosFft(N+1);
      // Check for optimised FFT sizes
      this->mNr = Transform::Fft::Tools::optimizeFft(this->mNr);

      // Get mixed dealiased FFT size
      this->mNt = Transform::Fft::Tools::dealiasMixedFft(M+1);
      // Check for optimised FFT sizes
      this->mNt = Transform::Fft::Tools::optimizeFft(this->mNt);

      // Get standard dealiased FFT size
      this->mNz = Transform::Fft::Tools::dealiasCosFft(Z+1);
      // Check for optimised FFT sizes
      this->mNz = Transform::Fft::Tools::optimizeFft(this->mNz);
   }

   int AFTMesher::nPhys1D() const
   {
      return this->mNz;
   }

   int AFTMesher::nPhys2D() const
   {
      return this->mNt;
   }

   int AFTMesher::nPhys3D() const
   {
      return this->mNz;
   }

   int AFTMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int AFTMesher::nSpec2D() const
   {
      const int& M = this->mDims.at(1);
      return M + 1;
   }

   int AFTMesher::nSpec3D() const
   {
      const int& Z = this->mDims.at(2);
      return Z + 1;
   }

   int AFTMesher::nDealias1D() const
   {
      return this->mNr;
   }

   int AFTMesher::nDealias2D() const
   {
      return this->mNt/2 + 1;
   }

   int AFTMesher::nDealias3D() const
   {
      return this->mNz;
   }

} // SpatialScheme
} // QuICC
