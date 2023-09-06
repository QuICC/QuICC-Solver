/** 
 * @file TFMesher.cpp
 * @brief Source of the TF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/2D/TFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   TFMesher::TFMesher(const GridPurpose::Id purpose)
      : I2DMesher(purpose), mNx(-1), mNy(-1)
   {
   }

   void TFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& I = this->mDims.at(0);
      int& J = this->mDims.at(1);

      // Get standard dealiased FFT size
      this->mNx = Transform::Fft::Tools::dealiasCosFft(I+1);
      // Check for optimised FFT sizes
      this->mNx = Transform::Fft::Tools::optimizeFft(this->mNx);

      // Get standard dealiased FFT size
      this->mNy = Transform::Fft::Tools::dealiasMixedFft(J+1);
      // Check for optimised FFT sizes
      this->mNy = Transform::Fft::Tools::optimizeFft(this->mNy);
   }

   int TFMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int TFMesher::nPhys2D() const
   {
      return this->mNy;
   }

   int TFMesher::nSpec1D() const
   {
      const int& I = this->mDims.at(0);
      return I + 1;
   }

   int TFMesher::nSpec2D() const
   {
      const int& J = this->mDims.at(1);
      return J + 1;
   }

   int TFMesher::nDealias1D() const
   {
      return this->mNx;
   }

   int TFMesher::nDealias2D() const
   {
      return this->mNy/2 + 1;
   }

} // SpatialScheme
} // QuICC
