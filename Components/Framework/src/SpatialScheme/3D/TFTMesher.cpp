/** 
 * @file TFTMesher.cpp
 * @brief Source of the TFT spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/TFTMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   TFTMesher::TFTMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNx(-1), mNy(-1), mNz(-1)
   {
   }

   void TFTMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& I = this->mDims.at(0);
      int& J = this->mDims.at(1);
      int& K = this->mDims.at(2);

      // Get standard dealiased FFT size
      this->mNx = Transform::Fft::Tools::dealiasCosFft(I+1);
      // Check for optimised FFT sizes
      this->mNx = Transform::Fft::Tools::optimizeFft(this->mNx);

      // Get mixed dealiased FFT size
      this->mNy = Transform::Fft::Tools::dealiasMixedFft(J+1);
      // Check for optimised FFT sizes
      this->mNy = Transform::Fft::Tools::optimizeFft(this->mNy);

      // Get standard dealiased FFT size
      this->mNz = Transform::Fft::Tools::dealiasCosFft(K+1);
      // Check for optimised FFT sizes
      this->mNz = Transform::Fft::Tools::optimizeFft(this->mNz);
   }

   int TFTMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int TFTMesher::nPhys2D() const
   {
      return this->mNy;
   }

   int TFTMesher::nPhys3D() const
   {
      return this->mNz;
   }

   int TFTMesher::nSpec1D() const
   {
      const int& I = this->mDims.at(0);
      return I + 1;
   }

   int TFTMesher::nSpec2D() const
   {
      const int& J = this->mDims.at(1);
      return J + 1;
   }

   int TFTMesher::nSpec3D() const
   {
      const int& K = this->mDims.at(2);
      return K + 1;
   }

   int TFTMesher::nDealias1D() const
   {
      return this->mNx;
   }

   int TFTMesher::nDealias2D() const
   {
      return this->mNy/2 + 1;
   }

   int TFTMesher::nDealias3D() const
   {
      return this->mNz;
   }

} // SpatialScheme
} // QuICC
