/** 
 * @file TFFMesher.cpp
 * @brief Source of the TFF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/TFFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   TFFMesher::TFFMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNx(-1), mNy(-1), mNz(-1)
   {
   }

   void TFFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& X = this->mDims.at(0);
      int& Y = this->mDims.at(1);
      int& Z = this->mDims.at(2);

      // Get standard dealiased cosine FFT size
      this->mNx = Transform::Fft::Tools::dealiasCosFft(X + 1);
      // Check for optimised FFT sizes
      this->mNx = Transform::Fft::Tools::optimizeFft(this->mNx);

      // Get standard dealiased FFT size
      this->mNy = Transform::Fft::Tools::dealiasFft(Y + 1);
      // Check for optimised FFT sizes
      this->mNy = Transform::Fft::Tools::optimizeFft(this->mNy);

      // Get standard dealiased FFT size
      this->mNz = Transform::Fft::Tools::dealiasMixedFft(Z + 1);
      // Check for optimised FFT sizes
      this->mNz = Transform::Fft::Tools::optimizeFft(this->mNz);
   }

   int TFFMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int TFFMesher::nPhys2D() const
   {
      return this->mNy;
   }

   int TFFMesher::nPhys3D() const
   {
      return this->mNz;
   }

   int TFFMesher::nSpec1D() const
   {
      const int& X = this->mDims.at(0);
      return X + 1;
   }

   int TFFMesher::nSpec2D() const
   {
      const int& Y = this->mDims.at(1);
      return Y + 1;
   }

   int TFFMesher::nSpec3D() const
   {
      const int& Z = this->mDims.at(2);
      return Z + 1;
   }

   int TFFMesher::nDealias1D() const
   {
      return this->mNx;
   }

   int TFFMesher::nDealias2D() const
   {
      return this->mNy;
   }

   int TFFMesher::nDealias3D() const
   {
      return this->mNz/2 + 1;
   }

} // SpatialScheme
} // QuICC
