/** 
 * @file FFFMesher.cpp
 * @brief Source of the FFF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/FFFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   FFFMesher::FFFMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNx(-1), mNy(-1), mNz(-1)
   {
   }

   void FFFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& I = this->mDims.at(0);
      int& J = this->mDims.at(1);
      int& K = this->mDims.at(2);

      // Get standard dealiased FFT size
      this->mNx = Transform::Fft::Tools::dealiasFft(I+1);
      // Check for optimised FFT sizes
      this->mNx = Transform::Fft::Tools::optimizeFft(this->mNx);

      // Get standard dealiased FFT size
      this->mNy = Transform::Fft::Tools::dealiasFft(J+1);
      // Check for optimised FFT sizes
      this->mNy = Transform::Fft::Tools::optimizeFft(this->mNy);

      // Get standard dealiased FFT size
      this->mNz = Transform::Fft::Tools::dealiasMixedFft(K+1);
      // Check for optimised FFT sizes
      this->mNz = Transform::Fft::Tools::optimizeFft(this->mNz);
   }

   int FFFMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int FFFMesher::nPhys2D() const
   {
      return this->mNy;
   }

   int FFFMesher::nPhys3D() const
   {
      return this->mNz;
   }

   int FFFMesher::nSpec1D() const
   {
      const int& I = this->mDims.at(0);
      return 2*(I + 1);
   }

   int FFFMesher::nSpec2D() const
   {
      const int& J = this->mDims.at(1);
      return 2*(J + 1);
   }

   int FFFMesher::nSpec3D() const
   {
      const int& K = this->mDims.at(2);
      return K + 1;
   }

   int FFFMesher::nDealias1D() const
   {
      return this->mNx;
   }

   int FFFMesher::nDealias2D() const
   {
      return this->mNy;
   }

   int FFFMesher::nDealias3D() const
   {
      return this->mNz/2 + 1;
   }

} // SpatialScheme
} // QuICC
