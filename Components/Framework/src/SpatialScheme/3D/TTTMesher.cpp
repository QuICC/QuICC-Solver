/** 
 * @file TTTMesher.cpp
 * @brief Source of the TTT spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/TTTMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   TTTMesher::TTTMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNx(-1), mNy(-1), mNz(-1)
   {
   }

   void TTTMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
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

      // Get standard dealiased FFT size
      this->mNy = Transform::Fft::Tools::dealiasCosFft(J+1);
      // Check for optimised FFT sizes
      this->mNy = Transform::Fft::Tools::optimizeFft(this->mNy);

      // Get standard dealiased FFT size
      this->mNz = Transform::Fft::Tools::dealiasCosFft(K+1);
      // Check for optimised FFT sizes
      this->mNz = Transform::Fft::Tools::optimizeFft(this->mNz);
   }

   int TTTMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int TTTMesher::nPhys2D() const
   {
      return this->mNy;
   }

   int TTTMesher::nPhys3D() const
   {
      return this->mNz;
   }

   int TTTMesher::nSpec1D() const
   {
      const int& I = this->mDims.at(0);
      return I + 1;
   }

   int TTTMesher::nSpec2D() const
   {
      const int& J = this->mDims.at(1);
      return J + 1;
   }

   int TTTMesher::nSpec3D() const
   {
      const int& K = this->mDims.at(2);
      return K + 1;
   }

   int TTTMesher::nDealias1D() const
   {
      return this->mNx;
   }

   int TTTMesher::nDealias2D() const
   {
      return this->mNy;
   }

   int TTTMesher::nDealias3D() const
   {
      return this->mNz;
   }

} // SpatialScheme
} // QuICC
