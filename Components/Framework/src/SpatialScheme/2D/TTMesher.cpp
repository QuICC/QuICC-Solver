/** 
 * @file TTMesher.cpp
 * @brief Source of the TT spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/2D/TTMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   TTMesher::TTMesher(const GridPurpose::Id purpose)
      : I2DMesher(purpose), mNx(-1), mNy(-1)
   {
   }

   void TTMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
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
      this->mNy = Transform::Fft::Tools::dealiasCosFft(J+1);
      // Check for optimised FFT sizes
      this->mNy = Transform::Fft::Tools::optimizeFft(this->mNy);
   }

   int TTMesher::nPhys1D() const
   {
      return this->mNx;
   }

   int TTMesher::nPhys2D() const
   {
      return this->mNy;
   }

   int TTMesher::nSpec1D() const
   {
      const int& I = this->mDims.at(0);
      return I + 1;
   }

   int TTMesher::nSpec2D() const
   {
      const int& J = this->mDims.at(1);
      return J + 1;
   }

   int TTMesher::nDealias1D() const
   {
      return this->mNx;
   }

   int TTMesher::nDealias2D() const
   {
      return this->mNy;
   }

} // SpatialScheme
} // QuICC
