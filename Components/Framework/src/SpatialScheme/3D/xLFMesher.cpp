/**
 * @file xLFMesher.cpp
 * @brief Source of the xLF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/xLFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   xLFMesher::xLFMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNdealias(-1), mNr(-1), mNt(-1), mNp(-1)
   {
   }

   void xLFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& L = this->mDims.at(1);
      int& M = this->mDims.at(2);

      // Get dealiased associated Legendre transform size
      this->mNt = Transform::Poly::Tools::dealias(L + 1);

      // Get standard dealiased FFT size
      this->mNp = Transform::Fft::Tools::dealiasMixedFft(M + 1);
      // Check for optimised FFT sizes
      this->mNp = Transform::Fft::Tools::optimizeFft(this->mNp);

      // Modify grid size for visualiation
      if(this->mPurpose == GridPurpose::VISUALIZATION)
      {
         // Make space for theta = 0 and  theta = pi
         this->mNt += 2;
      }
   }

   int xLFMesher::nPhys1D() const
   {
      return this->mNr;
   }

   int xLFMesher::nPhys2D() const
   {
      return this->mNt;
   }

   int xLFMesher::nPhys3D() const
   {
      return this->mNp;
   }

   int xLFMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int xLFMesher::nSpec2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int xLFMesher::nSpec3D() const
   {
      const int& M = this->mDims.at(2);
      return M + 1;
   }

   int xLFMesher::nDealias1D() const
   {
      return this->mNdealias;
   }

   int xLFMesher::nDealias2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int xLFMesher::nDealias3D() const
   {
      return this->mNp/2 + 1;
   }

} // SpatialScheme
} // QuICC
