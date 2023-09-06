/** 
 * @file SLFMesher.cpp
 * @brief Source of the SLF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/SLFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   SLFMesher::SLFMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNr(-1), mNt(-1), mNp(-1)
   {
   }

   void SLFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& N = this->mDims.at(0);
      int& L = this->mDims.at(1);
      int& M = this->mDims.at(2);

      // Get dealiased Worland transform size
      this->mNr = Transform::Poly::Tools::dealias(N + 1 + 8);
      // Check for optimised FFT sizes
      this->mNr = Transform::Fft::Tools::optimizeFft(this->mNr);

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
         // Make space for r = 0, r = 1
         this->mNr += 2;
      }
   }

   int SLFMesher::nPhys1D() const
   {
      return this->mNr;
   }

   int SLFMesher::nPhys2D() const
   {
      return this->mNt;
   }

   int SLFMesher::nPhys3D() const
   {
      return this->mNp;
   }

   int SLFMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int SLFMesher::nSpec2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int SLFMesher::nSpec3D() const
   {
      const int& M = this->mDims.at(2);
      return M + 1;
   }

   int SLFMesher::nDealias1D() const
   {
      return this->mNr;
   }

   int SLFMesher::nDealias2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int SLFMesher::nDealias3D() const
   {
      return this->mNp/2 + 1;
   }

} // SpatialScheme
} // QuICC
