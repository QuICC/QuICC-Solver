/**
 * @file JLFMesher.cpp
 * @brief Source of the JLF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/JLFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"

namespace QuICC {

namespace SpatialScheme {

   JLFMesher::JLFMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNdealias(-1), mNr(-1), mNt(-1), mNp(-1)
   {
   }

   void JLFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      // Call base implementation
      IMesher::init(dims, options);

      int& N = this->mDims.at(0);
      int& L = this->mDims.at(1);
      int& M = this->mDims.at(2);
      int& nN_ = this->mNdealias;

      // Safety check
      if(L < M)
      {
         throw std::logic_error("Max harmonic degree L cannot be smaller than max harmonic order M");
      }

      // radial spectral resolution
      nN_ = N + 1;
      // radial grid resolution
      int nR_ = 15*(2*N + L + 1)/8;
      // Get dealiased Bessel transform size
      this->mNr = Transform::Poly::Tools::dealias(nR_);

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

   int JLFMesher::nPhys1D() const
   {
      return this->mNr;
   }

   int JLFMesher::nPhys2D() const
   {
      return this->mNt;
   }

   int JLFMesher::nPhys3D() const
   {
      return this->mNp;
   }

   int JLFMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int JLFMesher::nSpec2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int JLFMesher::nSpec3D() const
   {
      const int& M = this->mDims.at(2);
      return M + 1;
   }

   int JLFMesher::nDealias1D() const
   {
      return this->mNdealias;
   }

   int JLFMesher::nDealias2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int JLFMesher::nDealias3D() const
   {
      return this->mNp/2 + 1;
   }

} // SpatialScheme
} // QuICC
