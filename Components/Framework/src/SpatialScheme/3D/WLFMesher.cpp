/**
 * @file WLFMesher.cpp
 * @brief Source of the WLF spatial scheme mesher
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFMesher.hpp"
#include "QuICC/Transform/Poly/Tools.hpp"
#include "QuICC/Transform/Fft/Tools.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Fft.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/Transform/Setup/Trapezoidal.hpp"
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"
#include "QuICC/SpatialScheme/Tools/TrapezoidalSH.hpp"
#include "QuICC/Debug/DebuggerMacro.h"

namespace QuICC {

namespace SpatialScheme {

   WLFMesher::WLFMesher(const GridPurpose::Id purpose)
      : IMesher(purpose), mNdealias(-1), mNr(-1), mNt(-1), mNp(-1)
   {
   }

   void WLFMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
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

      const auto& opt = options.at(0);

      int nR_ = 0;
      // Poly algorithm setup with triangular radial truncation
      if(
            std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end() &&
            std::find(opt.begin(), opt.end(), Transform::Setup::Triangular::id()) != opt.end())
      {
         DebuggerMacro_msg("Setting up triangular radial truncation", 1);

         Tools::TriangularSH t;

         // Check if resolution is compatible with triangular truncation
         if(!t.isOptimal(N + 1,L))
         {
            throw std::logic_error("Triangular truncation requires N - L/2 = " + std::to_string(t.min()-1));
         }

         // radial spectral resolution
         nN_ = N + 1;

         // radial grid resolution from l = 0
         const int Nat0 = t.truncationBwd(nN_+7,0,0) - 1;
         nR_ = (2*Nat0 + 1)/2;
         // ... check against l = L requirements
         const int NatL = t.truncationBwd(nN_ + 7, 0, L) - 1;
         nR_ = std::max(nR_, (2*NatL + L + 1)/2);
      }
      // Poly algorithm setup with trapezoidal radial truncation
      else if(
            std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end() &&
            std::find(opt.begin(), opt.end(), Transform::Setup::Trapezoidal::id()) != opt.end())
      {
         DebuggerMacro_msg("Setting up trapezoidal radial truncation", 1);

         Tools::TrapezoidalSH t;

         // Check if resolution is compatible with trapezoidal truncation
         if(!t.isOptimal(N + 1,L))
         {
            throw std::logic_error("Trapezoidal truncation requires N - L/2 >= " + std::to_string(t.min()-1));
         }

         // radial spectral resolution
         nN_ = N + 1;

         // radial grid resolution from l = 0
         const int Nat0 = t.truncationBwd(nN_ + 7, 0, 0) - 1;
         nR_ = (2*Nat0 + 1)/2;
         // ... check against l = L requirements
         const int NatL = t.truncationBwd(nN_ + 7, 0, L) - 1;
         nR_ = std::max(nR_, (2*NatL + L + 1)/2);
      }
      // Poly algorithm setup with uniform radial truncation
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::GaussianQuadrature::id()) != opt.end())
      {
         DebuggerMacro_msg("Setting up uniform radial truncation", 1);

         // radial spectral resolution
         nN_ = N + 1;

         // radial grid resolution
         nR_ = N + L/2 + 8;
      }
      // FFT algorithm setup with uniform radial truncation
      else if(std::find(opt.begin(), opt.end(), Transform::Setup::Fft::id()) != opt.end())
      {
         DebuggerMacro_msg("Setting up uniform radial truncation with FFT Worland transform" , 1);

         // radial spectral resolution
         nN_ = N + 1;

         // radial grid resolution
         nR_ = N + L/2 + 4;
      }

      // Get dealiased Worland transform size
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

   int WLFMesher::nPhys1D() const
   {
      return this->mNr;
   }

   int WLFMesher::nPhys2D() const
   {
      return this->mNt;
   }

   int WLFMesher::nPhys3D() const
   {
      return this->mNp;
   }

   int WLFMesher::nSpec1D() const
   {
      const int& N = this->mDims.at(0);
      return N + 1;
   }

   int WLFMesher::nSpec2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int WLFMesher::nSpec3D() const
   {
      const int& M = this->mDims.at(2);
      return M + 1;
   }

   int WLFMesher::nDealias1D() const
   {
      return this->mNdealias;
   }

   int WLFMesher::nDealias2D() const
   {
      const int& L = this->mDims.at(1);
      return L + 1;
   }

   int WLFMesher::nDealias3D() const
   {
      return this->mNp/2 + 1;
   }

} // SpatialScheme
} // QuICC
