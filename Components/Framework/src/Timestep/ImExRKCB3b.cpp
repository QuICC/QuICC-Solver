/**
 * @file ImExRKCB3b.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3b
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExRKCB3b.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   int ImExRKCB3b::steps() const
   {
      return 4;
   }

   // Scheme order
   int ImExRKCB3b::order() const
   {
      return 3;
   }

   // Scheme has embedded lower order scheme?
   bool ImExRKCB3b::hasEmbedded() const
   {
      return false;
   }

   // Name of the scheme
   std::string ImExRKCB3b::name() const
   {
      return "ImExRKCB3b";
   }

   ImExRKCB3b::ImExRKCB3b()
      : IImExRKCBScheme(),
      mAIm(Eigen::Array<MHDFloat,4,4>::Zero()),
      mBIm(Eigen::Array<MHDFloat,4,1>::Zero()),
      mAEx(Eigen::Array<MHDFloat,4,4>::Zero()),
      mBEx(Eigen::Array<MHDFloat,4,1>::Zero()),
      mCEx(Eigen::Array<MHDFloat,4,1>::Zero()),
      mBImErr(Eigen::Array<MHDFloat,4,1>::Zero()),
      mBExErr(Eigen::Array<MHDFloat,4,1>::Zero())
   {
   }

   ImExRKCB3b::~ImExRKCB3b()
   {
   }

   MHDFloat ImExRKCB3b::aIm(const int i, const int j) const
   {
      return ImExRKCB3b::mAIm(i,j);
   }

   MHDFloat ImExRKCB3b::bIm(const int i) const
   {
      return ImExRKCB3b::mBIm(i);
   }

   MHDFloat ImExRKCB3b::aEx(const int i, const int j) const
   {
      return ImExRKCB3b::mAEx(i,j);
   }

   MHDFloat ImExRKCB3b::bEx(const int i) const
   {
      return ImExRKCB3b::mBEx(i);
   }

   MHDFloat ImExRKCB3b::cEx(const int i) const
   {
      return ImExRKCB3b::mCEx(i);
   }

   MHDFloat ImExRKCB3b::bImErr(const int i) const
   {
      return ImExRKCB3b::mBImErr(i);
   }

   MHDFloat ImExRKCB3b::bExErr(const int i) const
   {
      return ImExRKCB3b::mBExErr(i);
   }

   void ImExRKCB3b::init()
   {
      MHDFloat g = 1./2. + std::sqrt(3.0)/6.0;
      MHDFloat c3 = 1./2. - std::sqrt(3.0)/6.0;
      MHDFloat c4 = 1./2. + std::sqrt(3.0)/6.0;

      // Initialize implicit a factors
      ImExRKCB3b::mAIm(1,1) = g;
      ImExRKCB3b::mAIm(2,1) = -std::sqrt(3.0)/3.0;
      ImExRKCB3b::mAIm(2,2) = g;
      ImExRKCB3b::mAIm(3,3) = g;

      // Initialize implicit b factors
      ImExRKCB3b::mBIm(2) = 1./2.;
      ImExRKCB3b::mBIm(3) = 1./2.;

      // Initialize explicit a factors
      ImExRKCB3b::mAEx(1,0) = g;
      ImExRKCB3b::mAEx(2,1) = c3;
      ImExRKCB3b::mAEx(3,2) = c4;

      // Initialize explicit b factors
      ImExRKCB3b::mBEx(2) = 1./2.;
      ImExRKCB3b::mBEx(3) = 1./2.;

      // Initialize explicit c factors
      ImExRKCB3b::mCEx(1) = g;
      ImExRKCB3b::mCEx(2) = c3;
      ImExRKCB3b::mCEx(3) = c4;
   }

}
}
