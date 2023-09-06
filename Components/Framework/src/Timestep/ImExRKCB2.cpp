/**
 * @file ImExRKCB2.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 2
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExRKCB2.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   int ImExRKCB2::steps() const
   {
      return 3;
   }

   // Scheme order
   int ImExRKCB2::order() const
   {
      return 2;
   }

   // Scheme has embedded lower order scheme?
   bool ImExRKCB2::hasEmbedded() const
   {
      return true;
   }

   // Name of the scheme
   std::string ImExRKCB2::name() const
   {
      return "ImExRKCB2";
   }

   ImExRKCB2::ImExRKCB2()
      : IImExRKCBScheme(),
      mAIm(Eigen::Array<MHDFloat,3,3>::Zero()),
      mBIm(Eigen::Array<MHDFloat,3,1>::Zero()),
      mAEx(Eigen::Array<MHDFloat,3,3>::Zero()),
      mBEx(Eigen::Array<MHDFloat,3,1>::Zero()),
      mCEx(Eigen::Array<MHDFloat,3,1>::Zero()),
      mBImErr(Eigen::Array<MHDFloat,3,1>::Zero()),
      mBExErr(Eigen::Array<MHDFloat,3,1>::Zero())
   {
      this->init();
   }

   MHDFloat ImExRKCB2::aIm(const int i, const int j) const
   {
      return this->mAIm(i,j);
   }

   MHDFloat ImExRKCB2::bIm(const int i) const
   {
      return this->mBIm(i);
   }

   MHDFloat ImExRKCB2::aEx(const int i, const int j) const
   {
      return this->mAEx(i,j);
   }

   MHDFloat ImExRKCB2::bEx(const int i) const
   {
      return this->mBEx(i);
   }

   MHDFloat ImExRKCB2::cEx(const int i) const
   {
      return this->mCEx(i);
   }

   MHDFloat ImExRKCB2::bImErr(const int i) const
   {
      return this->mBImErr(i);
   }

   MHDFloat ImExRKCB2::bExErr(const int i) const
   {
      return this->mBExErr(i);
   }

   void ImExRKCB2::init()
   {
      // Initialize implicit a factors
      this->mAIm(1,1) = 2./5.;
      this->mAIm(2,1) = 5./6.;
      this->mAIm(2,2) = 1./6.;

      // Initialize implicit b factors
      this->mBIm(1) = 5./6.;
      this->mBIm(2) = 1./6.;

      // Initialize explicit a factors
      this->mAEx(1,0) = 2./5.;
      this->mAEx(2,1) = 1.;

      // Initialize explicit b factors
      this->mBEx(1) = 5./6.;
      this->mBEx(2) = 1./6.;

      // Initialize explicit c factors
      this->mCEx(1) = 2./5.;
      this->mCEx(2) = 1.;

      // Initialize implicit embedded b factors
      this->mBImErr(1) = 4./5.;
      this->mBImErr(2) = 1./5.;

      // Initialize explicit embedded b factors
      this->mBExErr(1) = 4./5.;
      this->mBExErr(2) = 1./5.;
   }

}
}
