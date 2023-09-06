/**
 * @file ImExRKCB3a.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3a
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExRKCB3a.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   int ImExRKCB3a::steps() const
   {
      return 3;
   }

   // Scheme order
   int ImExRKCB3a::order() const
   {
      return 3;
   }

   // Scheme has embedded lower order scheme?
   bool ImExRKCB3a::hasEmbedded() const
   {
      return false;
   }

   // Name of the scheme
   std::string ImExRKCB3a::name() const
   {
      return "ImExRKCB3a";
   }

   ImExRKCB3a::ImExRKCB3a()
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

   MHDFloat ImExRKCB3a::aIm(const int i, const int j) const
   {
      return ImExRKCB3a::mAIm(i,j);
   }

   MHDFloat ImExRKCB3a::bIm(const int i) const
   {
      return ImExRKCB3a::mBIm(i);
   }

   MHDFloat ImExRKCB3a::aEx(const int i, const int j) const
   {
      return ImExRKCB3a::mAEx(i,j);
   }

   MHDFloat ImExRKCB3a::bEx(const int i) const
   {
      return ImExRKCB3a::mBEx(i);
   }

   MHDFloat ImExRKCB3a::cEx(const int i) const
   {
      return ImExRKCB3a::mCEx(i);
   }

   MHDFloat ImExRKCB3a::bImErr(const int i) const
   {
      return ImExRKCB3a::mBImErr(i);
   }

   MHDFloat ImExRKCB3a::bExErr(const int i) const
   {
      return ImExRKCB3a::mBExErr(i);
   }

   void ImExRKCB3a::init()
   {
      MHDFloat c2 = (27.0 + std::pow(2187.0 - 1458.0*std::sqrt(2.0), 1./3.) + 9.0*std::pow(3.0 + 2.0*std::sqrt(2.0), 1./3.))/54.0;
      MHDFloat c3 = c2/(6.0*std::pow(c2,2) - 3.0*c2 + 1.0);
      MHDFloat b2 = (3.0*c2 - 1.0)/(6.0*std::pow(c2,2));
      MHDFloat b3 = (6.0*std::pow(c2,2) - 3.0*c2 + 1.0)/(6.0*std::pow(c2,2));

      // Initialize implicit a factors
      ImExRKCB3a::mAIm(1,1) = c2;
      ImExRKCB3a::mAIm(2,2) = (1./6. - b2*std::pow(c2,2) - b3*c2*c3)/(b3*(c3 - c2));
      ImExRKCB3a::mAIm(2,1) = ImExRKCB3a::mAIm(2,2) - c3;

      // Initialize implicit b factors
      ImExRKCB3a::mBIm(1) = b2;
      ImExRKCB3a::mBIm(2) = b3;

      // Initialize explicit a factors
      ImExRKCB3a::mAEx(1,0) = c2;
      ImExRKCB3a::mAEx(2,1) = c3;

      // Initialize explicit b factors
      ImExRKCB3a::mBEx(1) = b2;
      ImExRKCB3a::mBEx(2) = b3;

      // Initialize explicit c factors
      ImExRKCB3a::mCEx(1) = c2;
      ImExRKCB3a::mCEx(2) = c3;
   }

}
}
