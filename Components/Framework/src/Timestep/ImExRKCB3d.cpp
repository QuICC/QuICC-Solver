/**
 * @file ImExRKCB3d.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3d
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExRKCB3d.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   int ImExRKCB3d::steps() const
   {
      return 4;
   }

   // Scheme order
   int ImExRKCB3d::order() const
   {
      return 3;
   }

   // Scheme has embedded lower order scheme?
   bool ImExRKCB3d::hasEmbedded() const
   {
      return true;
   }

   // Name of the scheme
   std::string ImExRKCB3d::name() const
   {
      return "ImExRKCB3d";
   }

   ImExRKCB3d::ImExRKCB3d()
      : IImExRKCBScheme(),
      mAIm(Eigen::Array<MHDFloat,4,4>::Zero()),
      mBIm(Eigen::Array<MHDFloat,4,1>::Zero()),
      mAEx(Eigen::Array<MHDFloat,4,4>::Zero()),
      mBEx(Eigen::Array<MHDFloat,4,1>::Zero()),
      mCEx(Eigen::Array<MHDFloat,4,1>::Zero()),
      mBImErr(Eigen::Array<MHDFloat,4,1>::Zero()),
      mBExErr(Eigen::Array<MHDFloat,4,1>::Zero())
   {
      this->init();
   }

   MHDFloat ImExRKCB3d::aIm(const int i, const int j) const
   {
      return ImExRKCB3d::mAIm(i,j);
   }

   MHDFloat ImExRKCB3d::bIm(const int i) const
   {
      return ImExRKCB3d::mBIm(i);
   }

   MHDFloat ImExRKCB3d::aEx(const int i, const int j) const
   {
      return ImExRKCB3d::mAEx(i,j);
   }

   MHDFloat ImExRKCB3d::bEx(const int i) const
   {
      return ImExRKCB3d::mBEx(i);
   }

   MHDFloat ImExRKCB3d::cEx(const int i) const
   {
      return ImExRKCB3d::mCEx(i);
   }

   MHDFloat ImExRKCB3d::bImErr(const int i) const
   {
      return ImExRKCB3d::mBImErr(i);
   }

   MHDFloat ImExRKCB3d::bExErr(const int i) const
   {
      return ImExRKCB3d::mBExErr(i);
   }

   void ImExRKCB3d::init()
   {
      // Initialize implicit b factors
      ImExRKCB3d::mBIm(0) = 0.0;
      ImExRKCB3d::mBIm(1) = 355931813527.0/1014712533305.0;
      ImExRKCB3d::mBIm(2) = 709215176366.0/1093407543385.0;
      ImExRKCB3d::mBIm(3) = 755675305.0/1258355728177.0;

      // Initialize implicit a factors
      ImExRKCB3d::mAIm(1,0) = 0.0;
      ImExRKCB3d::mAIm(1,1) = 418884414754.0/469594081263.0;
      ImExRKCB3d::mAIm(2,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAIm(2,1) = -304881946513433262434901.0/718520734375438559540570.0;
      ImExRKCB3d::mAIm(2,2) = 684872032315.0/962089110311.0;
      ImExRKCB3d::mAIm(3,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAIm(3,1) = ImExRKCB3d::mBIm(1);
      ImExRKCB3d::mAIm(3,2) = ImExRKCB3d::mBIm(2);
      ImExRKCB3d::mAIm(3,3) = ImExRKCB3d::mBIm(3);

      // Initialize explicit a factors
      ImExRKCB3d::mAEx(1,0) = 418884414754.0/469594081263.0;
      ImExRKCB3d::mAEx(2,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAEx(2,1) = 214744852859.0/746833870870.0;
      ImExRKCB3d::mAEx(3,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAEx(3,1) = ImExRKCB3d::mBIm(1);
      ImExRKCB3d::mAEx(3,2) = 658780719778.0/1014712533305.0;

      // Initialize explicit b factors
      ImExRKCB3d::mBEx(0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mBEx(1) = ImExRKCB3d::mBIm(1);
      ImExRKCB3d::mBEx(2) = ImExRKCB3d::mBIm(2);
      ImExRKCB3d::mBEx(3) = ImExRKCB3d::mBIm(3);

      // Initialize explicit c factors
      ImExRKCB3d::mCEx(1) = ImExRKCB3d::mAEx(1,0);
      ImExRKCB3d::mCEx(2) = ImExRKCB3d::mAEx(2,1);
      ImExRKCB3d::mCEx(3) = 1.0;

      // Initialize implicit embedded b factors
      ImExRKCB3d::mBImErr(0) = 0.0;
      ImExRKCB3d::mBImErr(1) = 226763370689.0/646029759300.0;
      ImExRKCB3d::mBImErr(2) = 1496839794860.0/2307829317197.0;
      ImExRKCB3d::mBImErr(3) = 353416193.0/889746336234.0;

      // Initialize explicit embedded b factors
      ImExRKCB3d::mBExErr(0) = 1226988580973.0/2455716303853.0;
      ImExRKCB3d::mBExErr(1) = 0.0;
      ImExRKCB3d::mBExErr(2) = 827818615.0/1665592077861.0;
      ImExRKCB3d::mBExErr(3) = 317137569431.0/634456480332.0;
   }

}
}
