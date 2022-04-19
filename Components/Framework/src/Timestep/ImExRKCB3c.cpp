/**
 * @file ImExRKCB3c.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3c
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExRKCB3c.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   int ImExRKCB3c::steps() const
   {
      return 4;
   }

   // Scheme order
   int ImExRKCB3c::order() const
   {
      return 3;
   }

   // Scheme has embedded lower order scheme?
   bool ImExRKCB3c::hasEmbedded() const
   {
      return true;
   }

   // Name of the scheme
   std::string ImExRKCB3c::name() const
   {
      return "ImExRKCB3c";
   }

   ImExRKCB3c::ImExRKCB3c()
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

   ImExRKCB3c::~ImExRKCB3c()
   {
   }

   MHDFloat ImExRKCB3c::aIm(const int i, const int j) const
   {
      return ImExRKCB3c::mAIm(i,j);
   }

   MHDFloat ImExRKCB3c::bIm(const int i) const
   {
      return ImExRKCB3c::mBIm(i);
   }

   MHDFloat ImExRKCB3c::aEx(const int i, const int j) const
   {
      return ImExRKCB3c::mAEx(i,j);
   }

   MHDFloat ImExRKCB3c::bEx(const int i) const
   {
      return ImExRKCB3c::mBEx(i);
   }

   MHDFloat ImExRKCB3c::cEx(const int i) const
   {
      return ImExRKCB3c::mCEx(i);
   }

   MHDFloat ImExRKCB3c::bImErr(const int i) const
   {
      return ImExRKCB3c::mBImErr(i);
   }

   MHDFloat ImExRKCB3c::bExErr(const int i) const
   {
      return ImExRKCB3c::mBExErr(i);
   }

   void ImExRKCB3c::init()
   {
      // Initialize implicit b factors
      ImExRKCB3c::mBIm(0) = 0.0;
      ImExRKCB3c::mBIm(1) = 673488652607.0/2334033219546.0;
      ImExRKCB3c::mBIm(2) = 493801219040.0/853653026979.0;
      ImExRKCB3c::mBIm(3) = 184814777513.0/1389668723319.0;

      // Initialize implicit a factors
      ImExRKCB3c::mAIm(1,0) = 0.0;
      ImExRKCB3c::mAIm(1,1) = 3375509829940.0/4525919076317.0;
      ImExRKCB3c::mAIm(2,0) = ImExRKCB3c::mBIm(0);
      ImExRKCB3c::mAIm(2,1) = -11712383888607531889907.0/32694570495602105556248.0;
      ImExRKCB3c::mAIm(2,2) = 566138307881.0/912153721139.0;
      ImExRKCB3c::mAIm(3,0) = ImExRKCB3c::mBIm(0);
      ImExRKCB3c::mAIm(3,1) = ImExRKCB3c::mBIm(1);
      ImExRKCB3c::mAIm(3,2) = ImExRKCB3c::mBIm(2);
      ImExRKCB3c::mAIm(3,3) = ImExRKCB3c::mBIm(3);

      // Initialize explicit a factors
      ImExRKCB3c::mAEx(1,0) = 3375509829940.0/4525919076317.0;
      ImExRKCB3c::mAEx(2,0) = ImExRKCB3c::mBIm(0);
      ImExRKCB3c::mAEx(2,1) = 272778623835.0/1039454778728.0;
      ImExRKCB3c::mAEx(3,0) = ImExRKCB3c::mBIm(0);
      ImExRKCB3c::mAEx(3,1) = ImExRKCB3c::mBIm(1);
      ImExRKCB3c::mAEx(3,2) = 1660544566939.0/2334033219546.0;

      // Initialize explicit b factors
      ImExRKCB3c::mBEx(0) = ImExRKCB3c::mBIm(0);
      ImExRKCB3c::mBEx(1) = ImExRKCB3c::mBIm(1);
      ImExRKCB3c::mBEx(2) = ImExRKCB3c::mBIm(2);
      ImExRKCB3c::mBEx(3) = ImExRKCB3c::mBIm(3);

      // Initialize explicit c factors
      ImExRKCB3c::mCEx(1) = ImExRKCB3c::mAEx(1,0);
      ImExRKCB3c::mCEx(2) = ImExRKCB3c::mAEx(2,1);
      ImExRKCB3c::mCEx(3) = 1.0;

      // Initialize implicit embedded b factors
      ImExRKCB3c::mBImErr(0) = 0.0;
      ImExRKCB3c::mBImErr(1) = 366319659506.0/1093160237145.0;
      ImExRKCB3c::mBImErr(2) = 270096253287.0/480244073137.0;
      ImExRKCB3c::mBImErr(3) = 104228367309.0/1017021570740.0;

      // Initialize explicit embedded b factors
      ImExRKCB3c::mBExErr(0) = 449556814708.0/1155810555193.0;
      ImExRKCB3c::mBExErr(1) = 0.0;
      ImExRKCB3c::mBExErr(2) = 210901428686.0/1400818478499.0;
      ImExRKCB3c::mBExErr(3) = 480175564215.0/1042748212601.0;
   }

}
}
