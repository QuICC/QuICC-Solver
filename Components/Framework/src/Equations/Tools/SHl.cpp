/**
 * @file SHl.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions with l spectral ordering
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/SHl.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> SHl::identifyIndexes(const Resolution& res, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      eigs.push_back(static_cast<MHDFloat>(res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)));

      return eigs;
   }


   int SHl::computeNMat(const Resolution& res) const
   {
      return res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

   void SHl::interpretTauN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHl::interpretGalerkinN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHl::interpretRhsN(ArrayI& rRhsCols, const Resolution& res) const
   {
      for(int l = 0; l < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); l++)
      {
         rRhsCols(l) = res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(l);
      }
   }

   void SHl::interpretSystemN(ArrayI& rSystemNs, const Resolution& res) const
   {
      // Python setup is sufficient
   }

}
}
}
