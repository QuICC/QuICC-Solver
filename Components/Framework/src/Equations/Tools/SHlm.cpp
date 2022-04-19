/**
 * @file SHlm.cpp
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
#include "QuICC/Equations/Tools/SHlm.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> SHlm::identifyIndexes(const Resolution& res, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = res.cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      eigs.push_back(static_cast<MHDFloat>(mode(2)));
      eigs.push_back(static_cast<MHDFloat>(mode(3)));

      return eigs;
   }

   int SHlm::computeNMat(const Resolution& res) const
   {
      int nMat = 0;

      for(int i = 0; i < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return nMat;
   }

   void SHlm::interpretTauN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHlm::interpretGalerkinN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHlm::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHlm::interpretSystemN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

}
}
}
