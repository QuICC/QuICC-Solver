/**
 * @file ISphericalCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalCflWrapper.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"

namespace QuICC {

namespace Diagnostics {

   ISphericalCflWrapper::ISphericalCflWrapper(const MHDFloat courant)
      : ICflWrapper(courant)
   {
   }

   void ISphericalCflWrapper::initMesh(const std::vector<Array>& mesh)
   {
      if(this->mFields.size() == 0 || this->mFields.begin()->second == nullptr)
      {
         throw std::logic_error("CFL calculation not setup correctly");
      }

      // Compute the mesh spacings
      this->mMeshSpacings.reserve(3);

      // Storage for radial grid
      this->mMeshSpacings.push_back(mesh.at(0));

      // Storage for radial grid
      this->mMeshSpacings.push_back(Array(mesh.at(0).size()));

      // Storage for horizontal average grid spacing
      this->mMeshSpacings.push_back(Array(mesh.at(0).size()));

      Array& r = this->mMeshSpacings.at(0);
      Array& dr = this->mMeshSpacings.at(1);
      Array& r_ll1 = this->mMeshSpacings.at(2);

      // Compute grid spacings
      for(int j = 0; j < r.size(); ++j)
      {
         // Get internal points grid spacing
         if(j > 0 && j < r.size() - 1)
         {
            dr(j) = std::min(std::abs(r(j) - r(j-1)), std::abs(r(j) - r(j+1)));
         }
         // Get left endpoint grid spacing
         else if(j > 0)
         {
            dr(j) = std::abs(r(j) - r(j-1));

         }
         // Get right endpoint grid spacing
         else
         {
            dr(j) = std::abs(r(j) - r(j+1));
         }

         // Compute average horizontal grid spacing
         MHDFloat effL = this->effectiveMaxL(r(j));
         r_ll1(j) = r(j)/std::sqrt(effL*(effL + 1.0));
         if(r_ll1(j) == 0)
         {
            r_ll1(j) = std::numeric_limits<MHDFloat>::max();
         }
      }
   }

}
}
