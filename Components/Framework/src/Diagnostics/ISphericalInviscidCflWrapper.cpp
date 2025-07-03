/**
 * @file ISphericalInviscidCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalInviscidCflWrapper.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"

namespace QuICC {

namespace Diagnostics {

   ISphericalInviscidCflWrapper::ISphericalInviscidCflWrapper(const MHDFloat courant)
      : ICflWrapper(courant)
   {
   }

   void ISphericalInviscidCflWrapper::initMesh(const std::vector<Array>& mesh)
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

      // Compute magnetic grid
      int nB = 3*this->mFields.begin()->second->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)/4 + 1;
      Internal::Array igrid, iweights;
      Polynomial::Quadrature::WorlandRule wquad;
      wquad.computeQuadrature(igrid, iweights, nB);
      Array rB = igrid.cast<MHDFloat>();
      Array drB = Array(rB.size());

      // Compute B grid spacings
      for(int j = 0; j < rB.size(); ++j)
      {
         // Get internal points grid spacing
         if(j > 0 && j < rB.size() - 1)
         {
            drB(j) = std::min(std::abs(rB(j) - rB(j-1)), std::abs(rB(j) - rB(j+1)));
         }
         // Get left endpoint grid spacing
         else if(j > 0)
         {
            drB(j) = std::abs(rB(j) - rB(j-1));

         }
         // Get right endpoint grid spacing
         else
         {
            drB(j) = std::abs(rB(j) - rB(j+1));
         }
      }

      // Compute grid spacings
      int jB = 0;
      for(int j = 0; j < r.size(); ++j)
      {
         // Grid position is smaller than B grid
         if(r(j) <= rB(jB) || jB == rB.size()-1)
         {
            dr(j) = drB(jB);
         }
         // Move to next
         else if(r(j) <= rB(jB+1))
         {
            jB++;
            dr(j) = drB(jB);
         }
         // At end of B grid
         else if(jB == rB.size()-1)
         {
            dr(j) = drB(jB);
         }
         else
         {
            throw std::logic_error("Mapping to radial B grid failed");
         }

         // Compute average horizontal grid spacing
         MHDFloat effL = this->effectiveMaxL(r(j));
         r_ll1(j) = r(j)/std::sqrt(effL*(effL + 1.0));
      }
   }

}
}
