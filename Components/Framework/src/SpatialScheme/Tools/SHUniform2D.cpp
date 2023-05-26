/** 
 * @file SHUniform2D.cpp
 * @brief Source of the tools for Spherical Harmonics + uniform + uniform for spatial schemes
 */

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/Tools/SHUniform2D.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   SHUniform2D::SHUniform2D(const int nL)
      : mNl(nL)
   {
   }

   int SHUniform2D::truncationFwd(const int nN, const int k)
   {
      return nN;
   }

   int SHUniform2D::truncationBwd(const int nN, const int k)
   {
      return nN - k;
   }

   int SHUniform2D::index(const int i, const int k)
   {
      return i + k;
   }

   bool SHUniform2D::isOptimal(const int nN, const int maxL)
   {
      return true;
   }

   int SHUniform2D::totalModes(const int n2D, const int n3D, Splitting::Locations::Id split)
   {
      int nModes = 0;
      if(split == Splitting::Locations::FIRST)
      {
         nModes = n2D;
      }
      else if(split == Splitting::Locations::SECOND)
      {
         nModes = n3D;
      }
      else if(split == Splitting::Locations::BOTH)
      {
         nModes = n3D;
      }
      else
      {
         throw std::logic_error("Splitting location is not implemented");
      }

      return nModes;
   }

   void SHUniform2D::buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>&, const std::vector<int>&)
   {
      auto nM = n3D;
      auto nL = this->mNl;

      int n0;
      int nN;

      // Get restricted list of harmonics
      std::vector<int> binOrders;
      this->binMLLoad(binOrders, nL, nM, id.at(1), bins.at(1));

      // Compute second distribution
      auto tot = this->totalModes(n2D, n3D, Splitting::Locations::FIRST);
      this->balancedSplit(n0, nN, tot, bins.at(0), id.at(0), false);

      // Fill with correct modes
      for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
      {
         for(int r = 0; r < nN; r++)
         {
            modes.insert(std::make_pair(*vIt, n0 + r));
         }
      }
   }

} // Tools
} // SpatialScheme
} // QuICC
