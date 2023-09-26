/** 
 * @file IBaseAllL.cpp
 * @brief Source of the tools for uniform + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/IBaseAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   IBaseAllL::IBaseAllL(const int min)
      : IBaseSH(min)
   {
   }

   int IBaseAllL::truncationFwd(const int nN, const int j, const int k)
   {
      return this->truncationBwd(nN, j, k);
   }

   void IBaseAllL::buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>&, const std::vector<int>&)
   {
      auto nL = n2D;
      auto nM = n3D;

      // Get restricted list of harmonic orders
      std::vector<int> binOrders;
      this->binMLLoad(binOrders, nL, nM, id.at(1), bins.at(1));

      // Create multimap of l,m modes
      std::multimap<int,int> mmap;
      for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
      {
         for(int l = *vIt; l < nL; l++)
         {
            mmap.insert(std::make_pair(*vIt, l));
         }
      }

      int i = 0;
      for(auto mapIt = mmap.begin(); mapIt != mmap.end(); ++mapIt)
      {
         // add mode
         if(i % bins.at(0) == id.at(0))
         {
            modes.insert(*mapIt);
         }

         ++i;
      }
   }

} // Tools
} // SpatialScheme
} // QuICC
