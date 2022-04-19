/** 
 * @file Regular.cpp
 * @brief Source of the tools for regular spatial schemes
 */

// System includes
//
#include <set>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/Tools/Regular.hpp"

// Project includes
//

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   void Regular::buildMap(std::multimap<int,int>& modes, const int k0, const int kN, const ArrayI& j0, const ArrayI& jN, const int c0, const int cN)
   {
      // Counter
      int c = 0;

      // Loop over third dimension
      for(int k = 0; k < kN; k++)
      {
         // Loop over second dimension
         for(int j = 0; j < jN(k); j++)
         {
            // Check for first mode
            if(c >= c0)
            {
               if(c >= cN)
               {
                  break;
               } else
               {
                  modes.insert(std::make_pair(k0 + k,j0(k) + j));
               }
            }
            c++;
         }
         if(c >= cN)
         {
            break;
         }
      }
   }

   void Regular::fillIndexes2D3D(std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::multimap<int,int>& modes)
   {
      // Set to extract the 3D indexes
      std::set<int>  filter;

      // Loop over all modes
      for(auto mapIt = modes.cbegin(); mapIt != modes.cend(); ++mapIt)
      {
         filter.insert(mapIt->first);
      }

      // Set third dimension
      idx3D.resize(filter.size());

      // Make full list of index in third dimension
      auto setIt = filter.begin();
      for(int k = 0; k < idx3D.size(); k++)
      {
         idx3D(k) = *setIt;
         ++setIt;
      }

      // Make full list of indexes for second dimension
      for(int k = 0; k < idx3D.size(); k++)
      {
         // Create storage for indexes
         idx2D.push_back(ArrayI(modes.count(idx3D(k))));

         // Get range
         auto mapRange = modes.equal_range(idx3D(k));

         // Loop over range
         int j = 0;
         for(auto mapIt = mapRange.first; mapIt != mapRange.second; ++mapIt)
         {
            idx2D.at(k)(j) = mapIt->second;
            j++;
         }
      }
   }

   void Regular::fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D)
   {
      // Make full list of indexes for first dimension
      for(int k = 0; k < idx3D.size(); k++)
      {
         // Create storage for indexes
         fwd1D.push_back(ArrayI(nF1D));

         // Fill array with indexes
         for(int i = 0; i < fwd1D.at(k).size(); i++)
         {
            fwd1D.at(k)(i) = i;
         }

         // Create storage for indexes
         bwd1D.push_back(ArrayI(nB1D));

         // Fill array with indexes
         for(int i = 0; i < bwd1D.at(k).size(); i++)
         {
            bwd1D.at(k)(i) = i;
         }
      }
   }

}
}
}
