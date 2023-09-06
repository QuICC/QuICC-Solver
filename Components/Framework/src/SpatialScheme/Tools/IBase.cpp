/** 
 * @file IBase.cpp
 * @brief Source of the base tools for spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   void IBase::buildMap(std::multimap<int,int>& modes, const int k0, const int kN, const std::vector<int>& j0, const std::vector<int>& jN, const int c0, const int cN)
   {
      // Counter
      int c = 0;

      // Loop over third dimension
      for(int k = 0; k < kN; k++)
      {
         // Loop over second dimension
         for(int j = 0; j < jN.at(k); j++)
         {
            // Check for first mode
            if(c >= c0)
            {
               if(c >= cN)
               {
                  break;
               }
               else
               {
                  modes.insert(std::make_pair(k0 + k, j0.at(k) + j));
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

   void IBase::fillIndexes2D3D(std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::multimap<int,int>& modes)
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

   void IBase::fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D)
   {
      // Make full list of indexes for first dimension
      for(int k = 0; k < idx3D.size(); k++)
      {
         // Create storage for indexes of forward transform
         fwd1D.push_back(ArrayI(this->truncationFwd(nF1D, idx3D(k))));

         // Fill array with indexes of forward transform
         for(int i = 0; i < fwd1D.at(k).size(); i++)
         {
            fwd1D.at(k)(i) = i;
         }

         // Create storage for indexes of backward transform
         bwd1D.push_back(ArrayI(this->truncationBwd(nB1D, idx3D(k))));

         // Fill array with indexes of backward transform
         for(int i = 0; i < bwd1D.at(k).size(); i++)
         {
            bwd1D.at(k)(i) = this->index(i, idx3D(k));
         }
      }
   }

   void IBase::binMLLoad(std::vector<int>& binOrders, const int nL, const int nM, const int id, const int bins)
   {
      // Distribute the heterogeneous harmonic order load by going back and forth
      for(int i = 0; i < nM; ++i)
      {
         // Is right bin while going left -> right
         bool isLeftToRight = ((i/bins) % 2 == 0 && (i % bins) == id);

         // Is right bin while going right -> left
         bool isRightToLeft = ((i/bins) % 2 == 1 && (i % bins) == (bins - id - 1));

         if(isLeftToRight || isRightToLeft)
         {
            // for a given harmonic order m there are (maxL + 1 - m) degrees
            binOrders.push_back(i);
         }
      }
   }

   void IBase::balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id, const bool allowEmpty)
   {
      // Avoid splitting with zero elements
      if(tot < parts && !allowEmpty)
      {
         throw std::logic_error("Number of parts is bigger than total!");
      }

      // Compute part assigned to id
      if(parts > 1)
      {
         nN = 0;
         n0 = 0;
         for(int i = 0; i < tot; i++)
         {
            if(i % parts == id)
            {
               nN++;
            }
            else if(i % parts < id)
            {
               n0++;
            }
         }

      // Single part, use total
      } else if(parts == 1)
      {
         n0 = 0;
         nN = tot;

      // Can't split into less than 1 part
      } else
      {
         throw std::logic_error("Number of parts < 1!");
      }
   }

} // Tools
} // SpatialScheme
} // QuICC
