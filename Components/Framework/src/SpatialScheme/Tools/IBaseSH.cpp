/** 
 * @file IBaseSH.cpp
 * @brief Source of the tools for generic + SH spatial schemes
 */

// System includes
//
#include <set>
#include <deque>
#include <queue>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/IBaseSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int IBaseSH::totalModes(const int n2D, const int n3D, Splitting::Locations::Id split)
   {
      int nModes = 0;
      if(split == Splitting::Locations::FIRST)
      {
         nModes = this->nHarmonics(n3D, n2D);
      }
      else if(split == Splitting::Locations::SECOND)
      {
         nModes = n2D;
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

   int IBaseSH::nM(const int l, const int nM)
   {
      return std::min(l+1, nM);
   }

   int IBaseSH::nL(const int m, const int nL)
   {
      return nL-m;
   }

   int IBaseSH::nHarmonics(const int nL, const int nM)
   {
      int tot = 0;

      // loop over all harmonics
      for(int l = 0; l < nL; l++)
      {
         // loop up to l, or if only a limited number of orders is required up to M
         for(int m = 0; m < this->nM(l, nM); m++)
         {
            tot++;
         }
      }

      return tot;
   }

   void IBaseSH::buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>&, const std::vector<int>&)
   {
      auto nL = n3D;
      auto nM = n2D;

      // Get restricted list of harmonic orders
      std::vector<int> binOrders;
      this->binMLLoad(binOrders, nL, nM, id.at(1), bins.at(1));

      // Create multimap of l,m modes
      std::multimap<int,int> lmap;
      for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
      {
         for(int l = *vIt; l < nL; l++)
         {
            lmap.insert(std::make_pair(l, *vIt));
         }
      }

      // Compute balanced harmonic load distribution
      std::vector<int> binLoad(bins.at(0), lmap.size()/bins.at(0));
      for(unsigned int i = 0; i < lmap.size() % bins.at(0); i++)
      {
         binLoad.at(i) += 1;
      }

      balanceOriginal(modes, lmap, id.at(0), binLoad);
      //auto keys = alternatingUpDown(lmap, binLoad.size());
      //balance(modes, lmap, id.at(0), binLoad, keys, 1);
   }

   void IBaseSH::balanceOriginal(std::multimap<int,int>& modes, std::multimap<int,int>& lmap, const int id, std::vector<int>& binLoad)
   {
      const int targetLoad = binLoad.at(id);
      auto& myLoad = binLoad.at(id);

      std::multimap<int,int> incomplete; // partially distributed modes: count -> l
      std::multimap<int,int>::iterator mapIt;
      std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> mapRange;
      int curId = 0;
      int counter = 0;
      int l = lmap.rbegin()->first;
      auto bins = binLoad.size();
      size_t maxLoop = lmap.size() + 5*bins;
      for(size_t loop = 0; loop < maxLoop; ++loop)
      {
         auto& curLoad = binLoad.at(curId);
         if(curLoad > 0)
         {
            // Extract modes from full list
            if(incomplete.size() == 0 || lmap.count(l) > static_cast<size_t>(incomplete.rbegin()->first))
            {
               // Get range of m for current l
               mapRange = lmap.equal_range(l);
               --l;
            }
            // Get modes from incomplete storage
            else
            {
               // Get range for largest incomplete set
               mapRange = lmap.equal_range(incomplete.rbegin()->second);

               // Delete it from incomplete map
               mapIt = incomplete.end();
               --mapIt;
               incomplete.erase(mapIt);
            }

            // Get number of harmonic orders in range
            int lCount = std::distance(mapRange.first, mapRange.second);

            // Extract partial m
            if(lCount > curLoad)
            {
               // Add information about remaining modes in incomplete list
               incomplete.insert(std::make_pair(lCount - curLoad, mapRange.first->first));
               // Create proper range
               mapRange.second = mapRange.first;
               std::advance(mapRange.second, curLoad);
               lCount = curLoad;
            }

            // Store the modes
            if(curId == id)
            {
               modes.insert(mapRange.first, mapRange.second);
            }

            // Delete used entries
            lmap.erase(mapRange.first, mapRange.second);

            // Substract used modes
            curLoad -= lCount;
         }

         // Load of this ID is zero or no modes left: break out of loop
         if(myLoad == 0 || lmap.size() == 0)
         {
            break;
         }

         // Update loop counter
         ++counter;

         // Fill by alternating directions
         curId = counter % bins;
         if(counter/bins % 2 == 1)
         {
            curId = bins - curId - 1;
         }

         assert(l >= 0);
      }

      assert(myLoad == 0);
      assert(modes.size() == static_cast<size_t>(targetLoad));
   }

   void IBaseSH::balanceAlternate(std::multimap<int,int>& modes, std::multimap<int,int>& lmap, const int id, std::vector<int>& binLoad)
   {
      const int targetLoad = binLoad.at(id);
      auto& myLoad = binLoad.at(id);

      // Check if there is an odd number of degrees
      int l = lmap.rbegin()->first;
      bool hasOddNumber = (l % 2 == 0);

      int counter = 0;
      for(int curId = 0; curId <= id; curId++)
      {
         auto& curLoad = binLoad.at(curId);
         while(curLoad > 0)
         {
            if(counter % 2 == 0)
            {
               l = lmap.begin()->first;
            }
            else
            {
               auto it = lmap.rbegin();
               // with odd numer of degrees, make pairs except for highest
               if(hasOddNumber && lmap.size() > lmap.count(it->first))
               {
                  std::advance(it, lmap.count(it->first));
               }
               l = it->first;
            }

            auto range = lmap.equal_range(l);
            int mcount = std::distance(range.first, range.second);
            range.second = range.first;
            std::advance(range.second, std::min(mcount, curLoad));

            // Store the modes
            if(curId == id)
            {
               modes.insert(range.first, range.second);
            }

            // Substract used modes
            auto usedModes = std::distance(range.first, range.second);
            assert(usedModes > 0);
            curLoad -= usedModes;

            // Delete used entries
            lmap.erase(range.first, range.second);

            if(lmap.count(l) == 0)
            {
               ++counter;
            }
         }
      }

      assert(myLoad == 0);
      assert(modes.size() == static_cast<size_t>(targetLoad));
   }

   std::list<int> IBaseSH::oldOrdering(std::multimap<int,int>& lmap, const std::size_t bins)
   {
      // Get list of harmonic degrees
      std::list<int> tmp;
      for(auto it = lmap.begin(); it != lmap.end(); )
      {
         tmp.push_back(it->first);
         std::advance(it, lmap.count(it->first));
      }
      tmp.reverse();
      std::vector<std::list<int> > bs(bins);
      int n = tmp.size();
      int sb = tmp.size()/bins + (tmp.size()%bins > 0);
      for(int i = 0; i < n; i++)
      {
         int bi = (i/sb)%bins;
         if(i%2 == 1)
         {
            bs.at(bi).push_back(tmp.front());
            tmp.pop_front();
         }
         else
         {
            bs.at(bi).push_back(tmp.back());
            tmp.pop_back();
         }
      }
      if(tmp.size() != 0)
      {
         throw std::logic_error("DID NOT USED ALL MODES");
      }

      std::list<int> keys;
      int shift = 0;
      for(int i = 0; i < n; i++)
      {
         int k = (i+shift)%bins;
         while(bs.at(k).size() == 0)
         {
            shift++;
            k = (i+shift)%bins;
         }
         keys.push_back(bs.at(k).front());
         bs.at(k).pop_front();
      }
      if(keys.size() != n)
      {
         throw std::logic_error("DID NOT reorder ALL MODES");
      }

      return keys;
   }

   std::list<int> IBaseSH::alternatingUpDown(std::multimap<int,int>& lmap, const std::size_t bins)
   {
      // Get list of harmonic degrees
      std::list<int> keys;
      for(auto it = lmap.begin(); it != lmap.end(); )
      {
         keys.push_back(it->first);
         std::advance(it, lmap.count(it->first));
      }

      // Reorder the keys
      keys.reverse();
      auto sorter = [=](const auto& a, const auto& b)
      {
         auto max = keys.front();
         auto idxA = (max - a);
         auto fwd = ((idxA/bins) % 2 == 0);
         auto va = (idxA/bins)*bins + fwd*(idxA%bins) + (!fwd)*(bins - idxA % bins);

         auto idxB = (max - b);
         fwd = ((idxB/bins) % 2 == 0);
         auto vb = (idxB/bins)*bins + fwd*(idxB%bins) + (!fwd)*(bins - idxB % bins);

         return (va < vb);
      };
      keys.sort(sorter);

      return keys;
   }

   void IBaseSH::balance(std::multimap<int,int>& modes, std::multimap<int,int>& lmap, const int id, std::vector<int>& binLoad, std::list<int>& keys, const int kind)
   {
      auto bins = binLoad.size();

      int splitPass;
      bool stopAtFirstSplit;

      // mode 1:
      // 1) Distribute full modes (ignoring big ones)
      // 2) Split up remaining modes
      if(kind == 1)
      {
         splitPass = 1;
         stopAtFirstSplit = false;
      }
      // mode 2:
      // 1) Distribute full modes until split is required
      // 2) Split up remaining modes
      else if(kind == 2)
      {
         splitPass = 1;
         stopAtFirstSplit = true;
      }
      // mode 3:
      // 1) Split up modes if required
      else if(kind == 3)
      {
         splitPass = 0;
         stopAtFirstSplit = false;
      }
      else
      {
         throw std::logic_error("Unknown split kind requested");
      }

      // First pass only uses full 3D modes
      // Second pass extracts partial 3D modes to fill bin
      for(int pass = 0; pass < splitPass+1; pass++)
      {
         int counter = 0;
         int curId = 0;
         for(auto it = keys.begin(); it != keys.end(); ++it)
         {
            bool usedFullMode = false;
            auto l = *it;

            // Try all bins
            for(unsigned int i = 0; i < bins; i++)
            {
               // Get bin ID and correspoonding work load
               curId = (counter + i) % bins;
               auto& curLoad = binLoad.at(curId);

               auto range = lmap.equal_range(-1);

               // process full mode
               if(lmap.count(l) <= curLoad)
               {
                  range = lmap.equal_range(l);
               }
               // split up mode
               else if(pass == splitPass)
               {
                  range = lmap.equal_range(l);
                  range.second = range.first;
                  std::advance(range.second, curLoad);
               }

               auto usedModes = std::distance(range.first, range.second);

               // Store modes and remove from list
               if(usedModes > 0)
               {
                  // Store the modes
                  if(curId == id)
                  {
                     modes.insert(range.first, range.second);
                  }

                  // Substract used modes
                  curLoad -= usedModes;

                  // Delete used entries
                  lmap.erase(range.first, range.second);

                  // If all modes used, delete key and break out
                  if(lmap.count(l) == 0)
                  {
                     usedFullMode = true;
                     it = keys.erase(it);
                     it--;
                     break;
                  }
               }
            }

            // stop-at-first-split mode is active and split is required
            if(stopAtFirstSplit && !usedFullMode)
            {
               break;
            }

            counter++;
         }
      }

      // Make sure work loads are full
      for(auto ld: binLoad)
      {
         if(ld != 0)
         {
            throw std::logic_error("Load is not fully distributed");
         }
      }

      // Make sure all modes have been distributed
      if(keys.size() != 0 || lmap.size() != 0)
      {
         throw std::logic_error("Not all modes distributed");
      }
   }

} // Tools
} // SpatialScheme
} // QuICC
