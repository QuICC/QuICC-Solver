/** 
 * @file SH.cpp
 * @brief Source of the Spherical Harmonics tools for spatial schemes
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/Tools/SH.hpp"

// Project includes
//

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int SH::nM(const int l, const int nM)
   {
      return std::min(l+1, nM);
   }

   int SH::nL(const int m, const int nL)
   {
      return nL-m;
   }

   int SH::nHarmonics(const int nL, const int nM)
   {
      int tot = 0;

      // loop over all harmonics
      for(int l = 0; l < nL; l++)
      {
         // loop up to l, or if only a limited number of orders is required up to M
         for(int m = 0; m < SH::nM(l, nM); m++)
         {
            tot++;
         }
      }

      return tot;
   }

   void SH::buildMMap(std::multimap<int,int>& harmonics, const int nL, const int nM)
   {
      // List of the harmonics
      harmonics.clear();

      // loop over all harmonic orders
      for(int m = 0; m < nM; m++)
      {
         // loop over the available harmonic degrees
         for(int l = 0; l < SH::nL(m, nL); l++)
         {
            harmonics.insert(std::make_pair(m, l+m));
         }
      }
   }

   void SH::buildLMap(std::multimap<int,int>& harmonics, const int nL, const int nM)
   {
      // List of the harmonics
      harmonics.clear();

      // loop over all harmonics degrees
      for(int l = 0; l < nL; l++)
      {
         // loop up the available hramonic orders
         for(int m = 0; m < SH::nM(l, nM); m++)
         {
            harmonics.insert(std::make_pair(l, m));
         }
      }
   }

   void SH::buildLHSortedMap(std::multimap<int,int>& harmonics, const int nL, const int nM, const int h0, const int nH)
   {
      // Make sure the list is empty at start
      harmonics.clear();

      // Treat the even and odd number of orders carefully
      int maxL;
      int maxM = nM;
      if(nM % 2 == 0)
      {
         maxL = maxM;
      } else
      {
         maxL = maxM - 1;
      }

      // Create first part of list as pairs
      std::queue<std::pair<int,int> >  modeQueue;
      for(int l = 0; l < maxL/2; ++l)
      {
         for(int m = 0; m < SH::nM(l,nM); ++m)
         {
            modeQueue.push(std::make_pair(l, m));
         }

         for(int m = 0; m < SH::nM(maxL-l-1, nM); ++m)
         {
            modeQueue.push(std::make_pair(maxL-l-1, m));
         }
      }

      // Add in the missing ones
      for(int l = maxL; l < nL; ++l)
      {
         for(int m = 0; m < SH::nM(l, nM); ++m)
         {
            modeQueue.push(std::make_pair(l, m));
         }
      }

      // Remove unused modes from queue
      for(int h = 0; h < h0; ++h)
      {
         modeQueue.pop();
      }

      // Fill harmonics map with requested harmonics
      for(int h = 0; h < nH; h++)
      {
         harmonics.insert(modeQueue.front());
         modeQueue.pop();
      }
   }

   void SH::buildMHSortedMap(std::multimap<int,int>& harmonics, const int nL, const int nM, const int h0, const int nH)
   {
      // Make sure the list is empty at start
      harmonics.clear();

      // Treat the even and odd number of orders carefully
      int maxM;
      int minM;
      if(nM % 2 == 0)
      {
         minM = 0;
         maxM = nM/2;
      } else
      {
         minM = 1;
         maxM = (nM+1)/2;
      }

      // Create first part of list as pairs
      std::queue<std::pair<int,int> >  modeQueue;
      for(int m = minM; m < maxM; ++m)
      {
         for(int l = 0; l < SH::nL(m,nL); ++l)
         {
            modeQueue.push(std::make_pair(m, l + m));
         }

         for(int l = 0; l < SH::nL(nM-1-m, nL); ++l)
         {
            modeQueue.push(std::make_pair(nM-1-m, l+nM-1-m));
         }
      }

      // Add in the first droped ones
      for(int m = 0; m < minM; ++m)
      {
         for(int l = 0; l < SH::nL(m, nL); ++l)
         {
            modeQueue.push(std::make_pair(m, l+m));
         }
      }

      // Remove unused modes from queue
      for(int h=0; h < h0; ++h)
      {
         modeQueue.pop();
      }

      // Fill harmonics map with requested harmonics
      for(int h = 0; h < nH; h++)
      {
         harmonics.insert(modeQueue.front());
         modeQueue.pop();
      }
   }

   void SH::binMLLoad(std::vector<int>& binOrders, const int nL, const int nM, const int id, const int bins)
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

   void SH::initMLLoad(std::deque<int>& loadList, std::vector<int>& binLoad, std::queue<int>& optimal, const int nL, const int nM, const int bins)
   {
      // Initialise total load
      double totalLoad = 0.0;

      // Compute total load and create list of loads
      for(int i = 0; i < nM; ++i)
      {
         // for a given harmonic order m there are (maxL + 1 - m) degrees
         loadList.push_back(nL - i);

         // increment total load
         totalLoad += static_cast<double>(nL - i);
      }

      // Compute the ideal load per CPU
      double ideal = totalLoad / static_cast<double>(bins);

      // Deal with part that can't be evenly shared 
      int extra = static_cast<int>(totalLoad) % bins;

      // Compute optimal load for all bins
      for(int i = 0; i < bins; ++i)
      {
         // add rounded ideal value + possible suboptimal part
         optimal.push(static_cast<int>(ideal) + (i < extra));
      }

      // Initialise the load per part
      for(int i = 0; i < bins; ++i)
      {
         binLoad.push_back(0);
      }
   }

   void SH::combineMPairs(std::deque<int>& list, std::multimap<int, int>& regular, const int nM, const int bins)
   {
      // Aiming for the same number of modes per part, compute the number of slots avaiable per part
      int slots = static_cast<int>(std::ceil(static_cast<double>(nM)/static_cast<double>(bins)));

      // Algorithm is based on pairing modes together, odd number of slots requires to be careful
      bool oddBins = slots % 2;
      // Take into account the possibly uneven share per part
      bool needCare = (nM % bins) + oddBins;

      // Number of pairs to build
      int nPairs;

      // Pairs a build from the two opposite ends of list, start index of reverse index
      int endIdx;

      // Need careful distribution ?
      if(needCare)
      {
         // Odd or even number of slots require different algorithm
         if(oddBins)
         {
            nPairs = slots/2 - 1;
         } else
         {
            nPairs = slots/2 - 2;
         }

         // Start index from other end of list
         endIdx = 2*nPairs*bins - 1;

      // Perfect split we simply build all the pairs
      } else
      {
         nPairs = slots/2;

         endIdx = nM - 1;
      }
      
      // Counter for the number of elements removed
      int k = 0;

      // Loop over the pairs
      for(int j = 0; j < nPairs; ++j)
      {
         // Loop over the number bins
         for(int i = 0; i < bins; ++i)
         {
            // Add load from beginning
            regular.insert(std::make_pair(i, list.at(k)));

            // ... and add load from end
            regular.insert(std::make_pair(i, list.at(endIdx - k)));

            // Increment counter
            k++;
         }
      }

      // Remove the used loads from the list
      list.erase(list.begin(), list.begin() + 2*k);
   }

   void SH::fillMPairsBins(std::deque<int>& list, std::multimap<int, int>& regular, std::vector<int>& load, const int bins)
   {
      // Obviously only required if any load is left
      if(list.size() > 0)
      {
         // This part of algorithm requires 2 modes per part
         if(list.size() >= static_cast<unsigned int>(2*bins))
         {
            // Aiming for the same number of modes per part, compute the number of slots avaiable per part
            int slots = static_cast<int>(std::ceil(static_cast<double>(list.size())/static_cast<double>(bins)));
            // Algorithm is based on pairing modes together, odd number of slots requires to be careful
            bool oddBins =  slots % 2;

            // Add an additional mode to the assigned loads (in all cases there are at least 2 spaces left)
            for(int j = 0; j < bins; ++j)
            {
               regular.insert(std::make_pair(j, list.at(j)));
            }

            // Remove just added loads from list
            list.erase(list.begin(), list.begin() + bins);

            // Special treatment is required if the number of rows is even
            if(! oddBins)
            {
               // Additional shifted mode (the idea is to get a similar outcome than in the odd case)
               //  Get "upper" half of the bins
               int uPart = bins/2 + (bins % 2);
               //  Get "lower" half of the bins
               int lPart = bins - uPart;

               // Assign loads to upper part
               for(int j = 0; j < uPart; ++j)
               {
                  regular.insert(std::make_pair(j, list.at(bins - 1 - 2*j)));
               }

               // Assign loads to lower part
               for(int j = 0; j < lPart; ++j)
               {
                  regular.insert(std::make_pair(uPart + j, list.at(bins - 2  - 2*j)));
               }

               // Remove just added loads from list
               list.erase(list.begin(), list.begin() + bins);
            }
         }

         // If the load list is empty now. Make sure all get at least one
         if(regular.size() == 0)
         {
            for(int i = 0; i < bins; i++)
            {
               //
               regular.insert(std::make_pair(i,list.at(i)));
               // Update related sum
               load.at(i) += list.at(i);
            }

            // Remove just added loads from list
            list.erase(list.begin(), list.begin() + bins);
         }
      }
   }

   void SH::fillMRestBins(std::deque<int>& list, std::multimap<int, int>& regular, std::vector<int>& load, std::queue<int>& optimal, const int bins)
   {
      // Obviously only required if any load is left
      if(list.size() > 0)
      {
         // Distribute the remaining loads
         unsigned int leftLoads = list.size();

         // Progression counters
         unsigned int idx = 0;
         int current;
         int curSum;

         // Algorithm is not perfect, so add margin to fit
         int margin = 0;

         // Loop until there is no remaining load
         while(idx < leftLoads)
         {
            // Reset selected part
            current = -1;
            curSum = -1;

            // Search for best place to assign load
            for(int j = 0; j < bins; ++j)
            {
               // Check if position is suitable
               if(load.at(j) +  list.at(idx) <= optimal.front() + margin && (load.at(j) +  list.at(idx) > curSum))
               {
                  // Set current best fit information
                  current = j;
                  curSum = load.at(j) +  list.at(idx);

                  // Reset margin
                  margin = 0;
               }
            }

            // If suitable position has be found assign load to it
            if(current != -1)
            {
               // Assign load
               regular.insert(std::make_pair(current,list.at(idx)));
               // Update related sum
               load.at(current) += list.at(idx);
               // increment counter
               idx++;

               // If part reached optimal load, remove it from load queue
               if(load.at(current) == optimal.front())
               {
                  optimal.pop();
               }
            // If no suitable position has be found, increase fit margin
            } else
            {
               margin++;
            }
         }
      }
   }

   void SH::convertLoadToOrders(std::multimap<int, int>& regular, const int nL)
   {
      for(auto it = regular.begin(); it != regular.end(); it++)
      {
         it->second = nL - it->second;
      }
   }

   void SH::fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D)
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
         bwd1D.push_back(ArrayI(nB1D - idx3D(k)));

         // Fill array with indexes
         for(int i = 0; i < bwd1D.at(k).size(); i++)
         {
            bwd1D.at(k)(i) = idx3D(k) + i;
         }
      }
   }
}
}
}
