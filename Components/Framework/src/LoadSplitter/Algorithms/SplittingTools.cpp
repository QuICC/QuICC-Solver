/** 
 * @file SplittingTools.cpp
 * @brief Source of the base of some tools used for the splitting
 */

// System includes
//
#include <algorithm>
#include <set>
#include <map>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   void SplittingTools::factorizeNCpu(std::list<int>& cpuFactors, const int nFactors, const int nCpu)
   {
      // Select factorisation algorithm depending on number of factors
      if(nFactors == 1)
      {
         // Add factor
         cpuFactors.push_back(nCpu);

      // Factorise CPUs into two groups
      } else if(nFactors == 2)
      {
         // Get the maximum factor
         int factor = static_cast<int>(std::sqrt(nCpu));

         // Compute smaller factors
         while(factor > 0)
         {
            if(nCpu % factor == 0)
            {
               // Add factor
               cpuFactors.push_back(factor);

               // Add nCpu / factor
               cpuFactors.push_back(nCpu/factor);

               // Add reversed splitting order
               if(factor != nCpu/factor)
               {
                  // Add factor
                  cpuFactors.push_back(nCpu/factor);

                  // Add nCpu / factor
                  cpuFactors.push_back(factor);
               }
            }
            --factor;
         }
      } else
      {
         throw std::logic_error("No factorisation algorithm available for requested factors!");
      }
   }

   void SplittingTools::filterFactors(std::list<int>& cpuFactors, const int nFactors, const int nCpu)
   {
      // Get iterator through known factors
      std::list<int>::iterator  it = cpuFactors.begin();
      std::list<int>::iterator  itF;

      ArrayI factors(nFactors);
      bool suitable;

      // Loop over all known factors
      while(it != cpuFactors.end())
      {
         // Extract factors to test
         itF = it;
         for(int i = 0; i < nFactors; i++)
         {
            factors(i) = *itF;
            itF++;
         }

         // Test if factors are usable splitting factors
         suitable = SplittingTools::confirmFactors(factors, nCpu);

         // Move to the next set of factors
         if(suitable)
         {
            std::advance(it, nFactors);

         // Factors are not usable
         } else
         {
            // Erase the unusable factors
            for(int i = 0; i < nFactors; i++)
            {
               it = cpuFactors.erase(it);
            }
         }
      }
   }

   bool SplittingTools::confirmFactors(const ArrayI& factors, const int nCpu)
   {
      // Loop over all factors
      for(int i =0; i < factors.size(); i++)
      {
         // We don't want the extrem cases (no splitting in one direction)
         if(factors(i) == 1 && nCpu > 1)
         {
            return false;
         }
      }

      // In all other cases accept the splitting
      return true;
   }

   void SplittingTools::balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id, const bool allowEmpty)
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

   void SplittingTools::splitMapped(const std::multimap<int, int>& mapped, ArrayI &rIdx, const int id)
   {
      // resize the array of indexes
      rIdx.resize(mapped.count(id));

      // Get range of indexes for given id
      auto range = mapped.equal_range(id);

      // Put indexes into a set to be sure to get them sorted
      std::set<int>  sorter;
      for(auto it = range.first; it != range.second; ++it)
      {
         sorter.insert(it->second);
      }

      // Extract the ordered indexes from set and store in output array
      int i = 0;
      for(auto setIt = sorter.begin(); setIt != sorter.end(); ++setIt, ++i)
      {
         rIdx(i) = *setIt;
      }
   }

   int SplittingTools::groupId(const ArrayI& factors, const int i, const int id)
   {
      // Assert on index of requested factor
      assert(i < factors.size());

      switch(i)
      {
         case(0):
            return id % factors(0);
            break;
         case(1):
            return id / factors(0);
            break;
         case(2):
            return id / (factors(0)*factors(1));
            break;
         default:
            throw std::logic_error("Unknown input provided!");
      }
   }

}
}
