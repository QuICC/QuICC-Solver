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
#include <numeric>
#include <functional>

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

      }
      // Factorise CPUs into two groups
      else if(nFactors == 2)
      {
         // Get the maximum factor
         int factor = static_cast<int>(std::sqrt(nCpu));

         // Compute smaller factors
         while(factor > 0 && cpuFactors.size() < static_cast<std::size_t>(nFactors*SplittingTools::mcMaxDecompositions))
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
      }
      else
      {
         throw std::logic_error("No factorisation algorithm available for requested factors!");
      }
   }

   void SplittingTools::filterFactors(std::list<int>& cpuFactors, const int nFactors, const int nCpu, const bool ignoreExtreme)
   {
      // Get iterator through known factors
      std::list<int>::iterator  it = cpuFactors.begin();
      std::list<int>::iterator  itF;

      std::vector<int> factors(nFactors);
      bool suitable;

      // Loop over all known factors
      while(it != cpuFactors.end())
      {
         // Extract factors to test
         itF = it;
         for(int i = 0; i < nFactors; i++)
         {
            factors.at(i) = *itF;
            itF++;
         }

         // Test if factors are usable splitting factors
         suitable = SplittingTools::confirmFactors(factors, nCpu, ignoreExtreme);

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

   bool SplittingTools::confirmFactors(const std::vector<int>& factors, const int nCpu, const bool ignoreExtreme)
   {
      bool status = true;

      // Loop over all factors
      for(std::size_t i = 0; i < factors.size(); i++)
      {
         // Check product of factors is nCpu
         status = status && (std::accumulate(factors.begin(), factors.end(), 1, std::multiplies<int>()) == nCpu);

         // We don't want the extrem cases (no splitting in one direction)
         if(ignoreExtreme)
         {
            status = status && !(factors.at(i) == 1 && nCpu > 1);
         }
      }

      return status;
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
