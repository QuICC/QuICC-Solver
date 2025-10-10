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
#include <fstream>

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

namespace details
{
void writeDot(std::string filename, const std::vector<int>& xnodes, const std::vector<int>& xadj, const std::vector<int>& adjncy, const std::vector<int>& part, const int stage)
{
   int n = xadj.size()-1;

   std::ofstream dot(filename);

   dot << "strict graph {" << std::endl;
   dot << "node [colorscheme=set19]" << std::endl;
   for(int s = 0; s < xnodes.size() - 1; s++)
   {
      for(int node = 0; node < xnodes[s+1] - xnodes[s]; node++)
      {
         dot << "\"s" << s << "_" << node << "\"" << " [style = filled, color=" << part[node] + 1 << "]" << std::endl;
      }
   }
   for(int node = 0; node < n; node++)
   {
      for(int i = xadj[node]; i < xadj[node+1]; i++)
      {
         int nLeft = node;
         int nRight = adjncy[i];
         int sLeft = -42;
         int sRight = -42;
         for(int k = 1; k < xnodes.size(); k++)
         {
            if(sLeft < 0 && nLeft - xnodes[k] < 0)
            {
               sLeft = k-1;
               nLeft -= xnodes[k-1];
            }
            if(sRight < 0 && nRight - xnodes[k] < 0)
            {
               sRight = k-1;
               nRight -= xnodes[k-1];
            }
         }

         if(stage < 0 || ((sLeft == stage && sRight == stage + 1) || (sLeft == stage + 1 && sRight == stage)))
         {
            dot << "\"s" << sLeft << "_" << nLeft << "\" -- \"s" << sRight << "_" << nRight << "\"" << std::endl;
         }
      }
   }

  // Close the file
  dot << "}" << std::endl;
  dot.close();
}

void writePartition(const std::string filename, const std::vector<int>& part)
{
   std::ofstream partition(filename);

   for(auto&& c: part)
   {
      partition << c << std::endl;
   }

   partition.close();
}

void writeMetis(const std::string filename, const std::vector<int>& xadj, const std::vector<int>& adjncy, const std::vector<int>& vwgt, const std::vector<int>& adjcwgt)
{
   std::ofstream metis(filename);

   metis << xadj.size() - 1 << " " << adjncy.size()/2 << " " << 11 << std::endl;

   for(int node = 0; node < xadj.size()-1; node++)
   {
      metis << vwgt[node];
      for(int i = xadj[node]; i < xadj[node+1]; i++)
      {
         metis << " " << adjncy[i] + 1 << " " << adjcwgt[i];
      }
      metis << std::endl;
   }

   metis.close();
}
}

}
}
