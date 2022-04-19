/** 
 * @file SplittingTools.hpp
 * @brief Base of the implementation of the load splitting algorithms
 */

#ifndef QUICC_PARALLEL_SPLITTINGTOOLS_HPP
#define QUICC_PARALLEL_SPLITTINGTOOLS_HPP

// Configuration includes
//

// System includes
//
#include <list>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Base of the implementation of the load splitting algorithms
    */
   class SplittingTools
   {
      public:
         /**
          * @brief Decompose \f$N_{cpu}\f$ into possible factors
          *
          * @param cpuFactors Storage for the CPU factors
          * @param nFactors   Number of factors in factorisation
          * @param nCpu       Number of CPUs
          */
         static void factorizeNCpu(std::list<int>& cpuFactors, const int nFactors, const int nCpu);

         /**
          * @brief Filter out the unusable factors
          *
          * @param cpuFactors CPU factors
          * @param nFactors   Number of factors in factorisation
          * @param nCpu       Number of CPUs
          */
         static void filterFactors(std::list<int>& cpuFactors, const int nFactors, const int nCpu);

         /**
          * @brief Compute a simple balanced split of the elements with regard to the given number of parts
          *
          * @param n0      Output start index
          * @param nN      Output number of indexes
          * @param tot     Total number of indexes
          * @param parts   Number of parts to split total into
          * @param id      ID of the CPU
          */
         static void balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id);

         /**
          * @brief Extract splitting from mapped indexes
          *
          * @param mapped  Mapped indexes to extract from
          * @param rIdx    Output storage to put the indexes in
          * @param id      ID of the CPU
          */
         static void splitMapped(const std::multimap<int, int>& mapped, ArrayI &rIdx, const int id);

         /**
          * @brief Convert ID to \f$F_{i}\f$ groupd ID
          *
          * @param      factors CPU factors
          * @param i    ID of the factorisation group
          * @param id   ID of CPU
          */
         static int groupId(const ArrayI& factors, const int i, const int id);
         
      protected:

      private:
         /**
          * @brief Test the splitting factors compatibility
          *
          * @param factors CPU factors
          * @param nCpu    Number of CPUs
          */
         static bool confirmFactors(const ArrayI& factors, const int nCpu);

         /**
          * @brief Constructor
          */
         SplittingTools();

         /**
          * @brief Destructor
          */
         ~SplittingTools();
   };

}
}

#endif // QUICC_PARALLEL_SPLITTINGTOOLS_HPP
