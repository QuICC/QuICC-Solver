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
          * @param ignoreExtreme Ignore end members in factorization
          */
         static void filterFactors(std::list<int>& cpuFactors, const int nFactors, const int nCpu, const bool ignoreExtreme);

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
          * @param ignoreExtreme Ignore end members in factorization
          */
         static bool confirmFactors(const std::vector<int>& factors, const int nCpu,  const bool ignoreExtreme);

         /**
          * @brief Maximum number of factors to test
          */
         static const int mcMaxDecompositions = 3;

         /**
          * @brief Constructor
          */
         SplittingTools() = default;

         /**
          * @brief Destructor
          */
         ~SplittingTools() = default;
   };

}
}

#endif // QUICC_PARALLEL_SPLITTINGTOOLS_HPP
