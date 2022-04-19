/**
 * @file SolverHasBc.hpp
 * @brief SolverHasBc ModelOperatorBoundary
 */

#ifndef QUICC_MODELOPERATORBOUNDARY_SOLVERHASBC_HPP
#define QUICC_MODELOPERATORBOUNDARY_SOLVERHASBC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/ModelOperatorBoundary/IRegisterId.hpp"

namespace QuICC {

namespace ModelOperatorBoundary {

   /**
    * @brief SolverHasBc ModelOperatorBoundary
    */
   class SolverHasBc: public IRegisterId<SolverHasBc>
   {
      public:
         /**
          * @brief Constructor
          */
         SolverHasBc();

         friend class IRegisterId<SolverHasBc>;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();

         /**
          * @brief Formatted name
          */
         static std::string sFormatted();
   };

}
}

#endif // QUICC_MODELOPERATORBOUNDARY_SOLVERHASBC_HPP
