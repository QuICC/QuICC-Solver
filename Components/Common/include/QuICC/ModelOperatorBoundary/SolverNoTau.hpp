/**
 * @file SolverNoTau.hpp
 * @brief SolverNoTau ModelOperatorBoundary
 */

#ifndef QUICC_MODELOPERATORBOUNDARY_SOLVERNOTAU_HPP
#define QUICC_MODELOPERATORBOUNDARY_SOLVERNOTAU_HPP

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
    * @brief SolverNoTau ModelOperatorBoundary
    */
   class SolverNoTau: public IRegisterId<SolverNoTau>
   {
      public:
         /**
          * @brief Constructor
          */
         SolverNoTau();

         friend class IRegisterId<SolverNoTau>;

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

#endif // QUICC_MODELOPERATORBOUNDARY_SOLVERNOTAU_HPP
