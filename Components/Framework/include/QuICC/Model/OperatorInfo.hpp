/**
 * @file OperatorInfo.hpp
 * @brief Information describing operators from backend
 */

#ifndef QUICC_MODEL_OPERATORINFO_HPP
#define QUICC_MODEL_OPERATORINFO_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Information describing operators from backend
    */
   struct OperatorInfo
   {
         /**
          * @brief Constructor
          */
         OperatorInfo(const int n)
            : tauN(n), galN(n), galShift(n,3), rhsCols(n), sysN(n)
         {};

         /**
          * @brief Destructor
          */
         ~OperatorInfo() = default;

         /**
          * @brief Tau dimensions
          */
         ArrayI tauN;

         /**
          * @brief Galerkin dimensions
          */
         ArrayI galN;

         /**
          * @brief Shift when using galerkin basis
          */
         MatrixI galShift;

         /**
          * @brief Number of cols in RHS
          */
         ArrayI rhsCols;

         /**
          * @brief Full system size
          */
         ArrayI sysN;
   };

}
}

#endif // QUICC_MODEL_OPERATORINFO_HPP
