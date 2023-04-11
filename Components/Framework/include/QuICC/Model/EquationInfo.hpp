/** 
 * @file EquationInfo.hpp
 * @brief Information describing equation setup from backend
 */

#ifndef QUICC_MODEL_EQUATIONINFO_HPP
#define QUICC_MODEL_EQUATIONINFO_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Information describing equation setup from backend
    */
   struct EquationInfo
   {
         /// List of field ids
         typedef std::vector<SpectralFieldId> SpectralFieldIds;

         /**
          * @brief Constructor
          */
         EquationInfo()
            : isComplex(false), isSplitEquation(false), indexMode(-1)
         {};

         /**
          * @brief Destructor
          */
         ~EquationInfo() = default;

         /**
          * @brief Equation is complex?
          */
         bool isComplex;

         /**
          * @brief Equation is split into lower order systems?
          */
         bool isSplitEquation;

         /**
          * @brief Index mode
          */
         int indexMode;

         /**
          * @brief Implicit coupled fields
          */
         SpectralFieldIds im;

         /**
          * @brief Explicit linear fields
          */
         SpectralFieldIds exL;

         /**
          * @brief Explicit nonlinear fields
          */
         SpectralFieldIds exNL;

         /**
          * @brief Explicit next step fields
          */
         SpectralFieldIds exNS;
   };

}
}

#endif // QUICC_MODEL_EQUATIONINFO_HPP
