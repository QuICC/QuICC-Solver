/**
 * @file VariableBase.hpp
 * @brief Base of the implementation of the variables
 */

#ifndef QUICC_DATATYPES_VARIABLEBASE_HPP
#define QUICC_DATATYPES_VARIABLEBASE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Base of the implementation of the variables
    */
   class VariableBase
   {
      public:
         /**
         * @brief Construct the general shared information for a physical variable
         *
         * @param spRes Resolution information
         */
         VariableBase(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~VariableBase();

         /**
          * @brief Get resolution information
          */
         const SharedResolution  spRes() const;

         /**
          * @brief Get resolution information
          */
         const Resolution&  res() const;

         /**
         * @brief Get the memory requirements
         */
         virtual MHDFloat requiredStorage() const;

      protected:

      private:
         /**
          * @brief Pointer to resolution information
          */
         SharedResolution   mspRes;
   };

   inline const SharedResolution VariableBase::spRes() const
   {
      return this->mspRes;
   }

   inline const Resolution& VariableBase::res() const
   {
      return *this->mspRes;
   }

}
}

#endif // QUICC_DATATYPES_VARIABLEBASE_HPP
