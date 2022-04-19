/**
 * @file Lambda.hpp
 * @brief Lambda number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_LAMBDA_HPP
#define QUICC_NONDIMENSIONAL_LAMBDA_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/BasicTypes.hpp"
#include "QuICC/NonDimensional/IRegisterId.hpp"

namespace QuICC {

namespace NonDimensional {

   /**
    * @brief Lambda number nondimensional number
    */
   class Lambda: public IRegisterId<Lambda>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Lambda number
          */
         Lambda(const MHDFloat value);

         friend class IRegisterId<Lambda>;

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

#endif // QUICC_NONDIMENSIONAL_LAMBDA_HPP
