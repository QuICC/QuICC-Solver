/**
 * @file Beta.hpp
 * @brief Beta number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_BETA_HPP
#define QUICC_NONDIMENSIONAL_BETA_HPP

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
    * @brief Beta number nondimensional number
    */
   class Beta: public IRegisterId<Beta>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Beta number
          */
         Beta(const MHDFloat value);

         friend class IRegisterId<Beta>;

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

#endif // QUICC_NONDIMENSIONAL_BETA_HPP
