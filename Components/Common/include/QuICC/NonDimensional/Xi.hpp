/**
 * @file Xi.hpp
 * @brief Xi number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_XI_HPP
#define QUICC_NONDIMENSIONAL_XI_HPP

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
    * @brief Xi number nondimensional number
    */
   class Xi: public IRegisterId<Xi>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Xi number
          */
         Xi(const MHDFloat value);

         friend class IRegisterId<Xi>;

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

#endif // QUICC_NONDIMENSIONAL_XI_HPP
