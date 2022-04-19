/**
 * @file Eady.hpp
 * @brief Eady number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_EADY_HPP
#define QUICC_NONDIMENSIONAL_EADY_HPP

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
    * @brief Eady number nondimensional number
    */
   class Eady: public IRegisterId<Eady>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Eady number
          */
         Eady(const MHDFloat value);

         friend class IRegisterId<Eady>;

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

#endif // QUICC_NONDIMENSIONAL_EADY_HPP
