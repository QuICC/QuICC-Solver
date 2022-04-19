/**
 * @file Ekman.hpp
 * @brief Ekman number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_EKMAN_HPP
#define QUICC_NONDIMENSIONAL_EKMAN_HPP

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
    * @brief Ekman number nondimensional number
    */
   class Ekman: public IRegisterId<Ekman>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Ekman number
          */
         Ekman(const MHDFloat value);

         friend class IRegisterId<Ekman>;

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

#endif // QUICC_NONDIMENSIONAL_EKMAN_HPP
