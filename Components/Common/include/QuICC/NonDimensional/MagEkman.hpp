/**
 * @file MagEkman.hpp
 * @brief magnetic Ekman number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_MAGEKMAN_HPP
#define QUICC_NONDIMENSIONAL_MAGEKMAN_HPP

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
    * @brief magnetic Ekman number nondimensional number
    */
   class MagEkman: public IRegisterId<MagEkman>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of magnetic Ekman number
          */
         MagEkman(const MHDFloat value);

         friend class IRegisterId<MagEkman>;

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

#endif // QUICC_NONDIMENSIONAL_MAGEKMAN_HPP
