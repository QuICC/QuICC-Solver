/**
 * @file Iota.hpp
 * @brief Iota number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_IOTA_HPP
#define QUICC_NONDIMENSIONAL_IOTA_HPP

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
    * @brief Iota number nondimensional number
    */
   class Iota: public IRegisterId<Iota>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Iota number
          */
         Iota(const MHDFloat value);

         friend class IRegisterId<Iota>;

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

#endif // QUICC_NONDIMENSIONAL_IOTA_HPP
