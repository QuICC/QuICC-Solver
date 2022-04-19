/**
 * @file Upsilon.hpp
 * @brief Upsilon number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_UPSILON_HPP
#define QUICC_NONDIMENSIONAL_UPSILON_HPP

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
    * @brief Upsilon number nondimensional number
    */
   class Upsilon: public IRegisterId<Upsilon>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Upsilon number
          */
         Upsilon(const MHDFloat value);

         friend class IRegisterId<Upsilon>;

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

#endif // QUICC_NONDIMENSIONAL_UPSILON_HPP
