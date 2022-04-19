/**
 * @file Pi.hpp
 * @brief Pi number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_PI_HPP
#define QUICC_NONDIMENSIONAL_PI_HPP

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
    * @brief Pi number nondimensional number
    */
   class Pi: public IRegisterId<Pi>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Pi number
          */
         Pi(const MHDFloat value);

         friend class IRegisterId<Pi>;

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

#endif // QUICC_NONDIMENSIONAL_PI_HPP
