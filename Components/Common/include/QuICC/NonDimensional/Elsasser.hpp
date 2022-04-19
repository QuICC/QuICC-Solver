/**
 * @file Elsasser.hpp
 * @brief Elsasser number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ELSASSER_HPP
#define QUICC_NONDIMENSIONAL_ELSASSER_HPP

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
    * @brief Elsasser number nondimensional number
    */
   class Elsasser: public IRegisterId<Elsasser>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Elsasser number
          */
         Elsasser(const MHDFloat value);

         friend class IRegisterId<Elsasser>;

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

#endif // QUICC_NONDIMENSIONAL_ELSASSER_HPP
