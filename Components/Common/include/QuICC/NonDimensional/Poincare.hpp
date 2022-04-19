/**
 * @file Poincare.hpp
 * @brief Poincare number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_POINCARE_HPP
#define QUICC_NONDIMENSIONAL_POINCARE_HPP

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
    * @brief Poincare number nondimensional number
    */
   class Poincare: public IRegisterId<Poincare>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Poincare number
          */
         Poincare(const MHDFloat value);

         friend class IRegisterId<Poincare>;

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

#endif // QUICC_NONDIMENSIONAL_POINCARE_HPP
