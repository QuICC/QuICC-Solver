/**
 * @file Tau.hpp
 * @brief Tau number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_TAU_HPP
#define QUICC_NONDIMENSIONAL_TAU_HPP

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
    * @brief Tau number nondimensional number
    */
   class Tau: public IRegisterId<Tau>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Tau number
          */
         Tau(const MHDFloat value);

         friend class IRegisterId<Tau>;

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

#endif // QUICC_NONDIMENSIONAL_TAU_HPP
