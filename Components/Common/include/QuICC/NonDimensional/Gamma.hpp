/**
 * @file Gamma.hpp
 * @brief Gamma number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_GAMMA_HPP
#define QUICC_NONDIMENSIONAL_GAMMA_HPP

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
    * @brief Gamma number nondimensional number
    */
   class Gamma: public IRegisterId<Gamma>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Gamma number
          */
         Gamma(const MHDFloat value);

         friend class IRegisterId<Gamma>;

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

#endif // QUICC_NONDIMENSIONAL_GAMMA_HPP
