/**
 * @file Sigma.hpp
 * @brief Sigma number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_SIGMA_HPP
#define QUICC_NONDIMENSIONAL_SIGMA_HPP

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
    * @brief Sigma number nondimensional number
    */
   class Sigma: public IRegisterId<Sigma>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Sigma number
          */
         Sigma(const MHDFloat value);

         friend class IRegisterId<Sigma>;

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

#endif // QUICC_NONDIMENSIONAL_SIGMA_HPP
