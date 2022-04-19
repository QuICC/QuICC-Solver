/**
 * @file Mu.hpp
 * @brief Mu number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_MU_HPP
#define QUICC_NONDIMENSIONAL_MU_HPP

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
    * @brief Mu number nondimensional number
    */
   class Mu: public IRegisterId<Mu>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Mu number
          */
         Mu(const MHDFloat value);

         friend class IRegisterId<Mu>;

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

#endif // QUICC_NONDIMENSIONAL_MU_HPP
