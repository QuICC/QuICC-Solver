/**
 * @file Rescaled.hpp
 * @brief Rescaled number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_RESCALED_HPP
#define QUICC_NONDIMENSIONAL_RESCALED_HPP

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
    * @brief Rescaled number nondimensional number
    */
   class Rescaled: public IRegisterId<Rescaled>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Rescaled number
          */
         Rescaled(const MHDFloat value);

         friend class IRegisterId<Rescaled>;

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

#endif // QUICC_NONDIMENSIONAL_RESCALED_HPP
