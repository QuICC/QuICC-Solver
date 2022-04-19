/**
 * @file Rayleigh.hpp
 * @brief Rayleigh number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_RAYLEIGH_HPP
#define QUICC_NONDIMENSIONAL_RAYLEIGH_HPP

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
    * @brief Rayleigh number nondimensional number
    */
   class Rayleigh: public IRegisterId<Rayleigh>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Rayleigh number
          */
         Rayleigh(const MHDFloat value);

         friend class IRegisterId<Rayleigh>;

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

#endif // QUICC_NONDIMENSIONAL_RAYLEIGH_HPP
