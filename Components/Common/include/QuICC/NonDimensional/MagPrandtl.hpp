/**
 * @file MagPrandtl.hpp
 * @brief magnetic Prandtl number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_MAGPRANDTL_HPP
#define QUICC_NONDIMENSIONAL_MAGPRANDTL_HPP

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
    * @brief magnetic Prandtl number nondimensional number
    */
   class MagPrandtl: public IRegisterId<MagPrandtl>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of magnetic Prandtl number
          */
         MagPrandtl(const MHDFloat value);

         friend class IRegisterId<MagPrandtl>;

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

#endif // QUICC_NONDIMENSIONAL_MAGPRANDTL_HPP
