/**
 * @file Prandtl.hpp
 * @brief Prandtl number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_PRANDTL_HPP
#define QUICC_NONDIMENSIONAL_PRANDTL_HPP

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
    * @brief Prandtl number nondimensional number
    */
   class Prandtl: public IRegisterId<Prandtl>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Prandtl number
          */
         Prandtl(const MHDFloat value);

         friend class IRegisterId<Prandtl>;

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

#endif // QUICC_NONDIMENSIONAL_PRANDTL_HPP
