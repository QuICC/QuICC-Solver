/**
 * @file CflAlfvenDamping.hpp
 * @brief Alfven damping for CFL number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_CFLALFVENDAMPING_HPP
#define QUICC_NONDIMENSIONAL_CFLALFVENDAMPING_HPP

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
    * @brief Alfven damping for CFL number nondimensional number
    */
   class CflAlfvenDamping: public IRegisterId<CflAlfvenDamping>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Alfven damping for CFL number
          */
         CflAlfvenDamping(const MHDFloat value);

         friend class IRegisterId<CflAlfvenDamping>;

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

#endif // QUICC_NONDIMENSIONAL_CFLALFVENDAMPING_HPP
