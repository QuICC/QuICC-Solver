/**
 * @file CflAlfvenScale.hpp
 * @brief Alfven scale for CFL number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_CFLALFVENSCALE_HPP
#define QUICC_NONDIMENSIONAL_CFLALFVENSCALE_HPP

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
    * @brief Alfven scale for CFL number nondimensional number
    */
   class CflAlfvenScale: public IRegisterId<CflAlfvenScale>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Alfven scale for CFL number
          */
         CflAlfvenScale(const MHDFloat value);

         friend class IRegisterId<CflAlfvenScale>;

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

#endif // QUICC_NONDIMENSIONAL_CFLALFVENSCALE_HPP
