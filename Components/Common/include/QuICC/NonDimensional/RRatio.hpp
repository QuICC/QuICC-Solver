/**
 * @file RRatio.hpp
 * @brief R ratio number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_RRATIO_HPP
#define QUICC_NONDIMENSIONAL_RRATIO_HPP

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
    * @brief R ratio number nondimensional number
    */
   class RRatio: public IRegisterId<RRatio>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of R ratio number
          */
         RRatio(const MHDFloat value);

         friend class IRegisterId<RRatio>;

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

#endif // QUICC_NONDIMENSIONAL_RRATIO_HPP
