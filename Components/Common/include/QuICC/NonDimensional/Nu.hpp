/**
 * @file Nu.hpp
 * @brief Nu number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_NU_HPP
#define QUICC_NONDIMENSIONAL_NU_HPP

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
    * @brief Nu number nondimensional number
    */
   class Nu: public IRegisterId<Nu>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Nu number
          */
         Nu(const MHDFloat value);

         friend class IRegisterId<Nu>;

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

#endif // QUICC_NONDIMENSIONAL_NU_HPP
