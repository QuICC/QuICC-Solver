/**
 * @file Chi.hpp
 * @brief Chi number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_CHI_HPP
#define QUICC_NONDIMENSIONAL_CHI_HPP

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
    * @brief Chi number nondimensional number
    */
   class Chi: public IRegisterId<Chi>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Chi number
          */
         Chi(const MHDFloat value);

         friend class IRegisterId<Chi>;

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

#endif // QUICC_NONDIMENSIONAL_CHI_HPP
