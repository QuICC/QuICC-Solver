/**
 * @file Roberts.hpp
 * @brief Roberts number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ROBERTS_HPP
#define QUICC_NONDIMENSIONAL_ROBERTS_HPP

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
    * @brief Roberts number nondimensional number
    */
   class Roberts: public IRegisterId<Roberts>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Roberts number
          */
         Roberts(const MHDFloat value);

         friend class IRegisterId<Roberts>;

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

#endif // QUICC_NONDIMENSIONAL_ROBERTS_HPP
