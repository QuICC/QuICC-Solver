/**
 * @file Upper2D.hpp
 * @brief Upper2D number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_UPPER2D_HPP
#define QUICC_NONDIMENSIONAL_UPPER2D_HPP

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
    * @brief Upper2D number nondimensional number
    */
   class Upper2D: public IRegisterId<Upper2D>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Upper2D number
          */
         Upper2D(const MHDFloat value);

         friend class IRegisterId<Upper2D>;

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

#endif // QUICC_NONDIMENSIONAL_UPPER2D_HPP
