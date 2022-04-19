/**
 * @file Upper1D.hpp
 * @brief Upper1D number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_UPPER1D_HPP
#define QUICC_NONDIMENSIONAL_UPPER1D_HPP

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
    * @brief Upper1D number nondimensional number
    */
   class Upper1D: public IRegisterId<Upper1D>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Upper1D number
          */
         Upper1D(const MHDFloat value);

         friend class IRegisterId<Upper1D>;

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

#endif // QUICC_NONDIMENSIONAL_UPPER1D_HPP
