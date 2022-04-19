/**
 * @file Lower1D.hpp
 * @brief Lower1D number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_LOWER1D_HPP
#define QUICC_NONDIMENSIONAL_LOWER1D_HPP

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
    * @brief Lower1D number nondimensional number
    */
   class Lower1D: public IRegisterId<Lower1D>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Lower1D number
          */
         Lower1D(const MHDFloat value);

         friend class IRegisterId<Lower1D>;

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

#endif // QUICC_NONDIMENSIONAL_LOWER1D_HPP
