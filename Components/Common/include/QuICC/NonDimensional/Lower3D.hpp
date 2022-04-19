/**
 * @file Lower3D.hpp
 * @brief Lower3D number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_LOWER3D_HPP
#define QUICC_NONDIMENSIONAL_LOWER3D_HPP

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
    * @brief Lower3D number nondimensional number
    */
   class Lower3D: public IRegisterId<Lower3D>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Lower3D number
          */
         Lower3D(const MHDFloat value);

         friend class IRegisterId<Lower3D>;

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

#endif // QUICC_NONDIMENSIONAL_LOWER3D_HPP
