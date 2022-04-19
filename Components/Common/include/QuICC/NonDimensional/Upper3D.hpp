/**
 * @file Upper3D.hpp
 * @brief Upper3D number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_UPPER3D_HPP
#define QUICC_NONDIMENSIONAL_UPPER3D_HPP

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
    * @brief Upper3D number nondimensional number
    */
   class Upper3D: public IRegisterId<Upper3D>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Upper3D number
          */
         Upper3D(const MHDFloat value);

         friend class IRegisterId<Upper3D>;

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

#endif // QUICC_NONDIMENSIONAL_UPPER3D_HPP
