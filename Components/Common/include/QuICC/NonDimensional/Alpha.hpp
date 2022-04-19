/**
 * @file Alpha.hpp
 * @brief Alpha number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ALPHA_HPP
#define QUICC_NONDIMENSIONAL_ALPHA_HPP

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
    * @brief Alpha number nondimensional number
    */
   class Alpha: public IRegisterId<Alpha>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Alpha number
          */
         Alpha(const MHDFloat value);

         friend class IRegisterId<Alpha>;

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

#endif // QUICC_NONDIMENSIONAL_ALPHA_HPP
