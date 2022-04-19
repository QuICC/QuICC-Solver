/**
 * @file Lower2D.hpp
 * @brief Lower2D number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_LOWER2D_HPP
#define QUICC_NONDIMENSIONAL_LOWER2D_HPP

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
    * @brief Lower2D number nondimensional number
    */
   class Lower2D: public IRegisterId<Lower2D>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Lower2D number
          */
         Lower2D(const MHDFloat value);

         friend class IRegisterId<Lower2D>;

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

#endif // QUICC_NONDIMENSIONAL_LOWER2D_HPP
