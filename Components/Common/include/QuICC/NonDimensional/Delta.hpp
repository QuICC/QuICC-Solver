/**
 * @file Delta.hpp
 * @brief Delta number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_DELTA_HPP
#define QUICC_NONDIMENSIONAL_DELTA_HPP

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
    * @brief Delta number nondimensional number
    */
   class Delta: public IRegisterId<Delta>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Delta number
          */
         Delta(const MHDFloat value);

         friend class IRegisterId<Delta>;

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

#endif // QUICC_NONDIMENSIONAL_DELTA_HPP
