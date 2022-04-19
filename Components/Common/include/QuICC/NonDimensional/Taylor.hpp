/**
 * @file Taylor.hpp
 * @brief Taylor number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_TAYLOR_HPP
#define QUICC_NONDIMENSIONAL_TAYLOR_HPP

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
    * @brief Taylor number nondimensional number
    */
   class Taylor: public IRegisterId<Taylor>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Taylor number
          */
         Taylor(const MHDFloat value);

         friend class IRegisterId<Taylor>;

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

#endif // QUICC_NONDIMENSIONAL_TAYLOR_HPP
