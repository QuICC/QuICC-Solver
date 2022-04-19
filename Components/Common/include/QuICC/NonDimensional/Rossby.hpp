/**
 * @file Rossby.hpp
 * @brief Rossby number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ROSSBY_HPP
#define QUICC_NONDIMENSIONAL_ROSSBY_HPP

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
    * @brief Rossby number nondimensional number
    */
   class Rossby: public IRegisterId<Rossby>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Rossby number
          */
         Rossby(const MHDFloat value);

         friend class IRegisterId<Rossby>;

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

#endif // QUICC_NONDIMENSIONAL_ROSSBY_HPP
