/**
 * @file Epsilon.hpp
 * @brief Epsilon number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_EPSILON_HPP
#define QUICC_NONDIMENSIONAL_EPSILON_HPP

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
    * @brief Epsilon number nondimensional number
    */
   class Epsilon: public IRegisterId<Epsilon>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Epsilon number
          */
         Epsilon(const MHDFloat value);

         friend class IRegisterId<Epsilon>;

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

#endif // QUICC_NONDIMENSIONAL_EPSILON_HPP
