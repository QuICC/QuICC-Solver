/**
 * @file Theta.hpp
 * @brief Theta number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_THETA_HPP
#define QUICC_NONDIMENSIONAL_THETA_HPP

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
    * @brief Theta number nondimensional number
    */
   class Theta: public IRegisterId<Theta>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Theta number
          */
         Theta(const MHDFloat value);

         friend class IRegisterId<Theta>;

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

#endif // QUICC_NONDIMENSIONAL_THETA_HPP
