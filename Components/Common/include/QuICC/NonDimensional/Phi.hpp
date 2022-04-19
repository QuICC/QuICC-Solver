/**
 * @file Phi.hpp
 * @brief Phi number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_PHI_HPP
#define QUICC_NONDIMENSIONAL_PHI_HPP

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
    * @brief Phi number nondimensional number
    */
   class Phi: public IRegisterId<Phi>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Phi number
          */
         Phi(const MHDFloat value);

         friend class IRegisterId<Phi>;

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

#endif // QUICC_NONDIMENSIONAL_PHI_HPP
