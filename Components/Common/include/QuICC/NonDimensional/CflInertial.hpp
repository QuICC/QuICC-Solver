/**
 * @file CflInertial.hpp
 * @brief Inertial CFL number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_CFLINERTIAL_HPP
#define QUICC_NONDIMENSIONAL_CFLINERTIAL_HPP

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
    * @brief Inertial CFL number nondimensional number
    */
   class CflInertial: public IRegisterId<CflInertial>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Inertial CFL number
          */
         CflInertial(const MHDFloat value);

         friend class IRegisterId<CflInertial>;

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

#endif // QUICC_NONDIMENSIONAL_CFLINERTIAL_HPP
