/**
 * @file CflTorsional.hpp
 * @brief Torsional CFL number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_CFLTORSIONAL_HPP
#define QUICC_NONDIMENSIONAL_CFLTORSIONAL_HPP

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
    * @brief Torsional CFL number nondimensional number
    */
   class CflTorsional: public IRegisterId<CflTorsional>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Torsional CFL number
          */
         CflTorsional(const MHDFloat value);

         friend class IRegisterId<CflTorsional>;

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

#endif // QUICC_NONDIMENSIONAL_CFLTORSIONAL_HPP
