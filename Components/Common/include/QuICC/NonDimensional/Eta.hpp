/**
 * @file Eta.hpp
 * @brief Eta number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_ETA_HPP
#define QUICC_NONDIMENSIONAL_ETA_HPP

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
    * @brief Eta number nondimensional number
    */
   class Eta: public IRegisterId<Eta>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Eta number
          */
         Eta(const MHDFloat value);

         friend class IRegisterId<Eta>;

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

#endif // QUICC_NONDIMENSIONAL_ETA_HPP
