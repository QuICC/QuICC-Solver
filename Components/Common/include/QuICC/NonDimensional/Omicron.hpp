/**
 * @file Omicron.hpp
 * @brief Omicron number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_OMICRON_HPP
#define QUICC_NONDIMENSIONAL_OMICRON_HPP

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
    * @brief Omicron number nondimensional number
    */
   class Omicron: public IRegisterId<Omicron>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Omicron number
          */
         Omicron(const MHDFloat value);

         friend class IRegisterId<Omicron>;

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

#endif // QUICC_NONDIMENSIONAL_OMICRON_HPP
