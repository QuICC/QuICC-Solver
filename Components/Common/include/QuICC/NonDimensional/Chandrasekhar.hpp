/**
 * @file Chandrasekhar.hpp
 * @brief Chandrasekhar number nondimensional number
 */

#ifndef QUICC_NONDIMENSIONAL_CHANDRASEKHAR_HPP
#define QUICC_NONDIMENSIONAL_CHANDRASEKHAR_HPP

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
    * @brief Chandrasekhar number nondimensional number
    */
   class Chandrasekhar: public IRegisterId<Chandrasekhar>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param value Value of Chandrasekhar number
          */
         Chandrasekhar(const MHDFloat value);

         friend class IRegisterId<Chandrasekhar>;

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

#endif // QUICC_NONDIMENSIONAL_CHANDRASEKHAR_HPP
