/**
 * @file Prognostic.hpp
 * @brief Prognostic SolveTiming
 */

#ifndef QUICC_SOLVETIMING_PROGNOSTIC_HPP
#define QUICC_SOLVETIMING_PROGNOSTIC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/SolveTiming/IRegisterId.hpp"

namespace QuICC {

namespace SolveTiming {

   /**
    * @brief Prognostic SolveTiming
    */
   class Prognostic: public IRegisterId<Prognostic>
   {
      public:
         /**
          * @brief Constructor
          */
         Prognostic();

         friend class IRegisterId<Prognostic>;

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

#endif // QUICC_SOLVETIMING_PROGNOSTIC_HPP
