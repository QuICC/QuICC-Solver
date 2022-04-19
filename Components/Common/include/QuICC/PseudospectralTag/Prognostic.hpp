/**
 * @file Prognostic.hpp
 * @brief Prognostic PseudospectralTag 
 */

#ifndef QUICC_PSEUDOSPECTRALTAG_PROGNOSTIC_HPP
#define QUICC_PSEUDOSPECTRALTAG_PROGNOSTIC_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PseudospectralTag/IRegisterId.hpp"

namespace QuICC {

namespace PseudospectralTag {

   /**
    * @brief Prognostic PseudospectralTag
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

#endif // QUICC_PSEUDOSPECTRALTAG_PROGNOSTIC_HPP
