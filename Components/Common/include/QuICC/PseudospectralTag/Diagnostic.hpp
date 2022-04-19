/**
 * @file Diagnostic.hpp
 * @brief Diagnostic PseudospectralTag 
 */

#ifndef QUICC_PSEUDOSPECTRALTAG_DIAGNOSTIC_HPP
#define QUICC_PSEUDOSPECTRALTAG_DIAGNOSTIC_HPP

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
    * @brief Diagnostic PseudospectralTag
    */
   class Diagnostic: public IRegisterId<Diagnostic>
   {
      public:
         /**
          * @brief Constructor
          */
         Diagnostic();

         friend class IRegisterId<Diagnostic>;
      
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

#endif // QUICC_PSEUDOSPECTRALTAG_DIAGNOSTIC_HPP
