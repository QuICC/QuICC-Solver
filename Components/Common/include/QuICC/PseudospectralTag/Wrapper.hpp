/**
 * @file Wrapper.hpp
 * @brief Wrapper PseudospectralTag 
 */

#ifndef QUICC_PSEUDOSPECTRALTAG_WRAPPER_HPP
#define QUICC_PSEUDOSPECTRALTAG_WRAPPER_HPP

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
    * @brief Wrapper PseudospectralTag
    */
   class Wrapper: public IRegisterId<Wrapper>
   {
      public:
         /**
          * @brief Constructor
          */
         Wrapper();

         friend class IRegisterId<Wrapper>;
      
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

#endif // QUICC_PSEUDOSPECTRALTAG_WRAPPER_HPP
