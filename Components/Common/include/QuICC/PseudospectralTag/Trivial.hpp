/**
 * @file Trivial.hpp
 * @brief Trivial PseudospectralTag 
 */

#ifndef QUICC_PSEUDOSPECTRALTAG_TRIVIAL_HPP
#define QUICC_PSEUDOSPECTRALTAG_TRIVIAL_HPP

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
    * @brief Trivial PseudospectralTag
    */
   class Trivial: public IRegisterId<Trivial>
   {
      public:
         /**
          * @brief Constructor
          */
         Trivial();

         friend class IRegisterId<Trivial>;
      
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

#endif // QUICC_PSEUDOSPECTRALTAG_TRIVIAL_HPP
