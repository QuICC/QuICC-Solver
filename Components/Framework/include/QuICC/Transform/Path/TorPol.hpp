/**
 * @file TorPol.hpp
 * @brief Path for Toroidal/Poloidal from physical values
 */

#ifndef QUICC_TRANSFORM_PATH_TORPOL_HPP
#define QUICC_TRANSFORM_PATH_TORPOL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Path/IId.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   /**
    * @brief Path for Toroidal/Poloidal from physical values
    */
   class TorPol: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         TorPol();

         /**
          * @brief Destructor
          */
         virtual ~TorPol();

         /**
          * @brief Unique id
          */
         static const std::size_t& id();
      
      protected:

      private:
         /**
          * @brief Unique tag
          */
         static std::string sTag();
   };

}
}
}

#endif // QUICC_TRANSFORM_PATH_TORPOL_HPP
