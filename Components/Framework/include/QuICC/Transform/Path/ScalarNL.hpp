/**
 * @file ScalarNL.hpp
 * @brief Path for scalar of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_SCALARNL_HPP
#define QUICC_TRANSFORM_PATH_SCALARNL_HPP

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
    * @brief Path for scalar of nonlinear terms 
    */
   class ScalarNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         ScalarNL();

         /**
          * @brief Destructor
          */
         virtual ~ScalarNL();

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

#endif // QUICC_TRANSFORM_PATH_SCALARNL_HPP
