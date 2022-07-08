/**
 * @file I2ScalarNL.hpp
 * @brief Path for second integral of scalar of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_I2SCALARNL_HPP
#define QUICC_TRANSFORM_PATH_I2SCALARNL_HPP

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
    * @brief Path for second integral of scalar of nonlinear terms 
    */
   class I2ScalarNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         I2ScalarNL();

         /**
          * @brief Destructor
          */
         virtual ~I2ScalarNL();

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

#endif // QUICC_TRANSFORM_PATH_I2SCALARNL_HPP
