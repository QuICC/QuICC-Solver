/**
 * @file I2CurlNL.hpp
 * @brief Path for second integral of curl of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_I2CURLNL_HPP
#define QUICC_TRANSFORM_PATH_I2CURLNL_HPP

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
    * @brief Path for second integral of curl of nonlinear terms 
    */
   class I2CurlNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         I2CurlNL();

         /**
          * @brief Destructor
          */
         virtual ~I2CurlNL();

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

#endif // QUICC_TRANSFORM_PATH_I2CURLNL_HPP
