/**
 * @file I2CurlCurlNL.hpp
 * @brief Path for second integral of curlcurl of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_I2CURLCURLNL_HPP
#define QUICC_TRANSFORM_PATH_I2CURLCURLNL_HPP

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
    * @brief Path for second integral of curlcurl of nonlinear terms 
    */
   class I2CurlCurlNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         I2CurlCurlNL();

         /**
          * @brief Destructor
          */
         virtual ~I2CurlCurlNL();

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

#endif // QUICC_TRANSFORM_PATH_I2CURLCURLNL_HPP
