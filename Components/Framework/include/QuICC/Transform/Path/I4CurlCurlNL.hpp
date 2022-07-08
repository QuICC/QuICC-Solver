/**
 * @file I4CurlCurlNL.hpp
 * @brief Path for fourth integral of curlcurl of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_I4CURLCURLNL_HPP
#define QUICC_TRANSFORM_PATH_I4CURLCURLNL_HPP

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
    * @brief Path for fourth integral of curlcurl of nonlinear terms 
    */
   class I4CurlCurlNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         I4CurlCurlNL();

         /**
          * @brief Destructor
          */
         virtual ~I4CurlCurlNL();

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

#endif // QUICC_TRANSFORM_PATH_I4CURLCURLNL_HPP
