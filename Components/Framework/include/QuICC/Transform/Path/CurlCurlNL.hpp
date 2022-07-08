/**
 * @file CurlCurlNL.hpp
 * @brief Path for curlcurl of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_CURLCURLNL_HPP
#define QUICC_TRANSFORM_PATH_CURLCURLNL_HPP

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
    * @brief Path for curlcurl of nonlinear terms 
    */
   class CurlCurlNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         CurlCurlNL();

         /**
          * @brief Destructor
          */
         virtual ~CurlCurlNL();

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

#endif // QUICC_TRANSFORM_PATH_CURLCURLNL_HPP
