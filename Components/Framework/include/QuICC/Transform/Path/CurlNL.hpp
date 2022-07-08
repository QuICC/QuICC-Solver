/**
 * @file CurlNL.hpp
 * @brief Path for curl of nonlinear terms
 */

#ifndef QUICC_TRANSFORM_PATH_CURLNL_HPP
#define QUICC_TRANSFORM_PATH_CURLNL_HPP

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
    * @brief Path for curl of nonlinear terms 
    */
   class CurlNL: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         CurlNL();

         /**
          * @brief Destructor
          */
         virtual ~CurlNL();

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

#endif // QUICC_TRANSFORM_PATH_CURLNL_HPP
