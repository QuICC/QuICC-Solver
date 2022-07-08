/**
 * @file Scalar.hpp
 * @brief Path for scalar of physical space values
 */

#ifndef QUICC_TRANSFORM_PATH_SCALAR_HPP
#define QUICC_TRANSFORM_PATH_SCALAR_HPP

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
    * @brief Path for scalar of physical values
    */
   class Scalar: public IId
   {
      public:
         /**
          * @brief Constructor
          */
         Scalar();

         /**
          * @brief Destructor
          */
         virtual ~Scalar();

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

#endif // QUICC_TRANSFORM_PATH_SCALAR_HPP
