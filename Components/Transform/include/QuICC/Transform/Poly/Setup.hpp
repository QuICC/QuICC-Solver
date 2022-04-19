/** 
 * @file Setup.hpp
 * @brief Implementation of the polynomial transform setup class
 */

#ifndef QUICC_TRANSFORM_POLY_SETUP_HPP
#define QUICC_TRANSFORM_POLY_SETUP_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

   /**
    * @brief Implementation of the polynomial transform setup class
    */ 
   class Setup: public TransformSetup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param specSize   Spectral output size (i.e without the padding)
          */
         Setup(const int size, const int specSize, const GridPurpose::Id purpose);

         /**
          * @brief Empty destructor
          */
         virtual ~Setup();

      protected:

      private:
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

}
}
}

#endif // QUICC_TRANSFORM_POLY_SETUP_HPP
