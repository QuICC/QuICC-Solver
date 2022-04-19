/** 
 * @file Setup.hpp
 * @brief Implementation of the polynomial transform setup class
 */

#ifndef QUICC_TRANSFORM_POLY_ALEGENDRE_SETUP_HPP
#define QUICC_TRANSFORM_POLY_ALEGENDRE_SETUP_HPP

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Poly/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace ALegendre {

   /**
    * @brief Specialization of the polynomial transform setup class for associated Legendre transforms
    */ 
   class Setup: public ::QuICC::Transform::Poly::Setup
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

         /**
          * @brief Add index with multiplicity and assume full fast index list
          */
         virtual void addIndex(const int slowIdx, const int mult = 1);

         using ::QuICC::Transform::Poly::Setup::addIndex;

      protected:

      private:
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

}
}
}
}

#endif // QUICC_TRANSFORM_POLY_ALEGENDRE_SETUP_HPP
