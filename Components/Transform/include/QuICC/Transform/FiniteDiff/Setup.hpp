/**
 * @file Setup.hpp
 * @brief Implementation of the finite differences transform setup class
 */

#ifndef QUICC_TRANSFORM_FINITEDIFF_SETUP_HPP
#define QUICC_TRANSFORM_FINITEDIFF_SETUP_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace Transform {

namespace FiniteDiff {

   /**
    * @brief Implementation of the finite differences transform setup class
    */
   class Setup: public TransformSetup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          */
         Setup(const int size, const GridPurpose::Id purpose);

         /**
          * @brief Empty destructor
          */
         virtual ~Setup() = default;

      protected:

      private:
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC

#endif // QUICC_TRANSFORM_FINITEDIFF_SETUP_HPP
