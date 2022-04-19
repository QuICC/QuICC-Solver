/**
 * @file MakeRandom.hpp
 * @brief Trivial kernel to make field constant
 */

#ifndef QUICC_PHYSICAL_KERNEL_MAKERANDOM_HPP
#define QUICC_PHYSICAL_KERNEL_MAKERANDOM_HPP

// First include
//

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   /**
    * @brief Trivial passthrough kernel
    */
   class MakeRandom: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit MakeRandom();

         /**
          * @brief Simple empty destructor
          */
         virtual ~MakeRandom();

         /**
          * @brief Initialize kernel
          *
          * @param value Scale of random values
          */
         void init(const MHDFloat scale);

         /**
          * @brief Compute the physical kernel
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const;
         
      protected:

      private:
         /**
          * @brief Scale random value
          */
         MHDFloat mScale;

   };

   /// Typedef for a smart MakeRandom
   typedef std::shared_ptr<MakeRandom> SharedMakeRandom;
   
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_MAKERANDOM_HPP
