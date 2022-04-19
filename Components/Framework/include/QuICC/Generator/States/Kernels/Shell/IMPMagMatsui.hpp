/**
 * @file IMPMagMatsui.hpp
 * @brief Magnetic IMP Matsui state generator kernel
 */

#ifndef QUICC_PHYSICAL_KERNEL_SHELL_IMPMAGMATSUI_HPP
#define QUICC_PHYSICAL_KERNEL_SHELL_IMPMAGMATSUI_HPP

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

namespace Shell {

   /**
    * @brief Magnetic IMP Matsui C1 generator kernel
    */
   class IMPMagMatsui: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit IMPMagMatsui();

         /**
          * @brief Simple empty destructor
          */
         virtual ~IMPMagMatsui();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDFloat ri, const MHDFloat ro);

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
          * @brief Inner radius of spherical shell
          */
         MHDFloat mRi;

         /**
          * @brief Outer radius of spherical shell
          */
         MHDFloat mRo;
   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SHELL_IMPMAGMATSUI_HPP
