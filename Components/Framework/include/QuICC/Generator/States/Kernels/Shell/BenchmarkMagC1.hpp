/**
 * @file BenchmarkMagC1.hpp
 * @brief Magnetic benchmark state from Christensen's C1 test case generator kernel
 *        DOI: https://doi.org/10.1016/S0031-9201(01)00275-8
 */

#ifndef QUICC_PHYSICAL_KERNEL_SHELL_BENCHMARKMAGC1_HPP
#define QUICC_PHYSICAL_KERNEL_SHELL_BENCHMARKMAGC1_HPP

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
    * @brief Magnetic benchmark state C1 generator kernel
    */
   class BenchmarkMagC1: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit BenchmarkMagC1();

         /**
          * @brief Simple empty destructor
          */
         virtual ~BenchmarkMagC1();

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

#endif // QUICC_PHYSICAL_KERNEL_SHELL_BENCHMARKMAGC1_HPP
