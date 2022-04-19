/**
 * @file ValidationTorPol.hpp
 * @brief Test case for validation Toroidal/Poloidal projection
 */

#ifndef QUICC_PHYSICAL_KERNEL_SPHERE_VALIDATIONTORPOL_HPP
#define QUICC_PHYSICAL_KERNEL_SPHERE_VALIDATIONTORPOL_HPP

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

namespace Sphere {

   /**
    * @brief Test case for validation Toroidal/Poloidal projection
    */
   class ValidationTorPol: public IPhysicalKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit ValidationTorPol();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ValidationTorPol();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init();

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
          * @brief Cartesian X in spherical coordinates
          */
         Array x(const int iR, const int iTh) const;

         /**
          * @brief Cartesian Y in spherical coordinates
          */
         Array y(const int iR, const int iTh) const;

         /**
          * @brief Cartesian Z in spherical coordinates
          */
         Array z(const int iR, const int iTh) const;

         /**
          * @brief X component of cartesian vector v
          */
         Array vx(const int iR, const int iTh) const;

         /**
          * @brief Y component of cartesian vector v
          */
         Array vy(const int iR, const int iTh) const;

         /**
          * @brief Z component of cartesian vector v
          */
         Array vz(const int iR, const int iTh) const;
   };

}
}
}
}

#endif // QUICC_PHYSICAL_KERNEL_SPHERE_VALIDATIONTORPOL_HPP
