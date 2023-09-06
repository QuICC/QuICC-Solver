/**
 * @file ConserveAngularMomentum.hpp
 * @brief Angular momentum conservation kernel to in sphere
 */

#ifndef QUICC_SPECTRAL_KERNEL_SPHERE_CONSERVEANGULARMOMENTUM_HPP
#define QUICC_SPECTRAL_KERNEL_SPHERE_CONSERVEANGULARMOMENTUM_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

namespace Sphere {

   /**
    * @brief Angular momentum conservation kernel in sphere
    */
   class ConserveAngularMomentum: public ISpectralKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit ConserveAngularMomentum(const bool isComplex);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ConserveAngularMomentum() = default;

         /**
          * @brief Initialize kernel
          */
         void init(const bool hasMOrdering);

         /**
          * @brief Compute the spectral kernel
          *
          * @param id   Component ID
          * @param i    Fastest index
          * @param j    Second index
          * @param k    Slowest index
          */
         virtual MHDVariant compute(const int i, const int j, const int k) const final;

         /**
          * @brief Apply kernel to field
          */
         virtual void apply(const std::size_t timeId) final;
         
      protected:

      private:
         /**
          * @brief m = 0 mode is present
          */
         bool mHasM0;

         /**
          * @brief m = 1 mode is present
          */
         bool mHasM1;

         /**
          * @brief j index of m = 0 mode
          */
         int mM0j;

         /**
          * @brief j index of m = 0 mode
          */
         int mM0k;

         /**
          * @brief j index of m = 1 mode
          */
         int mM1j;

         /**
          * @brief j index of m = 1 mode
          */
         int mM1k;

         /**
          * @brief Angular momentum conservation operator
          */
         Matrix mOp;

   };

   /// Typedef for a smart ConserveAngularMomentum
   typedef std::shared_ptr<ConserveAngularMomentum> SharedConserveAngularMomentum;
   
} // Sphere
} // Kernel
} // Spectral
} // QuICC

#endif // QUICC_SPECTRAL_KERNEL_SPHERE_CONSERVEANGULARMOMENTUM_HPP
