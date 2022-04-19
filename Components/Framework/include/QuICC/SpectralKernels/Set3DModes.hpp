/**
 * @file Set3DModes.hpp
 * @brief Set given spectral modes
 */

#ifndef QUICC_SPECTRAL_KERNEL_SET3DMODES_HPP
#define QUICC_SPECTRAL_KERNEL_SET3DMODES_HPP

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
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   /**
    * @brief Set given spectral modes
    */
   class Set3DModes: public ISpectralKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit Set3DModes(const bool isComplex);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Set3DModes();

         /**
          * @brief Initialize kernel
          */
         void init(const Real3DMapType modes);

         /**
          * @brief Initialize kernel
          */
         void init(const Complex3DMapType modes);

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDComplex modes);

         /**
          * @brief Compute the spectral kernel
          *
          * @param id   Component ID
          * @param i    Fastest index
          * @param j    Second index
          * @param k    Slowest index
          */
         virtual MHDVariant compute(const int i, const int j, const int k) const;
         
      protected:

      private:
         /**
          * @brief Real value to set the field to
          */
         SharedReal3DMapType mspRModes;

         /**
          * @brief Complex value to set the field to
          */
         SharedComplex3DMapType mspZModes;

   };

   /// Typedef for a smart Set3DModes
   typedef std::shared_ptr<Set3DModes> SharedSet3DModes;
   
}
}
}

#endif // QUICC_SPECTRAL_KERNEL_SET3DMODES_HPP
