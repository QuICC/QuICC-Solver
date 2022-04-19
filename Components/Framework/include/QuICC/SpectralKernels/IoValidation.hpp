/**
 * @file IoValidation.hpp
 * @brief Trivial kernel to make field constant
 */

#ifndef QUICC_SPECTRAL_KERNEL_IOVALIDATION_HPP
#define QUICC_SPECTRAL_KERNEL_IOVALIDATION_HPP

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

namespace QuICC {

namespace Spectral {

namespace Kernel {

   /**
    * @brief Trivial passthrough kernel
    */
   class IoValidation: public ISpectralKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit IoValidation(const bool isComplex);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IoValidation();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDFloat range);

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
          * @brief Range for indexe
          */
         MHDFloat mRange;

   };

   /// Typedef for a smart IoValidation
   typedef std::shared_ptr<IoValidation> SharedIoValidation;
   
}
}
}

#endif // QUICC_SPECTRAL_KERNEL_IOVALIDATION_HPP
