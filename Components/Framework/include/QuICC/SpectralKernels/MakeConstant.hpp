/**
 * @file MakeConstant.hpp
 * @brief Trivial kernel to make field constant
 */

#ifndef QUICC_SPECTRAL_KERNEL_MAKECONSTANT_HPP
#define QUICC_SPECTRAL_KERNEL_MAKECONSTANT_HPP

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
   class MakeConstant: public ISpectralKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit MakeConstant(const bool isComplex);

         /**
          * @brief Simple empty destructor
          */
         virtual ~MakeConstant();

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDFloat value);

         /**
          * @brief Initialize kernel
          *
          * @param value Value to set the field to
          */
         void init(const MHDComplex value);

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
         MHDFloat mRValue;

         /**
          * @brief Complex value to set the field to
          */
         MHDComplex mZValue;

   };

   /// Typedef for a smart MakeConstant
   typedef std::shared_ptr<MakeConstant> SharedMakeConstant;
   
}
}
}

#endif // QUICC_SPECTRAL_KERNEL_MAKECONSTANT_HPP
