/**
 * @file MakeRandom.hpp
 * @brief Trivial kernel to make field constant
 */

#ifndef QUICC_SPECTRAL_KERNEL_MAKERANDOM_HPP
#define QUICC_SPECTRAL_KERNEL_MAKERANDOM_HPP

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
   class MakeRandom: public ISpectralKernel
   {
      public:
         /**
          * @brief Simple constructor
          */
         explicit MakeRandom(const bool isComplex);

         /**
          * @brief Simple empty destructor
          */
         virtual ~MakeRandom();

         /**
          * @brief Initialize kernel
          *
          * @param value Scale of random values
          */
         void init(const MHDFloat minVal, const MHDFloat maxVal);

         /**
          * @brief Set the convergence ratio of spectrum
          */
         void setRatio(const std::vector<MHDFloat>& ratios);

         /**
          * @brief Set number of zeros at end of spectrum
          */
         void setZeros(const std::vector<int>& zeros);

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
         /**
          * @brief Get random value
          */
         void value(MHDFloat& val, const unsigned int seed = 1) const;

         /**
          * @brief Get random value
          */
         void value(MHDComplex& val, const unsigned int seed = 1) const;

         /**
          * @brief Get scaling ratio
          */
         virtual MHDFloat scalingRatio(const int i, const int n, const int dim) const;

         /**
          * @brief Compute mode value
          */
         template <typename T> void mode(T& val, const int i, const int j, const int k) const;

         /**
          * @brief Impose constraint on modes (ie. purely real entries, etc)
          */
         virtual void constrain(MHDFloat& val, const int i, const int j, const int k) const;

         /**
          * @brief Impose constraint modes (ie. purely real entries, etc)
          */
         virtual void constrain(MHDComplex& val, const int i, const int j, const int k) const;

      private:
         /**
          * @brief Starting seed used for random generator
          */
         unsigned int mStartSeed;

         /**
          * @brief Minimum value
          */
         MHDFloat mMin;

         /**
          * @brief Maximum value
          */
         MHDFloat mMax;

         /**
          * @brief Imposed convergence ratio on spectrum
          */
         std::vector<MHDFloat> mRatio;

         /**
          * @brief Number of zero modes at end of spectrum
          */
         std::vector<int> mZeros;

   };

   /// Typedef for a smart MakeRandom
   typedef std::shared_ptr<MakeRandom> SharedMakeRandom;
   
}
}
}

#endif // QUICC_SPECTRAL_KERNEL_MAKERANDOM_HPP
