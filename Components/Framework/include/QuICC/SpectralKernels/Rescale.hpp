/**
 * @file Rescale.hpp
 * @brief Rescale field by a constant
 */

#ifndef QUICC_SPECTRAL_KERNEL_RESCALE_HPP
#define QUICC_SPECTRAL_KERNEL_RESCALE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

/**
 * @brief Rescale field by a constant
 */
class Rescale : public ISpectralKernel
{
public:
   /**
    * @brief Simple constructor
    */
   explicit Rescale(const bool isComplex);

   /**
    * @brief Simple empty destructor
    */
   virtual ~Rescale() = default;

   /**
    * @brief Initialize kernel
    *
    * @param comp    Field component to act on
    * @param scale   Scaling constant
    */
   void init(const FieldComponents::Spectral::Id comp, const MHDFloat scale);

   /**
    * @brief Compute the spectral kernel
    *
    * @param id   Component ID
    * @param i    Fastest index
    * @param j    Second index
    * @param k    Slowest index
    */
   virtual MHDVariant compute(const int i, const int j,
      const int k) const final;

   /**
    * @brief Apply kernel to field
    */
   virtual void apply(const std::size_t timeId) final;

protected:
private:
   /**
    * @brief Spectral component to act on
    */
   FieldComponents::Spectral::Id mComp;

   /**
    * @brief Scaling constant
    */
   MHDFloat mScale;
};

} // namespace Kernel
} // namespace Spectral
} // namespace QuICC

#endif // QUICC_SPECTRAL_KERNEL_RESCALE_HPP
