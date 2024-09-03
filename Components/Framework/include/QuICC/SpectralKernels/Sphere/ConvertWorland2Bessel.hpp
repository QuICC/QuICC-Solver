/**
 * @file ConvertWorland2Bessel.hpp
 * @brief Convert Worland expansion to spherical Bessel basis
 */

#ifndef QUICC_SPECTRAL_KERNEL_SPHERE_CONVERTWORLAND2BESSEL_HPP
#define QUICC_SPECTRAL_KERNEL_SPHERE_CONVERTWORLAND2BESSEL_HPP

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
 * @brief Convert Worland expansion to spherical Bessel basis
 */
class ConvertWorland2Bessel : public ISpectralKernel
{
public:
   /// Enum for Worland kind
   enum class WorlandKind
   {
      Chebyshev = 0,
      Legendre,
      SphEnergy,
      CylEnergy,
   };

   /// Enum for Bessel kind
   enum class BesselKind
   {
      Value = 0,
      Insulating,
      NoSlip,
   };

   /**
    * @brief Simple constructor
    */
   explicit ConvertWorland2Bessel(const bool isComplex);

   /**
    * @brief Simple empty destructor
    */
   virtual ~ConvertWorland2Bessel() = default;

   /**
    * @brief Initialize kernel
    *
    * @param comp       Field component to act on
    * @param inWType    Input Worland type
    * @param outBType   Output Bessel type
    * @param scale      Scaling factor
    */
   void init(const FieldComponents::Spectral::Id comp, const WorlandKind inWType, const BesselKind outBType, const MHDFloat scale = 1.0);

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
    * @brief Input Worland type
    */
   WorlandKind mInWType;

   /**
    * @brief Output Bessel type
    */
   BesselKind mOutBType;

   /**
    * @brief Scaling factor
    */
   MHDFloat mScale;
};

/// Typedef for a smart ConvertWorland2Bessel
typedef std::shared_ptr<ConvertWorland2Bessel> SharedConvertWorland2Bessel;

} // namespace Sphere
} // namespace Kernel
} // namespace Spectral
} // namespace QuICC

#endif // QUICC_SPECTRAL_KERNEL_SPHERE_CONVERTWORLAND2BESSEL_HPP
