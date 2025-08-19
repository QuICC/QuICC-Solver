/**
 * @file Truncate.hpp
 * @brief Truncate Worland solution with minimal energy loss
 */

#ifndef QUICC_SPECTRAL_KERNEL_SPHERE_TRUNCATE_HPP
#define QUICC_SPECTRAL_KERNEL_SPHERE_TRUNCATE_HPP

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

namespace Sphere {

/**
 * @brief Truncate Worland solution with minimal energy loss
 */
class Truncate : public ISpectralKernel
{
public:
   /**
    * @brief Simple constructor
    */
   explicit Truncate(const bool isComplex);

   /**
    * @brief Simple empty destructor
    */
   virtual ~Truncate() = default;

   /**
    * @brief Initialize kernel
    *
    * @param comp       Field component to act on
    * @param wType      Input Worland type
    * @param bcId       Boundary condition
    * @param outN       Output radial truncation
    * @param outL       Output spherical degree
    * @param outM       Output spherical order
    */
   void init(const FieldComponents::Spectral::Id comp, const std::string& WType, const std::size_t bcId, const int outN, const int outL, const int outM);

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
    * @brief Output radial truncation
    */
   int mOutN;

   /**
    * @brief Output max spherical degree
    */
   int mOutL;

   /**
    * @brief Output max spherical order
    */
   int mOutM;

   /**
    * @brief Input Worland type
    */
   std::string mWType;

   /**
    * @brief Boundary condition ID
    */
   std::size_t mBcId;
};

} // namespace Sphere
} // namespace Kernel
} // namespace Spectral
} // namespace QuICC

#endif // QUICC_SPECTRAL_KERNEL_SPHERE_TRUNCATE_HPP
