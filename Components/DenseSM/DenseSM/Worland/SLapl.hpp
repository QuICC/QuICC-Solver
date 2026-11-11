/**
 * @file SLapl.hpp
 * @brief Implementation of the spherical laplacian tau operator
 */

#ifndef QUICC_DENSESM_WORLAND_SLAPL_HPP
#define QUICC_DENSESM_WORLAND_SLAPL_HPP

// System includes
//

// Project includes
//
#include "DenseSM/Worland/ITauOperator.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

/**
 * @brief Implementation of the dense spherical laplacian spectral operator
 */
class SLapl : public ITauOperator
{
public:
   enum class OrthoId: std::size_t {
      Worland = 0,
      ScaValueSphEnergyWorland,
      TorValueSphEnergyWorland,
      PolNoSlipSphEnergyWorland,
      PolInsulatingSphEnergyWorland,
      ScaValueSphEnergy,
      TorValueSphEnergy,
      PolNoSlipSphEnergy,
      PolInsulatingSphEnergy
   };

   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of columns
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    * @param l       Harmonic degree l
    * @param bId     Basis ID
    */
   SLapl(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l, const std::size_t bId);

   /**
    * @brief Constructor
    *
    * @param rows    Number of rows
    * @param cols    Number of columns
    * @param alpha   Jacobi alpha parameter
    * @param dBeta   Jacobi dBeta parameter
    * @param l       Harmonic degree l
    * @param bId     Basis ID
    */
   SLapl(const int rows, const int cols, const Scalar_t alpha,
      const Scalar_t dBeta, const int l, const OrthoId bId);

   /**
    * @brief Destructor
    */
   virtual ~SLapl() = default;

protected:
   /**
    * @brief Implementation of build dense matrix operator
    * @param output operator
    */
   virtual void buildOpImpl(Internal::Matrix& mat, const int rows,
      const int cols) const override;

   /**
    * @brief Basis base Jacobi alpha
    */
   static Scalar_t basisAlpha(const std::size_t bId);

   /**
    * @brief Basis base Jacobi alpha
    */
   static Scalar_t basisDBeta(const std::size_t bId);

private:
   /**
    * @brief Build dense tau operator
    */
   void buildTauOp(Internal::Matrix& mat, const int rows,
      const int cols) const;

   /**
    * @brief Normalization factor
    */
   Scalar_t norm(const OrthoId id, const int n, const int l) const;

   /**
    * @brief c1 expansion coefficient for ValueD1 PolSphEnergy basis
    */
   Scalar_t c1(const int n, const int l) const;

   /**
    * @brief c2 expansion coefficient for  ValueD1 PolSphEnergy basis
    */
   Scalar_t c2(const int n, const int l) const;

   /**
    * @brief c3 expansion coefficient for  ValueD1 PolSphEnergy basis
    */
   Scalar_t c3(const int n, const int l) const;

   /**
    * @brief d1 expansion coefficient for Insulating PolSphEnergy basis
    */
   Scalar_t d1(const int n, const int l) const;

   /**
    * @brief d2 expansion coefficient for Insulating PolSphEnergy basis
    */
   Scalar_t d2(const int n, const int l) const;
};

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC

#endif // QUICC_DENSESM_WORLAND_SLAPL_HPP
