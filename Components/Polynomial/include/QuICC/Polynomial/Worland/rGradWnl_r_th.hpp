/**
 * @file rGradWnl_r_th.hpp
 * @brief Implementation of the radial part of the tensor gradient component r,th
 * @brief Worland polynomials
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_RGRADWNL_R_TH_HPP
#define QUICC_POLYNOMIAL_WORLAND_RGRADWNL_R_TH_HPP

// System includes
//

// Project includes
//
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "Types/Internal/Literals.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

/**
 * @brief Implementation of the 1/r D Worland polynomial with recurrence
 * relation
 */
class rGradWnl_r_th: public WorlandBase
{
public:
   /**
    * @brief Default constructor
    */
   rGradWnl_r_th() = default;

   /**
    * @brief Constructor for specific alpha,beta pair
    *
    * @param alpha   Jacobi alpha
    * @param dBeta   Jacobi beta = l + dBeta
    */
   rGradWnl_r_th(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta) :
       WorlandBase(alpha, dBeta) {};

   template <typename T, typename TEvaluator>
   void compute(
      Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
      const int nPoly, const int l, const Internal::Array& igrid,
      const Internal::Array& scale, TEvaluator evaluator);
};

template <typename T, typename TEvaluator>
void rGradWnl_r_th::compute(
   Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> rOut,
   const int nPoly, const int l, const Internal::Array& igrid,
   const Internal::Array& scale, TEvaluator evaluator)
{
   using namespace Internal::Literals;
   int gN = igrid.rows();

   if (l < 0)
   {
      throw std::logic_error(
         "Tried to compute Worland 1/r D operator with l < 0");
   }

   if (nPoly < 1)
   {
      throw std::logic_error("Operator matrix should have at least 1 column");
   }

   if (gN != igrid.size())
   {
      throw std::logic_error("Operator matrix does not mach grid size");
   }

   Internal::MHDFloat a = this->alpha(l);
   Internal::MHDFloat b = this->beta(l);
   Internal::MHDFloat a1 = this->alpha(l) + 1.0_mp;
   Internal::MHDFloat b1 = this->beta(l) + 1.0_mp;
   Internal::MHDFloat dl = Internal::MHDFloat(l);

   // Make X grid in [-1, 1]
   Internal::Array ixgrid = 2.0_mp * igrid.array() * igrid.array() - 1.0_mp;

   // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
   Internal::Matrix ipnab(gN, 2);
   Internal::Matrix idpnab(gN, 2);

   // Compute P_0
   this->computeW0l(ipnab.col(0), l - 2, a, b, igrid, WorlandBase::normWP0ab());
   ipnab.col(0).array() *= (dl-1)*(dl+1);
   if (scale.size() > 0)
   {
      ipnab.col(0).array() *= scale.array();
   }

   // Compute DP_0
   idpnab.col(0).setZero();

   // Compute l P
   evaluator(rOut, ipnab.col(0), 0);

   if (nPoly > 1)
   {
      // Compute P_0
      ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid,
         WorlandBase::normWP1ab());

      // Compute DP_1
      this->computeW0l(idpnab.col(0), l, a1, b1, igrid,
         WorlandBase::normWDP0ab());
      if (scale.size() > 0)
      {
         idpnab.col(0).array() *= scale.array();
      }

      // Compute e P + 4r^2 DP
      evaluator(rOut, ipnab.col(1) - idpnab.col(0), 1);
   }

   if (nPoly > 2)
   {
      // Increment P_n
      ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0),
         ixgrid, WorlandBase::normWPnab());
      ipnab.col(0).swap(ipnab.col(1));

      // Compute DP_2
      ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid,
         WorlandBase::normWDP1ab());

      // Compute e P + 2(x+1) DP
      evaluator(rOut, ipnab.col(1) - idpnab.col(1), 2);
   }

   for (int i = 3; i < nPoly; ++i)
   {
      // Increment P_n
      ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0),
         ixgrid, WorlandBase::normWPnab());
      ipnab.col(0).swap(ipnab.col(1));

      // Increment DP_n
      ThreeTermRecurrence::Pn(idpnab.col(0), i - 1, a1, b1, idpnab.col(1),
         idpnab.col(0), ixgrid, WorlandBase::normWDPnab());
      idpnab.col(0).swap(idpnab.col(1));

      // Compute e P + 2(x+1) DP
      evaluator(rOut, ipnab.col(1) - idpnab.col(1), i);
   }
}

} // namespace Worland
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_WORLAND_RGRADWNL_R_TH_HPP
