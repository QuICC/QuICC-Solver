/** 
 * @file Pnab.hpp
 * @brief Implementation of the Jacobi polynomial
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_PNAB_HPP
#define QUICC_POLYNOMIAL_JACOBI_PNAB_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Precision.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Jacobi/JacobiBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Jacobi {

/**
 * @brief Implementation of the Jacobi polynomial
 */
template <typename NormType = Quadrature::natural_t>
class Pnab
{
   public:
      /**
       * @brief Compute \f$P_n^{(\alpha,\beta)} (x)\f$
       *
       * Internal computation can be done in multiple precision
       */
      template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator)
      {
         int gN = igrid.rows();

         if (alpha < MHD_MP(-1.0) || beta < MHD_MP(-1.0))
         {
            throw std::logic_error("Tried to compute Jacobi polynomial P_n^{(alpha,beta)} with alpha < -1 or beta < -1");
         }

         if (nPoly < 1)
         {
            throw std::logic_error("Operator matrix should have at least 1 column");
         }

         if (gN != igrid.size())
         {
            throw std::logic_error("Operator matrix does not mach grid size");
         }

         internal::Matrix ipoly(gN,2);

         ThreeTermRecurrence::P0(ipoly.col(0), alpha, beta, igrid,
            &JacobiBase::P0ab<normalizationTag_t>);
         if(scale.size() > 0)
         {
            ipoly.col(0).array() *= scale.array();
         }
         evaluator(rOut, ipoly.col(0), 0);

         if(nPoly > 0)
         {
            ThreeTermRecurrence::P1(ipoly.col(1), alpha, beta, igrid,
               &JacobiBase::P1ab<normalizationTag_t>);
            evaluator(rOut, ipoly.col(1), 1);
         }

         for(int i = 2; i < nPoly; ++i)
         {
            ThreeTermRecurrence::Pn(ipoly.col(0), i, alpha, beta, ipoly.col(1), ipoly.col(0), igrid,
               &JacobiBase::Pnab<normalizationTag_t>);
            ipoly.col(0).swap(ipoly.col(1));
            evaluator(rOut, ipoly.col(1), i);
         }
      }

   private:
      using normalizationTag_t = NormType;
};

} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_JACOBI_PNAB_HPP
