/** 
 * @file d3Pnab.hpp
 * @brief Implementation of the D^3 Jacobi polynomial
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_D3PNAB_HPP
#define QUICC_POLYNOMIAL_JACOBI_D3PNAB_HPP

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
 * @brief Implementation of the D^3 Jacobi polynomial
 */
template <typename NormType = Quadrature::natural_t>
class d3Pnab
{
   public:
      /**
       * @brief Compute \f$\frac{d}{dx}P_n^{(\alpha,\beta)} (x)\f$
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

         internal::Matrix idiff(gN,2);

         internal::MHDFloat a3 = alpha + MHD_MP(3.0);
         internal::MHDFloat b3 = beta + MHD_MP(3.0);

         idiff.col(0).setZero();
         evaluator(rOut, idiff.col(0), 0);

         if(nPoly > 1)
         {
            idiff.col(1).setZero();
            evaluator(rOut, idiff.col(1), 1);
         }

         if(nPoly > 2)
         {
            idiff.col(0).setZero();
            idiff.col(0).swap(idiff.col(1));
            evaluator(rOut, idiff.col(1), 2);
         }

         if(nPoly > 3)
         {
            ThreeTermRecurrence::P0(idiff.col(0), a3, b3, igrid, &JacobiBase::d3P0ab<normalizationTag_t>);
            if(scale.size() > 0)
            {
               idiff.col(0).array() *= scale.array();
            }
            idiff.col(0).swap(idiff.col(1));
            evaluator(rOut, idiff.col(1), 3);
         }

         if(nPoly > 4)
         {
            ThreeTermRecurrence::P1(idiff.col(0), a3, b3, idiff.col(1), igrid, &JacobiBase::d3P1ab<normalizationTag_t>);
            idiff.col(0).swap(idiff.col(1));
            evaluator(rOut, idiff.col(1), 4);
         }

         for(int i = 5; i < nPoly; ++i)
         {
            ThreeTermRecurrence::Pn(idiff.col(0), i-3, a3, b3, idiff.col(1), idiff.col(0), igrid, &JacobiBase::d3Pnab<normalizationTag_t>);
            idiff.col(0).swap(idiff.col(1));
            evaluator(rOut, idiff.col(1), i);
         }
      }
   private:
      using normalizationTag_t = NormType;

};

} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_JACOBI_D3PNAB_HPP
