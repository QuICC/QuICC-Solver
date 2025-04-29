/**
 * @file dPnab.hpp
 * @brief Implementation of the D Jacobi polynomial
 */

#ifndef QUICC_POLYNOMIAL_JACOBI_DPNAB_HPP
#define QUICC_POLYNOMIAL_JACOBI_DPNAB_HPP

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
#include "Types/Internal/BasicTypes.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Jacobi/JacobiBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Jacobi {

/**
 * @brief Implementation of the D Jacobi polynomial
 */
template <typename NormType = Quadrature::natural_t>
class dPnab
{
   public:
      /**
       * @brief Compute \f$\frac{d}{dx}P_n^{(\alpha,\beta)} (x)\f$
       *
       * Internal computation can be done in multiple precision
       */
      template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, const Internal::Array& scale, TEvaluator evaluator)
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

         Internal::Matrix idiff(gN,3);

         idiff.col(0).setZero();
         evaluator(rOut, idiff.col(0), 0);

         ThreeTermRecurrence::P0(idiff.col(0), alpha, beta, igrid, &JacobiBase::dP1ab<normalizationTag_t>);
         if(scale.size() > 0)
         {
            idiff.col(0).array() *= scale.array();
         }
         evaluator(rOut, idiff.col(0), 1);

         if(nPoly > 1)
         {
            ThreeTermRecurrence::P1(idiff.col(1), alpha, beta, igrid, &JacobiBase::dP2ab<normalizationTag_t>);
            if(scale.size() > 0)
            {
               idiff.col(1).array() *= scale.array();
            }
            evaluator(rOut, idiff.col(1), 2);
         }

         // Karniadakis and Sherwin page 586
         ThreeTermRecurrence::P0(idiff.col(0), alpha, beta, igrid, &JacobiBase::P0ab<normalizationTag_t>);
         ThreeTermRecurrence::P1(idiff.col(1), alpha, beta, igrid, &JacobiBase::P1ab<normalizationTag_t>);
         // P2
         ThreeTermRecurrence::Pn(idiff.col(0), 2, alpha, beta, idiff.col(1), idiff.col(0), igrid,
               &JacobiBase::Pnab<normalizationTag_t>);
         idiff.col(0).swap(idiff.col(1));

         for(int i = 3; i < nPoly+1; ++i)
         {
            // Pi
            ThreeTermRecurrence::Pn(idiff.col(0), i, alpha, beta, idiff.col(1), idiff.col(0), igrid,
               &JacobiBase::Pnab<normalizationTag_t>);
            idiff.col(0).swap(idiff.col(1));
            ThreeTermRecurrence::dPn(idiff.col(2), i, alpha, beta, idiff.col(1), idiff.col(0), igrid, &JacobiBase::dPnabPnab<normalizationTag_t>);
            if(scale.size() > 0)
            {
               idiff.col(2).array() *= scale.array();
            }
            evaluator(rOut, idiff.col(2), i);
         }
      }

   private:
      using normalizationTag_t = NormType;
   };

} // namespace Jacobi
} // namespace Polynomial
} // namespace QuICC

#endif // QUICC_POLYNOMIAL_JACOBI_DPNAB_HPP
