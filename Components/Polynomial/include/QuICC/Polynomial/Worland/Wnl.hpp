/**
 * @file Wnl.hpp
 * @brief Implementation of the Worland polynomial
 */

#ifndef QUICC_POLYNOMIAL_WORLAND_WNL_HPP
#define QUICC_POLYNOMIAL_WORLAND_WNL_HPP

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
#include "Types/Precision.hpp"
#include "QuICC/Polynomial/ThreeTermRecurrence.hpp"
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   /**
    * @brief Implementation of the Worland polynomial
    */
   class Wnl: public WorlandBase
   {
      public:
         /**
          * @brief Default constructorx
          */
         Wnl() = default;

         /**
          * @brief Constructor for specific alpha,beta pair
          */
         Wnl(const internal::MHDFloat alpha, const internal::MHDFloat dBeta);

         /**
          * @brief Compute worland polynomial
          *
          * @tparam TEvaluator The evaluator allows to change behavior from computing Matric operator, to On-the-fly transforms, etc
          */
         template <typename T, typename TEvaluator> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEvaluator evaluator);

   };

   template <typename T, typename TEval> inline void Wnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const internal::Array& igrid, const internal::Array& scale, TEval evaluator)
   {
      int gN = igrid.rows();

      if(l < 0)
      {
         throw std::logic_error("Tried to compute Worland polynomial W_n^l with l < 0");
      }

      if(nPoly < 1)
      {
         throw std::logic_error("Operator matrix should have at least 1 column");
      }

      if(gN != igrid.size())
      {
         throw std::logic_error("Operator matrix does not mach grid size");
      }

      internal::Matrix ipoly(gN,2);

      internal::MHDFloat a = this->alpha(l);
      internal::MHDFloat b = this->beta(l);

      this->computeW0l(ipoly.col(0), l, a, b, igrid, WorlandBase::normWP0ab());
      if(scale.size() > 0)
      {
         ipoly.col(0).array() *= scale.array();
      }
      evaluator(rOut, ipoly.col(0), 0);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), a, b, ipoly.col(0), ixgrid, WorlandBase::normWP1ab());
         evaluator(rOut, ipoly.col(1), 1);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(0), i, a, b, ipoly.col(1), ipoly.col(0), ixgrid, WorlandBase::normWPnab());
         ipoly.col(0).swap(ipoly.col(1));
         evaluator(rOut, ipoly.col(1), i);
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_WORLAND_WNL_HPP
