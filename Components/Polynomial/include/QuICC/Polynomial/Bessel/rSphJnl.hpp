/**
 * @file rSphJnl.hpp
 * @brief Implementation of the spherical Bessel basis
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_RSPHJNL_HPP
#define QUICC_POLYNOMIAL_BESSEL_RSPHJNL_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Bessel/details/Operators.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

   /**
    * @brief Implementation of the spherical Bessel basis
    */
   class rSphJnl
   {
      public:
         /**
          * @brief Needed additional modes for computation
          */
         static const int EXTRA_POLY = 0;

         /**
          * @brief Additional roots with different harmonic degree
          */
         static const int EXTRA_L = 0;

         /**
          * @brief Default constructor
          */
         rSphJnl() = default;

         /**
          * @brief Compute spherical bessel basis
          */
         template <typename T> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const std::vector<Internal::MHDFloat>& roots, const int l, const Internal::Array& igrid, const Internal::Array& scale);
   };

   template <typename T> inline void rSphJnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const std::vector<Internal::MHDFloat>& roots, const int l, const Internal::Array& igrid, const Internal::Array& scale)
   {
      const int nPoly = roots.size();

      for(int j = 0; j < nPoly; j++)
      {
         auto k = roots.at(j);
         for(int i = 0; i < igrid.size(); i++)
         {
            rOut(i,j) = details::rSphJnl(k, l, igrid(i));
         }

         if(scale.size() > 0)
         {
            rOut.col(j).array() *= scale.array();
         }
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_BESSEL_RSPHJNL_HPP
