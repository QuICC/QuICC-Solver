/**
 * @file dTorSphJnl.hpp
 * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_DTORSPHJNL_HPP
#define QUICC_POLYNOMIAL_BESSEL_DTORSPHJNL_HPP

// System includes
//
#include <boost/math/special_functions/bessel.hpp>

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

   /**
    * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
    */
   class dTorSphJnl
   {
      public:
         /**
          * @brief Default constructor
          */
         dTorSphJnl() = default;

         /**
          * @brief Compute spherical bessel basis for magnetic toroidal boundary condition
          */
         template <typename T> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale);
   };

   template <typename T> inline void dTorSphJnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int lIn, const Internal::Array& igrid, const Internal::Array& scale)
   {
      using namespace Internal::Literals;
      std::vector<Internal::MHDFloat> roots;
      Internal::MHDFloat nu = static_cast<Internal::MHDFloat>(lIn) + 0.5_mp;
      boost::math::cyl_bessel_j_zero(nu, 1, nPoly, std::back_inserter(roots));

      auto dl = static_cast<Internal::MHDFloat>(lIn);
      for(int j = 0; j < nPoly; j++)
      {
         auto k = roots.at(j);
         for(int i = 0; i < igrid.size(); i++)
         {
            const auto& r = igrid(i);

            // dJnl = l Jnl(k,l,r) / r - k Jnl(k,l+1,r)
            rOut(i,j) = (dl/r)*boost::math::sph_bessel(lIn, k*r) - k*boost::math::sph_bessel(lIn + 1, k*r);
         }
      }
   }

}
}
}

#endif // QUICC_POLYNOMIAL_BESSEL_DTORSPHJNL_HPP
