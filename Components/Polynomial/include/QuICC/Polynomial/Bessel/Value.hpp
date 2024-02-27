/**
 * @file Value.hpp
 * @brief Implementation of the spherical Bessel basis for zero value boundary conditions
 */

#ifndef _VALUE_
#define _VALUE_

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
    * @brief Implementation of the spherical Bessel basis for zero value boundary conditions
    */
   template <typename TOp> class Value: private TOp
   {
      public:
         /**
          * @brief Default constructor
          */
         Value() = default;

         /**
          * @brief Compute spherical bessel basis for magnetic toroidal boundary condition
          */
         template <typename T> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale);
   };

   template <typename TOp> template <typename T> inline void Value<TOp>::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale)
   {
      std::vector<Internal::MHDFloat> roots;
      details::getRoots(roots, l, nPoly + TOp::EXTRA_POLY, details::Value_dNu());

      if constexpr(TOp::EXTRA_L == 0)
      {
         TOp::compute(rOut, roots, l, igrid, scale, details::Value_dNu());
      }
      else
      {
         std::vector<Internal::MHDFloat> roots_extra;
         details::getRoots(roots_extra, l + TOp::EXTRA_L, nPoly + TOp::EXTRA_POLY, details::Value_dNu());

         TOp::compute(rOut, roots, roots_extra, l, igrid, scale, details::Value_dNu());
      }
   }

}
}
}

#endif // _VALUE_
