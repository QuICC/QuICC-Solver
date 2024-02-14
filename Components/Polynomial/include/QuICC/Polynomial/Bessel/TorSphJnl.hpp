/**
 * @file TorSphJnl.hpp
 * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
 */

#ifndef QUICC_POLYNOMIAL_BESSEL_TORSPHJNL_HPP
#define QUICC_POLYNOMIAL_BESSEL_TORSPHJNL_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Bessel/Operators.hpp"

namespace QuICC {

namespace Polynomial {

namespace Bessel {

   /**
    * @brief Implementation of the spherical Bessel basis for magnetic toroidal boundary conditions
    */
   class TorSphJnl
   {
      public:
         /**
          * @brief Default constructor
          */
         TorSphJnl() = default;

         /**
          * @brief Compute spherical bessel basis for magnetic toroidal boundary condition
          */
         template <typename T> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int l, const Internal::Array& igrid, const Internal::Array& scale);
   };

   template <typename T> inline void TorSphJnl::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int nPoly, const int lIn, const Internal::Array& igrid, const Internal::Array& scale)
   {
      std::vector<Internal::MHDFloat> roots;
      getTorRoots(roots, lIn, nPoly);

      for(int j = 0; j < nPoly; j++)
      {
         auto k = roots.at(j);
         for(int i = 0; i < igrid.size(); i++)
         {
            rOut(i,j) = Jnl(k, lIn, igrid(i));
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

#endif // QUICC_POLYNOMIAL_BESSEL_TORSPHJNL_HPP
