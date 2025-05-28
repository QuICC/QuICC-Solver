/**
 * @file InsulatingSphere.hpp
 * @brief Implementation of InsulatingSphere boundary operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_BOUNDARY_INSULATINGSPHERE_HPP
#define QUICC_FINITEDIFF_SPHERE_BOUNDARY_INSULATINGSPHERE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Literals.hpp"
#include "FiniteDiff/Sphere/Operator.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   /**
    * @brief Implementation of a InsulatingSphere boundary operator
    */
   class InsulatingSphere: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         InsulatingSphere(const size_t order);

         /**
          * @brief Constructor
          */
         InsulatingSphere();

         /**
          * @brief Destructor
          */
         ~InsulatingSphere() = default;

         /**
          * @brief Compute operator on grid
          *
          * @param rOut    Sparse operator
          * @param p       Position of boundary condition
          * @param l       Harmonic degree l
          * @param igrid   Radial grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut, const int p, const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void InsulatingSphere::compute(Eigen::SparseMatrix<T>& rOut, const int p, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, 1);
      auto& dMat = wMat.at(1);

      std::vector<Eigen::Triplet<T>> coeffs;
      for(int k = 0; k < dMat.outerSize(); ++k)
      {
         for(Internal::SparseMatrix::InnerIterator it(dMat,k); it; ++it)
         {
            if(it.row() == p)
            {
               // \partial_r
               auto v = it.value();
               // + (l+1)/r
               if(it.col() == p)
               {
                  v += static_cast<Internal::MHDFloat>(l + 1)/igrid(p);
               }
               coeffs.emplace_back(it.row(), it.col(), v);
            }
         }
      }

      rOut.setFromTriplets(coeffs.begin(), coeffs.end());
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_BOUNDARY_INSULATINGSPHERE_HPP
