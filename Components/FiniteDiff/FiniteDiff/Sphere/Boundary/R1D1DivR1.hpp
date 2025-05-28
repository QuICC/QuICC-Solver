/**
 * @file R1D1DivR1.hpp
 * @brief Implementation of R1D1DivR1 boundary operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_BOUNDARY_R1D1DIVR1_HPP
#define QUICC_FINITEDIFF_SPHERE_BOUNDARY_R1D1DIVR1_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Literals.hpp"
#include "FiniteDiff/Sphere/Operator.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   /**
    * @brief Implementation of value boundary operator operator
    */
   class R1D1DivR1: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         R1D1DivR1(const size_t order);

         /**
          * @brief Constructor
          */
         R1D1DivR1();

         /**
          * @brief Destructor
          */
         ~R1D1DivR1() = default;

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
   void R1D1DivR1::compute(Eigen::SparseMatrix<T>& rOut, const int p, const int l, const Internal::Array& igrid)
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
               // - 1/r
               if(it.col() == p)
               {
                  v -= 1_mp/igrid(p);
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

#endif // QUICC_FINITEDIFF_SPHERE_BOUNDARY_R1D1DIVR1_HPP
