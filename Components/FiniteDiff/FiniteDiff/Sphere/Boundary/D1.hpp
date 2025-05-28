/**
 * @file D1.hpp
 * @brief Implementation of D1 boundary operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_BOUNDARY_D1_HPP
#define QUICC_FINITEDIFF_SPHERE_BOUNDARY_D1_HPP

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
    * @brief Implementation of D1 boundary operator
    */
   class D1: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         D1(const size_t order);

         /**
          * @brief Constructor
          */
         D1();

         /**
          * @brief Destructor
          */
         ~D1() = default;

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
   void D1::compute(Eigen::SparseMatrix<T>& rOut, const int p, const int l, const Internal::Array& igrid)
   {
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
               coeffs.emplace_back(it.row(), it.col(), it.value());
            }
         }
      }

      rOut.setFromTriplets(coeffs.begin(), coeffs.end());
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_BOUNDARY_D1_HPP
