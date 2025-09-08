/**
 * @file Overr1.hpp
 * @brief Implementation of division by R operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_OVERR1_HPP
#define QUICC_FINITEDIFF_SPHERE_OVERR1_HPP

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

   /**
    * @brief Implementation of division by R operator
    */
   class Overr1: public Operator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param order   Order of accuracy
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         Overr1(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          *
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         Overr1(const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          */
         Overr1();

         /**
          * @brief Destructor
          */
         ~Overr1() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void Overr1::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      int i0 = std::max(static_cast<int>(this->mZtop), 1);
      std::vector<Eigen::Triplet<Internal::MHDFloat>> coeffs;
      for(int i = i0; i < nR - static_cast<int>(this->mZbot); i++)
      {
         coeffs.emplace_back(i,i, 1_mp/igrid(i));
      }

      Internal::SparseMatrix tmp(nR, nR);
      tmp.setFromTriplets(coeffs.begin(), coeffs.end());

      if(l == 1 && this->mZtop == 0)
      {
         std::vector<Internal::SparseMatrix> wMat;
         this->fdMatrices(wMat, igrid, this->mOrder, 1);
         tmp += this->zeroTopBottom(nR, 0, nR-1)*wMat.back();
      }

      rOut = tmp.cast<T>();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_OVERR1_HPP
