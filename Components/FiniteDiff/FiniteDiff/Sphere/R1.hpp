/**
 * @file R1.hpp
 * @brief Implementation of multiplication by R operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_R1_HPP
#define QUICC_FINITEDIFF_SPHERE_R1_HPP

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
    * @brief Implementation of multiplication by R operator
    */
   class R1: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         R1(const size_t order);

         /**
          * @brief Constructor
          */
         R1();

         /**
          * @brief Destructor
          */
         ~R1() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void R1::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Eigen::Triplet<T>> coeffs;
      for(int i = 0; i < nR; i++)
      {
         coeffs.emplace_back(i,i, igrid(i));
      }

      rOut.setFromTriplets(coeffs.begin(), coeffs.end());

      // Zero r = 0 and r = 1
      Internal::Array qid = Internal::Array::Ones(igrid.size());
      qid(0) = 0;
      qid(nR-1) = 0;
      rOut = qid.asDiagonal() * rOut;
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_R1_HPP
