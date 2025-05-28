/**
 * @file Overr1.hpp
 * @brief Implementation of division by R operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_OVERRR1_HPP
#define QUICC_FINITEDIFF_SPHERE_OVERRR1_HPP

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
          */
         Overr1(const size_t order);

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

      std::vector<Eigen::Triplet<T>> coeffs;
      coeffs.emplace_back(0,0, 0_mp);
      for(int i = 1; i < nR; i++)
      {
         coeffs.emplace_back(i,i, 1_mp/igrid(i));
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

#endif // QUICC_FINITEDIFF_SPHERE_OVERRR1_HPP
