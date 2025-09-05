/**
 * @file D4.hpp
 * @brief Implementation of first derivative operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_D4_HPP
#define QUICC_FINITEDIFF_SPHERE_D4_HPP

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
    * @brief Implementation of first derivative operator
    */
   class D4: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         D4(const size_t order);

         /**
          * @brief Constructor
          */
         D4();

         /**
          * @brief Destructor
          */
         ~D4() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void D4::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, 4);

      rOut = wMat.back();

      Internal::Array qid = Internal::Array::Ones(igrid.size());
      qid(0) = 0;
      qid(nR-1) = 0;

      // Zero r = 0 and r = 1
      rOut = qid.asDiagonal() * rOut;
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_D4_HPP
