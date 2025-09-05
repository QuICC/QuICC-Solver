/**
 * @file D1R1.hpp
 * @brief Implementation of d r  operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_D1R1_HPP
#define QUICC_FINITEDIFF_SPHERE_D1R1_HPP

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
    * @brief Implementation of d r operator
    */
   class D1R1: public Operator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param order   Order of accuracy
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         D1R1(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          *
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         D1R1(const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          */
         D1R1();

         /**
          * @brief Destructor
          */
         ~D1R1() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void D1R1::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, 1);

      Internal::SparseMatrix tmp = wMat.at(1)*igrid.asDiagonal();
      rOut = (this->zeroTopBottom(nR) * tmp).cast<T>();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_D1R1_HPP
