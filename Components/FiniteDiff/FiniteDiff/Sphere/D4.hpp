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
          *
          * @param order   Order of accuracy
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         D4(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          *
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         D4(const std::size_t zTop, const std::size_t zBot);

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

      rOut = (this->zeroTopBottom(nR)*wMat.back()).cast<T>();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_D4_HPP
