/**
 * @file D3.hpp
 * @brief Implementation of first derivative operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_D3_HPP
#define QUICC_FINITEDIFF_SPHERE_D3_HPP

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
   class D3: public Operator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param order   Order of accuracy
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         D3(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          *
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         D3(const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          */
         D3();

         /**
          * @brief Destructor
          */
         ~D3() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void D3::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, 3);

      int zTop = std::max(static_cast<int>(this->mZtop), static_cast<int>(l == 0 || l == 2 || l > 3));
      int zBot = static_cast<int>(this->mZbot);
      rOut = (this->zeroTopBottom(nR, zTop, zBot)*wMat.back()).cast<T>();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_D3_HPP
