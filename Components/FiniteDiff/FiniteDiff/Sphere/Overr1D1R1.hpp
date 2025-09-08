/**
 * @file Overr1D1R1.hpp
 * @brief Implementation of 1/r d r  operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_OVERR1D1R1_HPP
#define QUICC_FINITEDIFF_SPHERE_OVERR1D1R1_HPP

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
    * @brief Implementation of 1/r d r operator
    */
   class Overr1D1R1: public Operator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param order   Order of accuracy
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         Overr1D1R1(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          *
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         Overr1D1R1(const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          */
         Overr1D1R1();

         /**
          * @brief Destructor
          */
         ~Overr1D1R1() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void Overr1D1R1::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, 1);

      Internal::Array invgrid = Internal::Array::Zero(igrid.size());
      for(int i = 1; i < igrid.size(); i++)
      {
         invgrid(i) = 1_mp/igrid(i);
      }

      Internal::SparseMatrix tmp = wMat.back();
      tmp += invgrid.asDiagonal();

      if(l == 1 && this->mZtop == 0)
      {
         tmp += this->zeroTopBottom(nR, 0, nR-1)*wMat.back();
      }

      int zTop = std::max(static_cast<int>(this->mZtop), static_cast<int>(l == 0 || l > 1));
      int zBot = static_cast<int>(this->mZbot);
      rOut = (this->zeroTopBottom(nR, zTop, zBot)*tmp).cast<T>();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_OVERR1D1R1_HPP
