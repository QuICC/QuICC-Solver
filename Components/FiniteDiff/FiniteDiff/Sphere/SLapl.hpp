/**
 * @file SLapl.hpp
 * @brief Implementation of spherical laplacian operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_SLAPL_HPP
#define QUICC_FINITEDIFF_SPHERE_SLAPL_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "Types/Internal/Literals.hpp"
#include "FiniteDiff/Sphere/Operator.hpp"
#include "FiniteDiff/Sphere/UniformRadialGrid.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   /**
    * @brief Implementation of spherical laplacian operator
    */
   class SLapl: public Operator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param order   Order of accuracy
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         SLapl(const size_t order, const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          *
          * @param zTop    Zero rows at top
          * @param zBot    Zero rows at bottom
          */
         SLapl(const std::size_t zTop, const std::size_t zBot);

         /**
          * @brief Constructor
          */
         SLapl();

         /**
          * @brief Destructor
          */
         ~SLapl() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void SLapl::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, 2);

      Internal::Array invgrid = Internal::Array::Zero(igrid.size());
      for(int i = 1; i < igrid.size(); i++)
      {
         invgrid(i) = 1_mp/igrid(i);
      }

      Internal::SparseMatrix tmp = wMat.at(2);
      tmp += 2_mp*invgrid.asDiagonal()*wMat.at(1);
      tmp -= static_cast<Internal::MHDFloat>(l*(l+1))*invgrid.array().pow(2).matrix().asDiagonal();

      int zTop = std::max(static_cast<int>(this->mZtop), static_cast<int>(l > 0));
      int zBot = static_cast<int>(this->mZbot);
      rOut = (this->zeroTopBottom(nR, zTop, zBot)*tmp).cast<T>();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_SLAPL_HPP
