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
          */
         SLapl(const size_t order);

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

      rOut = wMat.at(2);
      rOut += 2_mp*igrid.array().pow(-1).matrix().asDiagonal()*wMat.at(1);
      rOut -= static_cast<Internal::MHDFloat>(l*(l+1))*igrid.array().pow(-2).matrix().asDiagonal();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_SLAPL_HPP
