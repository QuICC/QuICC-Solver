/**
 * @file Integral.hpp
 * @brief Implementation of integral operator
 */

#ifndef QUICC_FINITEDIFF_SPHERE_INTEGRAL_HPP
#define QUICC_FINITEDIFF_SPHERE_INTEGRAL_HPP

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
    * @brief Implementation of integral operator
    */
   class Integral: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         Integral(const size_t order);

         /**
          * @brief Constructor
          */
         Integral();

         /**
          * @brief Destructor
          */
         ~Integral() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut , const int l, const Internal::Array& igrid);

      private:
         /**
          * @brief Order of Taylor expansion
          */
         const unsigned int mcTaylorN;

   };

   template <typename T>
   void Integral::compute(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      assert(rOut.rows() == 1);
      assert(rOut.cols() == nR);

      std::vector<Internal::SparseMatrix> wMat;
      this->fdMatrices(wMat, igrid, this->mOrder, this->mcTaylorN);

      // Compute grid delta: x_{n+1} - x_{n}
      Internal::Array dgp = -igrid;
      dgp.topRows(igrid.size()-1) += igrid.bottomRows(igrid.size()-1);
      dgp.bottomRows(1).setZero();
      // Compute grid delta: x_{n-1} - x_{n}
      Internal::Array dgm = -igrid;
      dgm.bottomRows(igrid.size()-1) += igrid.topRows(igrid.size()-1);
      dgm.topRows(1).setZero();

      Internal::Array integral = Internal::Array::Zero(igrid.size());
      int c = 2;
      for(int i = 0; i <= static_cast<int>(this->mcTaylorN); i++)
      {
         c *= (i+1);
         integral += (1_mp/static_cast<Internal::MHDFloat>(c))*(dgp.array().pow(i+1) - dgm.array().pow(i+1)).matrix().transpose()*wMat.at(i);
      }
      rOut = integral.cast<T>().transpose();
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_INTEGRAL_HPP
