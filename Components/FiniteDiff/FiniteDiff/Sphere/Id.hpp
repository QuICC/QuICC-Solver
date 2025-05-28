/**
 * @file Id.hpp
 * @brief Implementation of identity operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_ID_HPP
#define QUICC_FINITEDIFF_SPHERE_ID_HPP

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
    * @brief Implementation of a identity operator
    */
   class Id: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         Id(const size_t order);

         /**
          * @brief Constructor
          */
         Id();

         /**
          * @brief Destructor
          */
         ~Id() = default;

         /**
          * @brief Compute operator on grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut , const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void Id::compute(Eigen::SparseMatrix<T>& rOut, const int l, const Internal::Array& igrid)
   {
      using namespace Internal::Literals;
      const int nR = igrid.size();
      rOut.resize(nR, nR);

      std::vector<Eigen::Triplet<T>> coeffs;
      for(int i = 0; i < nR; i++)
      {
         coeffs.emplace_back(i,i, 1_mp);
      }

      rOut.setFromTriplets(coeffs.begin(), coeffs.end());
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_ID_HPP
