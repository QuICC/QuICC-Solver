/**
 * @file Value.hpp
 * @brief Implementation of Value boundary operator 
 */

#ifndef QUICC_FINITEDIFF_SPHERE_BOUNDARY_VALUE_HPP
#define QUICC_FINITEDIFF_SPHERE_BOUNDARY_VALUE_HPP

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

namespace Boundary {

   /**
    * @brief Implementation of value boundary operator operator
    */
   class Value: public Operator
   {
      public:
         /**
          * @brief Constructor
          */
         Value(const size_t order);

         /**
          * @brief Constructor
          */
         Value();

         /**
          * @brief Destructor
          */
         ~Value() = default;

         /**
          * @brief Compute operator on grid
          *
          * @param rOut    Sparse operator
          * @param p       Position of boundary condition
          * @param l       Harmonic degree l
          * @param igrid   Radial grid
          */
         template <typename T> void compute(Eigen::SparseMatrix<T>& rOut, const int p, const int l, const Internal::Array& igrid);

      private:

   };

   template <typename T>
   void Value::compute(Eigen::SparseMatrix<T>& rOut, const int p, const int l, const Internal::Array& igrid)
   {
      std::vector<Eigen::Triplet<T>> coeffs;
      coeffs.emplace_back(p, p, 1);

      rOut.setFromTriplets(coeffs.begin(), coeffs.end());
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC

#endif // QUICC_FINITEDIFF_SPHERE_BOUNDARY_VALUE_HPP
