/**
 * @file PrueferAlgorithm.hpp
 * @brief Implementation of a Pruefer algorithm to compute a quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_PRUEFERALGORITHM_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_PRUEFERALGORITHM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Precision.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   /**
    * @brief Implementation of an iterative quadrature rule computation.
    *
    * This implementation relies (in general) on a coarse Runge-Kutta solver to find an initial guess. Newton
    * iterations are then used to reach high accuracy. It does not require any linear algebra solvers, removing the
    * issues coming from the lack of special matrix solvers in Eigen. Computation can be done in multiple precision.
    *
    * The implemented algorithms are based on the paper by Glaser, Liu & Rohklin, 2007. "A fast algorithm for the
    * calculation of the roots of special functions"
    */
   class PrueferAlgorithm
   {
      public:
         /**
          * @brief Number of Newton iterations in refinement loop
          */
         static const int NEWTON_ITERATIONS;

         /**
          * @brief Maximum order in the Taylor expansion
          */
         static const int TAYLOR_ORDER;

         /**
          * @brief Number of Runge-Kutta steps
          */
         static const int RK_STEPS;

      protected:
         /**
          * @brief Compute the Taylor expansion
          *
          * @param taylor  Storage for the taylor expansion
          * @param size    Grid size
          * @param u_0     Zeroth order Taylor expansion
          * @param u_1     First order Taylor expansion
          * @param xi_1    */
         void computeTaylor(internal::ArrayL& taylor, const int size, const internal::MHDLong u_0, const internal::MHDLong u_1, const internal::MHDLong xi_1);

         /**
          * @brief Runge-Kutta solver for approximations to roots
          */
         void rungekutta(internal::ArrayL& grid, const int i, const int size, const bool isZero);

         /**
          * @brief Theta expression
          */
         internal::MHDLong theta(const internal::MHDLong x, const int n);

         /**
          * @brief dX/dTheta differential equation RHS
          */
         internal::MHDLong diffeq_f(const internal::MHDLong x, const internal::MHDLong y, const int n);

         /**
          * @brief Refine node value through Newton iteration
          *
          * @param i Index of the node
          */
         void refineNode(internal::ArrayL& grid, internal::ArrayL& weights, const int i, const internal::ArrayL& taylor);

         /**
          * @brief Sort the grid and weights from the quadrature
          */
         void sortQuadrature(internal::Array& grid, internal::Array& weights);

         /**
          * @brief Get p polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   p(const internal::MHDLong xi, const int diff) = 0;

         /**
          * @brief Get q polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   q(const internal::MHDLong xi, const int diff) = 0;

         /**
          * @brief Get r polynomial
          *
          * @param size Size of the grid
          * @param diff Order of the derivative
          */
         virtual internal::MHDLong   r(const int n, const int diff) = 0;

         /**
          * @brief Compute polynomial value
          */
         virtual internal::MHDLong   u(const internal::MHDLong x, const int n);

         /**
          * @brief Compute first derivative value
          */
         virtual internal::MHDLong  du(const internal::MHDLong x, const int n);

      private:
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_PRUEFERALGORITHM_HPP
