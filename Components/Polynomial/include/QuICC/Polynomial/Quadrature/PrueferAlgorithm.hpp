/**
 * @file PrueferAlgorithm.hpp
 * @brief Implementation of a Pruefer algorithm to compute a quadrature rule
 */

#ifndef QUICC_POLYNOMIAL_QUADRATURE_PRUEFERALGORITHM_HPP
#define QUICC_POLYNOMIAL_QUADRATURE_PRUEFERALGORITHM_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

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
         void computeTaylor(Internal::ArrayL& taylor, const int size, const Internal::MHDLong u_0, const Internal::MHDLong u_1, const Internal::MHDLong xi_1);

         /**
          * @brief Runge-Kutta solver for approximations to roots
          */
         void rungekutta(Internal::ArrayL& grid, const int i, const int size, const bool isZero);

         /**
          * @brief Theta expression
          */
         Internal::MHDLong theta(const Internal::MHDLong x, const int n);

         /**
          * @brief dX/dTheta differential equation RHS
          */
         Internal::MHDLong diffeq_f(const Internal::MHDLong x, const Internal::MHDLong y, const int n);

         /**
          * @brief Refine node value through Newton iteration
          *
          * @param i Index of the node
          */
         void refineNode(Internal::ArrayL& grid, Internal::ArrayL& weights, const int i, const Internal::ArrayL& taylor);

         /**
          * @brief Sort the grid and weights from the quadrature
          */
         void sortQuadrature(Internal::Array& grid, Internal::Array& weights);

         /**
          * @brief Get p polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   p(const Internal::MHDLong xi, const int diff) = 0;

         /**
          * @brief Get q polynomial
          *
          * @param xi   Grid value
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   q(const Internal::MHDLong xi, const int diff) = 0;

         /**
          * @brief Get r polynomial
          *
          * @param size Size of the grid
          * @param diff Order of the derivative
          */
         virtual Internal::MHDLong   r(const int n, const int diff) = 0;

         /**
          * @brief Compute polynomial value
          */
         virtual Internal::MHDLong   u(const Internal::MHDLong x, const int n);

         /**
          * @brief Compute first derivative value
          */
         virtual Internal::MHDLong  du(const Internal::MHDLong x, const int n);

      private:
   };

}
}
}

#endif // QUICC_POLYNOMIAL_QUADRATURE_PRUEFERALGORITHM_HPP
