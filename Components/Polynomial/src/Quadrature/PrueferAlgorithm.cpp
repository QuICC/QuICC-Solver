/**
 * @file PrueferAlgorithm.cpp
 * @brief Source of the implementation of the Pruefer algorithm for the computation of a quadrature rule
 */

// Configuration includes
//

// System includes
//
#include <map>
#include <stdexcept>
#include <iostream>
#include <iomanip>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/Quadrature/PrueferAlgorithm.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

#ifdef QUICC_MULTPRECISION
   const int PrueferAlgorithm::NEWTON_ITERATIONS = 10;

   const int PrueferAlgorithm::TAYLOR_ORDER = 60;

   const int PrueferAlgorithm::RK_STEPS = 20;
#else
   const int PrueferAlgorithm::NEWTON_ITERATIONS = 10;

   const int PrueferAlgorithm::TAYLOR_ORDER = 30;

   const int PrueferAlgorithm::RK_STEPS = 20;
#endif //QUICC_MULTPRECISION

   void PrueferAlgorithm::refineNode(Internal::ArrayL& grid, Internal::ArrayL& weights, const int i, const Internal::ArrayL& taylor)
   {
      // Storage for (x_{k+1} - x_{k})
      Internal::MHDLong   h;
      // Storage for the poylnomial value
      Internal::MHDLong   f;
      // Constants to break out of Newton iterations
      auto epsilon = std::numeric_limits<Internal::MHDLong>::epsilon();
      static constexpr int ulp = 2;

      // Compute Newton refinement iterations
      auto grid_i = grid(i);
      const auto grid_i_mone = grid(i-1);
      for(int n = 0; n < PrueferAlgorithm::NEWTON_ITERATIONS; n++)
      {
         // Initialise starting values
         h = MHD_MP_LONG(1.0);
         f = MHD_MP_LONG(0.0);
         weights(i) = MHD_MP_LONG(0.0);

         // Loop over all Taylor coefficients
         for(int k = 0; k < taylor.size()-1; k++)
         {

            // Add next order to function value
            f += taylor(k)*h;

            // Add next order to derivative value
            weights(i) += taylor(k+1)*h;

            // Increment the monomial (including factorial part)
            h *= (grid_i - grid_i_mone)/static_cast<Internal::MHDLong>(k+1);
         }
         // Add last order to function value
         f += taylor(taylor.size()-1)*h;

         auto delta = f/weights(i);

         // Check
         auto tol = epsilon * Internal::Math::abs(grid_i) * ulp;
         if( Internal::Math::abs(delta) <= tol )
         {
            break;
         }

         // update the Newton iteration value
         grid_i -= delta;
      }
      // store
      grid(i) = grid_i;

   }

   void PrueferAlgorithm::sortQuadrature(Internal::Array& grid, Internal::Array& weights)
   {
      // Safety assert
      assert(grid.size() == weights.size());

      // Create a map to sort elements
      std::map<Internal::MHDFloat, Internal::MHDFloat> sorter;

      // fill map with grid/weights pairs
      for(int i = 0; i < grid.size(); ++i)
      {
         sorter.insert(std::make_pair(grid(i), weights(i)));
      }

      // Check that no point got lost
      if(sorter.size() != static_cast<size_t>(grid.size()))
      {
         throw std::logic_error("PrueferAlgorithm::sortQuadrature: Lost grid points during sorting!");
      }

      // Replace grid point values with reordered version
      int i = 0;
      for(auto it = sorter.cbegin(); it != sorter.cend(); ++it, ++i)
      {
         grid(i) = it->first;
         weights(i) = it->second;
      }
   }

   void PrueferAlgorithm::computeTaylor(Internal::ArrayL& taylor, const int size, const Internal::MHDLong u_0, const Internal::MHDLong u_1, const Internal::MHDLong xi)
   {
      // Make sure to reset to zero
      taylor.setZero();

      // Fill zeroth and first order
      taylor(0) = u_0;
      taylor(1) = u_1;

      for(int k = 0; k < taylor.size()-2; k++)
      {
         Internal::MHDLong dk = static_cast<Internal::MHDLong>(k);
         Internal::MHDLong dkk_1 = static_cast<Internal::MHDLong>(k*(k-1));

         // First term
         taylor(k+2) = -(dk*this->p(xi,1) + this->q(xi,0))*taylor(k+1);

         // Second term
         taylor(k+2) -= ((dkk_1/MHD_MP_LONG(2.0))*this->p(xi,2) + dk*this->q(xi,1) + this->r(size,0))*taylor(k);

         // Third term (if applicable)
         if(k > 0)
         {
            taylor(k+2) -= ((dkk_1/MHD_MP_LONG(2.0))*this->q(xi,2) + dk*this->r(size,1))*taylor(k-1);
         }

         // Fourth term (if applicable)
         if(k > 1)
         {
            taylor(k+2) -= ((dkk_1/MHD_MP_LONG(2.0))*this->r(size,2))*taylor(k-2);
         }

         taylor(k+2) /= this->p(xi,0);
      }
   }

   Internal::MHDLong PrueferAlgorithm::theta(const Internal::MHDLong x, const int n)
   {
      return Internal::Math::atan(MHD_MP_LONG(1.0)/Internal::Math::sqrt(this->r(n,0)*this->p(x,0))*(this->p(x,0)*this->du(x,n)/this->u(x,n)));
   }

   Internal::MHDLong PrueferAlgorithm::diffeq_f(const Internal::MHDLong x, const Internal::MHDLong y, const int n)
   {
      return -MHD_MP_LONG(1.0)/(Internal::Math::sqrt(this->r(n,0)/this->p(y,0)) + ((this->r(n,1)*this->p(y,0) - this->p(y,1)*this->r(n,0) + MHD_MP_LONG(2.0)*this->r(n,0)*this->q(y,0))/(MHD_MP_LONG(2.0)*this->r(n,0)*this->p(y,0)))*(Internal::Math::sin(MHD_MP_LONG(2.0)*x)/MHD_MP_LONG(2.0)));
   }

   void PrueferAlgorithm::rungekutta(Internal::ArrayL& grid, const int i, const int n, const bool isZero)
   {
      Internal::MHDLong y = grid(i);
      Internal::MHDLong x;

      if(isZero)
      {
         x = Internal::Math::PI_long/MHD_MP_LONG(2.0);
      } else
      {
         x = this->theta(y, n);
      }

      Internal::MHDLong h = (-Internal::Math::PI_long/MHD_MP_LONG(2.0)-x)/static_cast<Internal::MHDLong>(PrueferAlgorithm::RK_STEPS);
      Internal::MHDLong kold,k;
      kold = h*this->diffeq_f(x,y,n);
      for(int i = 0; i < PrueferAlgorithm::RK_STEPS; i++)
      {
         x += h;
         k = h*this->diffeq_f(x,y + kold,n);
         y += MHD_MP_LONG(0.5)*(kold + k);
         kold = k;
      }

      assert(!Internal::Math::isnan(y));

      grid(i) = y;
   }

   Internal::MHDLong PrueferAlgorithm::u(const Internal::MHDLong, const int)
   {
      throw std::logic_error("PrueferAlgorithm::u needs to be implemented in quadrature rule");
   }

   Internal::MHDLong PrueferAlgorithm::du(const Internal::MHDLong, const int)
   {
      throw std::logic_error("PrueferAlgorithm::u needs to be implemented in quadrature rule");
   }
}
}
}
