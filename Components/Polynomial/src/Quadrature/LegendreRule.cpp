/** 
 * @file LegendreRule.cpp
 * @brief Source of the Legendre quadrature
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/PrueferAlgorithm.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   void LegendreRule::computeQuadrature(internal::Array& igrid, internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      igrid.setZero();
      iweights.resize(size);
      iweights.setZero();

      internal::ArrayL ig(size);
      ig.setZero();
      internal::ArrayL iw(size);
      iw.setZero();

      internal::ArrayL   taylor(std::min(size+1,PrueferAlgorithm::TAYLOR_ORDER+1));

      // Need to treat differently odd or even sizes
      if(size % 2 == 0)
      {
         // Start from extremum value 
         ig(0) = MHD_MP_LONG(0.0);
         ig(1) = this->estimateNode(0, size);
         // Compute accurate node and derivative
         this->computeTaylor(taylor, size, this->zeroPoly(size), MHD_MP_LONG(0.0), ig(0));
         this->refineNode(ig, iw, 1, taylor);

         // Set obtained node as first node and its derivative
         ig(0) = ig(1);
         iw(0) = iw(1);

         // Start filling opposite end
         ig(size-1) = -ig(0);
         iw(size-1) = -iw(0);

      } else
      {
         // 0 is a node
         ig(0) = MHD_MP_LONG(0.0);
         iw(0) = this->zeroDiff(size);
      }

      // Compute grid
      for(int i = 1; i < (size+1)/2; i++)
      {
         // Estimate node position with formulae
         ig(i) = this->estimateNode(i - size%2, size);
         this->computeTaylor(taylor, size, MHD_MP_LONG(0.0), iw(i-1), ig(i-1));
         this->refineNode(ig, iw, i, taylor);

         // If solution is too far from estimate, redo full loop starting from solution
         // This should only be required for "small" number of grid points (estimate is asymptotic formulae)
         if(precision::abs((ig(i) - this->estimateNode(i - size%2, size))/ig(i)) > 1.0e-8)
         {
            this->computeTaylor(taylor, size, MHD_MP_LONG(0.0), iw(i-1), ig(i-1));
            this->refineNode(ig, iw, i, taylor);
         }

         // Fill other half of symmetric nodes
         ig(size - 1 - i + (size%2)) = -ig(i);
         iw(size - 1 - i + (size%2)) = iw(i);
      }

      // Convert derivative to weights
      igrid = ig.cast<internal::MHDFloat>();
      iweights = (MHD_MP_LONG(2.0)*((MHD_MP_LONG(1.0)-ig.array().square()).array()*iw.array().square()).inverse()).cast<internal::MHDFloat>();

      // Sort the grid and weights
      this->sortQuadrature(igrid, iweights);
   }

   internal::MHDLong   LegendreRule::p(const internal::MHDLong xi, const int diff)
   {
      // Safety asserts
      assert(diff >= 0);

      // Get p polynomial
      if(diff == 0)
      {
         return MHD_MP_LONG(1.0)-xi*xi;

      // Get first derivative of p polynomial
      } else if(diff == 1)
      {
         return MHD_MP_LONG(-2.0)*xi;

      // Get second derivative of p polynomial
      } else if(diff == 2)
      {
         return MHD_MP_LONG(-2.0);

      } else
      {
         return MHD_MP_LONG(0.0);
      }
   }

   internal::MHDLong   LegendreRule::q(const internal::MHDLong xi, const int diff)
   {
      // Safety asserts
      assert(diff >= 0);

      // Get q polynomial
      if(diff == 0)
      {
         return MHD_MP_LONG(-2.0)*xi;

      // Get first derivative of q polynomial
      } else if(diff == 1)
      {
         return MHD_MP_LONG(-2.0);

      // Get second derivative of q polynomial
      } else if(diff == 2)
      {
         return MHD_MP_LONG(0.0);

      } else
      {
         return MHD_MP_LONG(0.0);
      }
   }

   internal::MHDLong   LegendreRule::r(const int size, const int diff)
   {
      // Safety asserts
      assert(diff >= 0);

      // Get r polynomial
      if(diff == 0)
      {
         return static_cast<internal::MHDLong>(size*(size+1));

      // Get first derivative of r polynomial
      } else if(diff == 1)
      {
         return MHD_MP_LONG(0.0);

      // Get second derivative of r polynomial
      } else if(diff == 2)
      {
         return MHD_MP_LONG(0.0);

      } else
      {
         return MHD_MP_LONG(0.0);
      }
   }

   internal::MHDLong LegendreRule::estimateNode(const int k, const int size)
   {
      // Storage for the node estimate
      internal::MHDLong   x;

      // Cast grid size to floating value
      internal::MHDLong   rN = static_cast<internal::MHDLong>(size);

      // Get theta value
      internal::MHDLong   th = static_cast<internal::MHDLong>(4*((size/2)-k)-1)/(MHD_MP_LONG(4.0)*rN+MHD_MP_LONG(2.0))*Precision::PI_long;

      x = (MHD_MP_LONG(1.0) - (rN-MHD_MP_LONG(1.0))/(MHD_MP_LONG(8.0)*rN*rN*rN) - MHD_MP_LONG(1.0)/(MHD_MP_LONG(384.0)*rN*rN*rN*rN)*(MHD_MP_LONG(39.0)-MHD_MP_LONG(28.0)/(precision::sin(th)*precision::sin(th))))*precision::cos(th);

      return x;
   }

   internal::MHDLong LegendreRule::zeroPoly(const int size)
   {
      // Initialise start value of recurrence
      internal::MHDLong p = MHD_MP_LONG(1.0);

      for(int i = 0; i < size/2; i++)
      {
         p *= -static_cast<internal::MHDLong>(2*i+1)/static_cast<internal::MHDLong>((2*i+1) + 1);
      }

      return p;
   }

   internal::MHDLong LegendreRule::zeroDiff(const int size)
   {
      // Initialise start value of recurrence
      internal::MHDLong p = MHD_MP_LONG(1.0);

      // Initialise start value of recurrence
      internal::MHDLong dp = MHD_MP_LONG(1.0);

      for(int i = 0; i < size/2; i++)
      {
         p *= -static_cast<internal::MHDLong>(2*i+1)/static_cast<internal::MHDLong>((2*i+1) + 1);

         dp *= -static_cast<internal::MHDLong>(2*i+2)/static_cast<internal::MHDLong>((2*i+2) + 1);
         dp += static_cast<internal::MHDLong>(2*(2*i+2)+1)/static_cast<internal::MHDLong>((2*i+2) + 1)*p;
      }

      return dp;
   }

}
}
}
