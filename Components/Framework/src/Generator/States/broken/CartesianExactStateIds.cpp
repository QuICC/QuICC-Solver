/**
 * @file CartesianExactStateIds.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Generator/States/CartesianExactStateIds.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Equations {

   const MHDFloat CartesianExactStateIds::PCOS = 99999;

   const MHDFloat CartesianExactStateIds::PSIN = -99999;

   MHDFloat CartesianExactStateIds::cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::cos(mode*theta);
   }

   MHDFloat CartesianExactStateIds::sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta)
   {
      return amplitude*std::sin(mode*theta);
   }

   MHDFloat CartesianExactStateIds::poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      MHDFloat val;

      if(mode == CartesianExactStateIds::PCOS)
      {
         val = CartesianExactStateIds::cos(amplitude,mode,Math::PI*(x-1.0)/2.0);
      } else if(mode == CartesianExactStateIds::PSIN)
      {
         val = CartesianExactStateIds::sin(amplitude,mode,Math::PI*(x-1.0)/2.0);
      } else
      {
         val = amplitude*std::pow(x,mode);
      }

      return val;
   }

   MHDFloat CartesianExactStateIds::chebyshev(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      return amplitude*std::cos(mode*std::acos(x));
   }

   MHDFloat CartesianExactStateIds::zero(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x)
   {
      MHDFloat val;

      if(mode == CartesianExactStateIds::PCOS)
      {
         val = CartesianExactStateIds::cos(amplitude,mode,Math::PI*x/2.0);
      } else if(mode == CartesianExactStateIds::PSIN)
      {
         val = CartesianExactStateIds::sin(amplitude,mode,Math::PI*(x-1.0)/2.0);
      } else
      {
         val = amplitude*std::pow((1-x*x),mode);
      }

      return val;
   }

   MHDFloat CartesianExactStateIds::exact2D(const CartesianExactStateIds::Id id, const Array& amplitude, const Array& mode, const Array& x)
   {
      MHDFloat val;

      if(id == POLYPOLY)
      {
         val = CartesianExactStateIds::poly(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::poly(amplitude(1), mode(1), x(1));
      } else if(id == POLYCOS)
      {
         val = CartesianExactStateIds::poly(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::cos(amplitude(1), mode(1), x(1));
      } else if(id == POLYSIN)
      {
         val = CartesianExactStateIds::poly(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::sin(amplitude(1), mode(1), x(1));
      } else if(id == COSCOS)
      {
         val = CartesianExactStateIds::cos(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::cos(amplitude(1), mode(1), x(1));
      } else if(id == SINSIN)
      {
         val = CartesianExactStateIds::sin(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::sin(amplitude(1), mode(1), x(1));
      } else if(id == COSSIN)
      {
         val = CartesianExactStateIds::cos(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::sin(amplitude(1), mode(1), x(1));
      } else if(id == SINCOS)
      {
         val = CartesianExactStateIds::sin(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::cos(amplitude(1), mode(1), x(1));
      } else if(id == ZEROCOS)
      {
         val = CartesianExactStateIds::zero(amplitude(0), mode(0), x(0));
         val *= CartesianExactStateIds::cos(amplitude(1), mode(1), x(1));
      } else
      {
         throw std::logic_error("Unknown exact 2D state");
      }

      return val;
   }

   MHDFloat CartesianExactStateIds::exact3D(const CartesianExactStateIds::Id id, const Array& amplitude, const Array& mode, const Array& x)
   {
      MHDFloat val;

      if(id == POLYPOLYPOLY)
      {
         val = CartesianExactStateIds::exact2D(POLYPOLY, amplitude, mode, x);
         val *= CartesianExactStateIds::poly(amplitude(2), mode(2), x(2));
      } else if(id == POLYCOSPOLY)
      {
         val = CartesianExactStateIds::exact2D(POLYCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::poly(amplitude(2), mode(2), x(2));
      } else if(id == POLYSINPOLY)
      {
         val = CartesianExactStateIds::exact2D(POLYSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::poly(amplitude(2), mode(2), x(2));
      } else if(id == POLYCOSCOS)
      {
         val = CartesianExactStateIds::exact2D(POLYCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else if(id == POLYSINSIN)
      {
         val = CartesianExactStateIds::exact2D(POLYSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::sin(amplitude(2), mode(2), x(2));
      } else if(id == POLYSINCOS)
      {
         val = CartesianExactStateIds::exact2D(POLYSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else if(id == POLYCOSSIN)
      {
         val = CartesianExactStateIds::exact2D(POLYCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::sin(amplitude(2), mode(2), x(2));
      } else if(id == COSCOSCOS)
      {
         val = CartesianExactStateIds::exact2D(COSCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else if(id == SINSINSIN)
      {
         val = CartesianExactStateIds::exact2D(SINSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::sin(amplitude(2), mode(2), x(2));
      } else if(id == COSSINCOS)
      {
         val = CartesianExactStateIds::exact2D(COSSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else if(id == SINCOSSIN)
      {
         val = CartesianExactStateIds::exact2D(SINCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::sin(amplitude(2), mode(2), x(2));
      } else if(id == SINSINCOS)
      {
         val = CartesianExactStateIds::exact2D(SINSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else if(id == COSCOSSIN)
      {
         val = CartesianExactStateIds::exact2D(COSCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::sin(amplitude(2), mode(2), x(2));
      } else if(id == SINCOSCOS)
      {
         val = CartesianExactStateIds::exact2D(SINCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else if(id == COSSINSIN)
      {
         val = CartesianExactStateIds::exact2D(COSSIN, amplitude, mode, x);
         val *= CartesianExactStateIds::sin(amplitude(2), mode(2), x(2));
      } else if(id == ZEROCOSCOS)
      {
         val = CartesianExactStateIds::exact2D(ZEROCOS, amplitude, mode, x);
         val *= CartesianExactStateIds::cos(amplitude(2), mode(2), x(2));
      } else
      {
         throw std::logic_error("Unknown exact 3D state");
      }

      return val;
   }

}
}
