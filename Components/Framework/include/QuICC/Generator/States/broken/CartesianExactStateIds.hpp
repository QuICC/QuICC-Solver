/**
 * @file CartesianExactStateIds.hpp
 * @brief Class to hold the list of possible exact states
 */

#ifndef QUICC_EQUATIONS_CARTESIANEXACTSTATEIDS_HPP
#define QUICC_EQUATIONS_CARTESIANEXACTSTATEIDS_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Class to hold the list of possible exact states
    */
   struct CartesianExactStateIds
   {
      /// Polynomial approximation to Cosine
      static const MHDFloat PCOS;

      /// Polynomial approximation to Sine
      static const MHDFloat PSIN;

      /**
       * @brief Enums for the avaialable exact states
       */
      enum Id {
         // Special states
         CONSTANT = 0,  // All constant
         // ---------------------------------------------------
         // 3D modes
         // ---------------------------------------------------
         // TTT states
         POLYPOLYPOLY = 10,  // Polynomial, Polynomial, Polynomial
         // TFT states
         POLYCOSPOLY = 20, // Polynomial, Cosine, Polynomial
         POLYSINPOLY,      // Polynomial, Sine, Polynomial
         // TFF states
         POLYCOSCOS = 30,  // Polynomial, Cosine, Cosine
         POLYSINSIN,       // Polynomial, Sine, Sine
         POLYSINCOS,       // Polynomial, Sine, Cosine
         POLYCOSSIN,       // Polynomial, Cosine, Sine
         BXHELICOIDAL,     // Helicoidal Bx (Stellmach & Hansen, 2004)
         BYHELICOIDAL,     // Helicoidal By (Stellmach & Hansen, 2004)
         NULLFIELD,
         CONSTANTFIELD,
         // FFF states
         COSCOSCOS = 50,   // Cosine, Cosine, Cosine
         SINSINSIN,        // Sine, Sine, Sine
         COSCOSSIN,        // Cosine, Cosine, Sine
         SINSINCOS,        // Sine, Sine, Cosine
         COSSINSIN,        // Cosine, Sine, Sine
         SINCOSCOS,        // Sine, Cosine, Cosine
         COSSINCOS,        // Cosine, Sine, Cosine
         SINCOSSIN,        // Sine, Cosine, Sine
         // ---------------------------------------------------
         // Special 3D states (benchmarks, etc)
         // ---------------------------------------------------
         ZEROCOSCOS = 90,  // Galerkin, Cosine, Cosine
         PEYRET1DA,        // Kind of a place holder for tests
         TORPOLTFF,        // Divergence free state for Toroidal/Poloidal test in TFF scheme
         TORPOLCNST,       // Unit field for Toroidal/Poloidal test in TFF scheme
         PLANFORMSQUARES,  // Planform squares
         // ---------------------------------------------------
         // Test states
         // ---------------------------------------------------
         TEST_DIFFUSION,         // Initial state for diffusion test (T = T' = 0 z=+-1)
         TEST_BIDIFFUSION,       // Initial state for bidiffusion test (T = T'= Lapl T = 0 z=+-1)
         TEST_BIDIFFUSION_SPLIT, // Initial state for bidiffusion split (\phi = lapl T)
         // ---------------------------------------------------
         // 2D modes
         // ---------------------------------------------------
         // TT states
         POLYPOLY = 110,  // Polynomial, Polynomial
         // TF states
         POLYCOS = 120, // Polynomial, Cosine
         POLYSIN,       // Polynomial, Sine
         // FF states
         COSCOS = 130,  // Cosine, Cosine
         SINSIN,        // Sine, Sine
         COSSIN,        // Cosine, Sine
         SINCOS,        // Sine, Cosine
         // ---------------------------------------------------
         // Special 2D states (benchmarks, etc)
         // ---------------------------------------------------
         ZEROCOS = 190,  // Kind of a place holder for tests
      };

      /**
       * @brief Compute even periodic mode
       */
      static MHDFloat cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta);

      /**
       * @brief Compute odd periodic mode
       */
      static MHDFloat sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta);

      /**
       * @brief Compute polynomial mode
       */
      static MHDFloat poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x);

      /**
       * @brief Compute Chebyshev mode
       */
      static MHDFloat chebyshev(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x);

      /**
       * @brief Compute galerkin polynomial mode (zero at boundary)
       */
      static MHDFloat zero(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x);

      /**
       * @brief Compute 2D mode
       */
      static MHDFloat exact2D(const Id id, const Array& amplitude, const Array& mode, const Array& x);

      /**
       * @brief Compute 3D mode
       */
      static MHDFloat exact3D(const Id id, const Array& amplitude, const Array& mode, const Array& x);
   };

}
}

#endif // QUICC_EQUATIONS_CARTESIANEXACTSTATEIDS_HPP
