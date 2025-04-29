/**
 * @file CylinderExactStateIds.hpp
 * @brief Class to hold the list of possible exact states
 */

#ifndef QUICC_EQUATIONS_CYLINDEREXACTSTATEIDS_HPP
#define QUICC_EQUATIONS_CYLINDEREXACTSTATEIDS_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Class to hold the list of possible exact states
    */
   struct CylinderExactStateIds
   {
      /// Polynomial approximation to Cosine
      static const MHDFloat PCOS;

      /// Polynomial approximation to Sine
      static const MHDFloat PSIN;

      /**
       * @brief Enums for the avaialable exact states
       */
      enum Id {
         NOTUSED = -1,  // Initialisation state (do NOT use)
         // Special states
         CONSTANT = 0,  // All constant
         POLYCOSPOLY = 10, // Polynomial, Cosine, Polynomial
         POLYSINPOLY,      // Polynomial, Sine, Polynomial
         SPEC_UNIT, // Spectral space solution: unit spectrum of max 5 modes
         // ---------------------------------------------------
         // Test states
         // ---------------------------------------------------
         TEST_DIFFUSION,   // Initial state for diffusion test
         TEST_BIDIFFUSION_DIRECT,      // Initial state for bidiffusion test (direct solve)
         TEST_BIDIFFUSION_INFLUENCE,   // Initial state for diffusion test (influence matrix solve)
         TEST_BIDIFFUSION_SPLIT,       // Initial state for diffusion test (split solve)
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
   };

}
}

#endif // QUICC_EQUATIONS_CYLINDEREXACTSTATEIDS_HPP
