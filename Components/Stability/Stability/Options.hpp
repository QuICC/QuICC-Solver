/**
 * @file Options.hpp
 * @brief Small struct holding options
 */

#ifndef QUICC_STABILITY_OPTIONS_HPP
#define QUICC_STABILITY_OPTIONS_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Stability {

/**
 * @brief Options for linear stability calculations
 */
struct Options
{
   /// Tolerance for EPS solver
   double tolerance = 1e-8;

   /// Max iteration for EPS solver
   int maxIteration = 2000;

   /// Write matrices as MatrixMarket files
   bool writeMtx = false;

   /// Show verbose diagnostic
   bool verboseDiagnostics = false;

   /// Use custom initial guess
   bool useCustomGuess = false;

   /// Custom initial guess type
   int guessType = 0;

   /// scaling type:
   /// 0: Nothing special
   /// 1: Make first nonzero coefficient real and 1
   /// 2: Make first nonzero coefficient real and normalized to 1
   int scalingType = 0;

   /// Force m=0 coefficients to be real
   bool makeM0Real = false;
};

} // namespace Stability
} // namespace QuICC

#endif // QUICC_STABILITY_OPTIONS_HPP
