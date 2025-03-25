/**
 * @file Options.hpp
 * @brief Small struct holding options
 */

#ifndef QUICC_STABILITY_OPTIONS_HPP
#define QUICC_STABILITY_OPTIONS_HPP

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Stability {

/**
 * @brief Options for linear stability calculations
 */
struct Options
{
   /// Solver mode
   int solver_mode = 0;

   /// Eigen solver type
   /// 0: Krylov-Schur
   /// 1: CISS
   int eigensolver_type = 0;

   /// Number of eigenvalues to compute
   std::size_t nev = 3;

   /// Tolerance for EPS solver
   double tolerance = 1e-8;

   /// Max iteration for EPS solver
   int maxIteration = 2000;

   /// Write original matrices as MatrixMarket files
   bool writeMtx = false;

   /// Write PETSc matrices as binary Petsc files
   bool writePetsc = false;

   /// Show verbose diagnostic
   bool verboseDiagnostics = false;

   /// Use custom initial guess
   bool useCustomGuess = false;

   /// Custom initial guess type
   int guessType = 0;

   /// scaling type:
   /// 0: Nothing special
   /// 1: Make coefficient with largest amplitude real positive
   /// 2: Make coefficient with largest amplitude 1
   int scalingType = 1;

   /// List of fields to get reference from
   std::vector<std::pair<std::size_t, FieldComponents::Spectral::Id>>
      scalingRef;

   /// Force m=0 coefficients to be real
   bool makeM0Real = false;
};

} // namespace Stability
} // namespace QuICC

#endif // QUICC_STABILITY_OPTIONS_HPP
