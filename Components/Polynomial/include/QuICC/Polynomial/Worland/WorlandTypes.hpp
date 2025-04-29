/**
 * @file WorlandTypes.hpp
 * @brief Worland types
 */
#pragma once

// External includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandChebyshevRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandLegendreRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandCylEnergyRule.hpp"
#include "QuICC/Polynomial/Quadrature/WorlandSphEnergyRule.hpp"

namespace QuICC {
namespace Polynomial {
namespace Worland {

//
// Worland types
//

/// @brief Chebyshev type
struct worland_chebyshev_t
{
   /// Typedef for quadrature rule
   typedef Polynomial::Quadrature::WorlandChebyshevRule Rule;

   /// Jacobi alpha parameter
   static const Internal::MHDFloat ALPHA;
   /// Jacobi beta = l + dBeta parameter
   static const Internal::MHDFloat DBETA;
};

/// @brief Legendre type
struct worland_legendre_t
{
   /// Typedef for quadrature rule
   typedef Polynomial::Quadrature::WorlandLegendreRule Rule;

   /// Jacobi alpha parameter
   static const Internal::MHDFloat ALPHA;
   /// Jacobi beta = l + dBeta parameter
   static const Internal::MHDFloat DBETA;
};

/// @brief CylEnergy type
struct worland_cylenergy_t
{
   /// Typedef for quadrature rule
   typedef Polynomial::Quadrature::WorlandCylEnergyRule Rule;

   /// Jacobi alpha parameter
   static const Internal::MHDFloat ALPHA;
   /// Jacobi beta = l + dBeta parameter
   static const Internal::MHDFloat DBETA;
};

/// @brief SphEnergy type
struct worland_sphenergy_t
{
   /// Typedef for quadrature rule
   typedef Polynomial::Quadrature::WorlandSphEnergyRule Rule;

   /// Jacobi alpha parameter
   static const Internal::MHDFloat ALPHA;
   /// Jacobi beta = l + dBeta parameter
   static const Internal::MHDFloat DBETA;
};

/// @brief Default type selected at CMake setup
struct worland_default_t
{
   /// Typedef for quadrature rule
   typedef Polynomial::Quadrature::WorlandRule Rule;

   /// Jacobi alpha parameter
   static const Internal::MHDFloat ALPHA;
   /// Jacobi beta = l + dBeta parameter
   static const Internal::MHDFloat DBETA;
};

} // namespace Worland
} // namespace Polynomial
} // namespace QuICC
