/**
 * @file Feature.hpp
 * @brief Spatial scheme features enum
 */

#ifndef QUICC_SPATIALSCHEME_FEATURE_HPP
#define QUICC_SPATIALSCHEME_FEATURE_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Possible features for spatial schemes
    */
   enum class Feature {
      /// Regular data for spectrum
      RegularSpectrum,
      /// Spectral matrices are 1D
      SpectralMatrix1D,
      /// Spectral matrices are 2D
      SpectralMatrix2D,
      /// Spectral matrices are 3D
      SpectralMatrix3D,
      /// Spectral matrix ordering is 1,2,3 (SLFm -> R, L, M -> ordering R, L, M)
      SpectralOrdering123,
      /// Spectral matrix ordering is 1,3,2 (SLFl -> R, L, M -> ordering R, M, L)
      SpectralOrdering132,
      /// Spectral ordering in transform is 1,2,3 (SLFm -> R, L, M -> ordering R, L, M)
      TransformSpectralOrdering123,
      /// Spectral ordering in transform is 1,3,2 (SLFl -> R, L, M -> ordering R, M, L)
      TransformSpectralOrdering132,
      /// Cartesian geometry
      CartesianGeometry,
      /// Full sphere geometry
      SphereGeometry,
      /// Spherical shell geometry
      ShellGeometry,
      /// Annulus geometry
      AnnulusGeometry,
      /// Full cylinder geometry
      CylinderGeometry,
      /// Fourier as 2nd index
      FourierIndex2,
      /// Fourier as 3rd index
      FourierIndex3,
      /// Fourier as 2nd and 3rd index
      FourierIndex23,
      /// Fourier as 1st, 2nd and 3rd index
      FourierIndex123,
      /// Use galerkin basis
      GalerkinBasis,
      /// Real valued spectrum
      RealSpectrum,
      /// Complex valued spectrum
      ComplexSpectrum,
      /// Solve fourth order equation as is
      SolveFourthOrder,
      /// Split fourth order equation into 2 second order
      SplitFourthOrder,
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_FEATURE_HPP
