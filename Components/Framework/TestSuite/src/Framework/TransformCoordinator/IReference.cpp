/**
 * @file IReference.cpp
 * @brief Base class for input and output reference provider
 */

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/IReference.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   MHDComplex unitReference(const Test& test, const int i, const int j, const int k)
   {
      auto&& ss = test.spRes->sim().ss();
      if(ss.has(SpatialScheme::Feature::SphereGeometry) || ss.has(SpatialScheme::Feature::ShellGeometry))
      {
         if(ss.has(SpatialScheme::Feature::SpectralOrdering123))
         {
            return unitReferenceSH(i,j,k);
         }
         else
         {
            return unitReferenceSH(i,k,j);
         }
      }
      else if(ss.has(SpatialScheme::Feature::CartesianGeometry) && ss.has(SpatialScheme::Feature::FourierIndex23))
      {
         return unitReferenceFF(i,j,k);
      }
      else
      {
         throw std::logic_error("Unit spectrum for this geometry has not been implemented");
      }
   }

   MHDComplex unitReferenceSH(const int n, const int l, const int m)
   {
      MHDComplex ref(std::sqrt(2.0),-std::sqrt(2.0));

      if(l == 0)
      {
         ref = 0.0;
      }
      else if(m == 0)
      {
         ref.imag(0.0);
      }

      return ref;
   }

   MHDComplex unitReferenceFF(const int n, const int k1, const int k2)
   {
      MHDComplex ref(std::sqrt(2.0),-std::sqrt(2.0));

      if(k1 == 0 && k2 == 0)
      {
         ref.imag(0.0);
      }

      return ref;
   }
}
}
}
}
