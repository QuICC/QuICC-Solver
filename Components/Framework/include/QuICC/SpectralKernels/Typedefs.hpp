/**
 * @file Set3DModes.hpp
 * @brief Set given spectral modes
 */

#ifndef QUICC_SPECTRAL_KERNEL_TYPEDEFS_HPP
#define QUICC_SPECTRAL_KERNEL_TYPEDEFS_HPP

// First include
//

// Configuration includes
//
#include <memory>
#include <map>

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   /// Typedef to simplify notations for real valued fast index mode map
   typedef std::map<int,MHDFloat> RealFastMapType;

   /// Typedef to simplify notations for complex valued fast index mode map
   typedef std::map<int,MHDComplex> ComplexFastMapType;

   /// Typedef to simplify notations for 3D mode
   typedef std::map<std::pair<int, int>, RealFastMapType > Real3DMapType;

   /// Typedef to simplify notations for 3D mode map
   typedef std::map<std::pair<int, int>, ComplexFastMapType > Complex3DMapType;

   /// Typedef to simplify notations for a shared real 3D mode map
   typedef std::shared_ptr<Real3DMapType> SharedReal3DMapType;

   /// Typedef to simplify notations for a shared complex 3D mode map
   typedef std::shared_ptr<Complex3DMapType> SharedComplex3DMapType;

   /// Typedef to simplify notations for a component map for shared real modes map
   typedef std::map<FieldComponents::Spectral::Id, SharedReal3DMapType> VectSharedR3DMapType;

   /// Typedef to simplify notations for a component map for shared complex modes map
   typedef std::map<FieldComponents::Spectral::Id, SharedComplex3DMapType> VectSharedZ3DMapType;
}
}
}

#endif // QUICC_SPECTRAL_KERNEL_TYPEDEFS_HPP
