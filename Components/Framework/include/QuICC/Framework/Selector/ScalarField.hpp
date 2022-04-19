/** 
 * @file Datatypes.hpp
 * @brief Selector for implementation of generic scalar field
 */

#ifndef QUICC_FRAMEWORK_SELECTOR_DATATYPES_HPP
#define QUICC_FRAMEWORK_SELECTOR_DATATYPES_HPP

#include <variant>

#include "QuICC/ScalarFields/FlatScalarField.hpp"
#include "QuICC/Variables/Variable.hpp"
#include "QuICC/Variables/Spectral/ScalarVariable.hpp"
#include "QuICC/Variables/Spectral/VectorVariable.hpp"

namespace QuICC {

namespace Framework {

namespace Selector {

   /// Template typedef for scalar field implementation
   template <typename TData>
      using ScalarField = QuICC::Datatypes::FlatScalarField<TData>;

   /// Template typedef for shared scalar field implementation
   template <typename TData>
      using SharedScalarField = std::shared_ptr<ScalarField<TData> >;

   /// Typedef for complex valued scalar field type
   typedef ScalarField<MHDComplex> ComplexScalarField;

   /// Typedef for real valued scalar field type
   typedef ScalarField<MHDFloat> RealScalarField;

   /// Typedef for variant shared scalar field
   typedef std::variant<SharedScalarField<MHDFloat>, SharedScalarField<MHDComplex> > VariantSharedScalarField;

   // Typedef for real spectral vector field
   typedef Datatypes::VectorField<RealScalarField,FieldComponents::Spectral::Id> SpectralRealVectorField;

   // Typedef for complex spectral vector field
   typedef Datatypes::VectorField<ComplexScalarField,FieldComponents::Spectral::Id> SpectralComplexVectorField;

   // Typedef for variant shared spectral vector field
   typedef std::variant<std::shared_ptr<SpectralRealVectorField>, std::shared_ptr<SpectralComplexVectorField> > VariantSharedSpectralVectorField;

   /// Typdef for physical space scalar field type
   typedef RealScalarField PhysicalScalarField;

   /// Intermediate scalar variable type to specialize template
   template <typename T> struct ScalVarType
   {
      typedef Datatypes::Variable<Datatypes::ScalarVariable<ScalarField<T>,RealScalarField>, 1> Type;
   };

   /// Intermediate vector variable type to specialize template
   template <typename T> struct VectVarType
   {
      typedef Datatypes::Variable<Datatypes::VectorVariable<ScalarField<T>,RealScalarField>, 1> Type;
   };

   /// Templated typedef for scalar variable
   template <typename T>
      using ScalarVariable = typename ScalVarType<T>::Type;

   /// Templated typedef for shared scalar variable
   template <typename T>
      using SharedScalarVariable = typename std::shared_ptr<typename ScalVarType<T>::Type>;

   /// Templated typedef or vector variable
   template <typename T>
      using VectorVariable = typename VectVarType<T>::Type;

   /// Templated typedef or vector variable
   template <typename T>
      using SharedVectorVariable = typename std::shared_ptr<typename VectVarType<T>::Type>;

   typedef SharedScalarVariable<MHDFloat> SharedRealScalarVariable;
   typedef SharedScalarVariable<MHDComplex> SharedComplexScalarVariable;
   typedef std::variant<SharedRealScalarVariable,SharedComplexScalarVariable> VariantSharedScalarVariable;

   typedef SharedVectorVariable<MHDFloat> SharedRealVectorVariable;
   typedef SharedVectorVariable<MHDComplex> SharedComplexVectorVariable;
   typedef std::variant<SharedRealVectorVariable,SharedComplexVectorVariable> VariantSharedVectorVariable;

}
}
}

#endif // QUICC_FRAMEWORK_SELECTOR_DATATYPES_HPP
