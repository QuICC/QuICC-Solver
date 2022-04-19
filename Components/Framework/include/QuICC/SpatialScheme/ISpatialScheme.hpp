/**
 * @file ISpatialScheme.hpp
 * @brief Interface for a spatial scheme
 */

#ifndef QUICC_SPATIALSCHEME_ISPATIALSCHEME_HPP
#define QUICC_SPATIALSCHEME_ISPATIALSCHEME_HPP

// System includes
//
#include <string>
#include <memory>
#include <variant>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Tools/ComponentAlias.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"

// forward declarations
namespace QuICC {

namespace Transform {

   class ITransform;
   class TransformSetup;
}

namespace Parallel {

   class IIndexConv;
}

namespace Equations {

namespace Tools {

   class ICoupling;
}
}

class Resolution;

namespace SpatialScheme {

   class IBuilder;
}
}


namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Interface for a spatial scheme name 
    */
   class ISpatialScheme
   {
      public:
         /// Typedef for real transform data type
         typedef Framework::Selector::RealScalarField RealTransformDataType;
         /// Typedef for complex transform data type
         typedef Framework::Selector::ComplexScalarField ComplexTransformDataType;
         /// Typedef for real/complex variant transform data type
         typedef std::variant<RealTransformDataType *,ComplexTransformDataType *> VariantTransformDataPointer;
         /// Typedef for scalar variable
         typedef Framework::Selector::VariantSharedScalarVariable ScalarVariable;
         /// Typedef for complex transform data type
         typedef Framework::Selector::VariantSharedVectorVariable VectorVariable;

         /**
          * @brief Constructor
          */
         ISpatialScheme(const VectorFormulation::Id formulation, const GridPurpose::Id purpose, const int dimension, const std::size_t id, const std::string tag, const std::string formatted);

         /**
          * @brief Destructor
          */
         virtual ~ISpatialScheme();

         /**
          * @brief Formulation used for vector fields
          */
         VectorFormulation::Id formulation() const;

         /**
          * @brief Main purpose of the grid values
          */
         GridPurpose::Id purpose() const;

         /**
          * @brief ID of the spatial scheme
          */
         std::size_t id() const;

         /**
          * @brief Dimension of the spatial scheme
          */
         int dimension() const;

         /**
          * @brief Tag of the spatial scheme
          */
         std::string tag() const;

         /**
          * @brief Formatted name of the spatial scheme
          */
         std::string formatted() const;

         /**
          * @brief Physical field component aliases
          */
         const Tools::ComponentAlias<FieldComponents::Physical::Id>& physical() const;

         /**
          * @brief Physical field component aliases
          */
         const Tools::ComponentAlias<FieldComponents::Spectral::Id>& spectral() const;

         /**
          * @brief Check for additional features
          */
         virtual bool has(const Feature f) const;

         /**
          * @brief Enable additional feature
          */
         virtual void enable(const Feature f);

         /**
          * @brief Enable additional feature
          */
         virtual void enable(const std::set<Feature>& f);

         /**
          * @brief Disable  feature
          */
         virtual void disable(const Feature f);

         /**
          * @brief Enable additional feature
          */
         virtual const std::set<Feature>& features() const;

         /**
          * @brief Set implementation type
          */
         void setImplementation(const ArrayI& type);

         /**
          * @brief Create the scheme builder
          */
         virtual std::shared_ptr<IBuilder> createBuilder(ArrayI& dim, const bool needInterpretation) const = 0;

         /**
          * @brief Create transform for given dimension
          */
         virtual std::shared_ptr<Transform::ITransform> createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const = 0;

         /**
          * @brief Create transform for given dimension
          */
         virtual std::shared_ptr<Parallel::IIndexConv> createIndexConv(const Dimensions::Transform::Id id) const = 0;

         /**
          * @brief Create spectral decomposition tools
          */
         virtual std::shared_ptr<Equations::Tools::ICoupling> createCouplingTools(const Equations::CouplingIndexType indexType) const = 0;

         /**
          * @brief Get variant forward transform data type with correct type ininitialized
          */
         virtual VariantTransformDataPointer fwdPtr(const Dimensions::Transform::Id id) const = 0;

         /**
          * @brief Get variant forward transform data type with correct type ininitialized
          */
         virtual VariantTransformDataPointer bwdPtr(const Dimensions::Transform::Id id) const = 0;

         /**
          * @brief Create variant scalar variable
          */
         virtual ScalarVariable createSVar(std::shared_ptr<Resolution> spRes) const = 0;

         /**
          * @brief Create variant vector variable
          */
         virtual VectorVariable createVVar(std::shared_ptr<Resolution> spRes) const = 0;
         
      protected:
         /**
          * @brief Implementation type of transforms
          */
         ArrayI mImplType;

         /**
          * @brief Physical component aliases
          */
         Tools::ComponentAlias<FieldComponents::Physical::Id>  mPhys;

         /**
          * @brief Physical component aliases
          */
         Tools::ComponentAlias<FieldComponents::Spectral::Id>  mSpec;

      private:
         /**
          * @brief Formulation used for vector fields
          */
         const VectorFormulation::Id mFormulation;

         /**
          * @brief Main purpose of the grid values
          */
         const GridPurpose::Id mPurpose;

         /**
          * @brief Dimension of spatial scheme
          */
         const int mDimension;

         /**
          * @brief ID of spatial scheme
          */
         const size_t mId;

         /**
          * @brief Tag of spatial scheme
          */
         const std::string mTag;

         /**
          * @brief Formatted name of spatial scheme
          */
         const std::string mFormatted;

         /**
          * @brief Enabled features
          */
         std::set<Feature> mFeatures;
   };

   /// Typedef for shared_pointer spatial scheme
   typedef std::shared_ptr<ISpatialScheme> SharedISpatialScheme;

   /// Typedef for shared_pointer spatial scheme
   typedef std::shared_ptr<const ISpatialScheme> SharedCISpatialScheme;
}
}

#endif // QUICC_SPATIALSCHEME_ISPATIALSCHEME_HPP
