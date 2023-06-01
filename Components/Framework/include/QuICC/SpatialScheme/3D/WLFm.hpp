/**
 * @file WLFm.hpp
 * @brief ID of the sphere Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_WLFM_HPP
#define QUICC_SPATIALSCHEME_3D_WLFM_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief ID of the sphere Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
    */
   class WLFm: public ISpatialScheme
   {
      public:
         /**
          * @brief Constructor
          */
         explicit WLFm(const VectorFormulation::Id formulation, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~WLFm() = default;

         /**
          * @brief Unique id
          */
         static const std::size_t sId;

         /**
          * @brief Set implementation type
          */
         void setImplementation(const std::map<std::size_t,std::vector<std::size_t>>& type) final;

         /**
          * @brief Create the scheme builder
          */
         std::shared_ptr<IBuilder> createBuilder(ArrayI& dim, const bool needInterpretation) const final;

         /**
          * @brief Create transform for given dimension
          */
         std::shared_ptr<Transform::ITransform> createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const final;

         /**
          * @brief Create index converter for given dimension
          */
         std::shared_ptr<Parallel::IIndexConv> createIndexConv(const Dimensions::Transform::Id id) const final;

         /**
          * @brief Create spectral decomposition tools
          */
         std::shared_ptr<Equations::Tools::ICoupling> createCouplingTools(const Equations::CouplingIndexType indexType) const final;

         /**
          * @brief Get variant forward transform data type with correct type ininitialized
          */
         VariantTransformDataPointer fwdPtr(const Dimensions::Transform::Id id) const final;

         /**
          * @brief Get variant forward transform data type with correct type ininitialized
          */
         VariantTransformDataPointer bwdPtr(const Dimensions::Transform::Id id) const final;

         /**
          * @brief Create variant scalar variable
          */
         ScalarVariable createSVar(std::shared_ptr<Resolution> spRes) const final;

         /**
          * @brief Create variant vector variable
          */
         VectorVariable createVVar(std::shared_ptr<Resolution> spRes) const final;

      protected:

      private:
         /**
          * @brief Unique tag
          */
         static const std::string sTag;

         /**
          * @brief Formatted name
          */
         static const std::string sFormatted;
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_WLFM_HPP
