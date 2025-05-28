/**
 * @file xLF.hpp
 * @brief ID of the generic radial + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral l ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_XLF_HPP
#define QUICC_SPATIALSCHEME_3D_XLF_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief ID of the generic radial + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral l ordering
    */
   class xLF: public ISpatialScheme
   {
      public:
         /**
          * @brief Constructor
          */
         explicit xLF(const VectorFormulation::Id formulation, const GridPurpose::Id purpose, const std::size_t id, const std::string tag, const std::string formatted);

         /**
          * @brief Destructor
          */
         ~xLF() = default;

         /**
          * @brief Set implementation type
          */
         void setImplementation(const std::map<std::size_t,std::vector<std::size_t>>& type) override;

         /**
          * @brief Create transform for given dimension
          */
         std::shared_ptr<Transform::ITransform> createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const override;

         /**
          * @brief Create index converter for given dimension
          */
         std::shared_ptr<Parallel::IIndexConv> createIndexConv(const Dimensions::Transform::Id id) const override;

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
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_XLF_HPP
