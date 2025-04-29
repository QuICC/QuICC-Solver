/**
 * @file SLFm.hpp
 * @brief ID of the shell Chebyshev(FFT) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_SLFM_HPP
#define QUICC_SPATIALSCHEME_3D_SLFM_HPP

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
    * @brief ID of the shell Chebyshev(FFT) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
    */
   class SLFm: public ISpatialScheme
   {
      public:
         /**
          * @brief Constructor
          */
         explicit SLFm(const VectorFormulation::Id formulation, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~SLFm() = default;

         /**
          * @brief Unique id
          */
         static const std::size_t sId;

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

#endif // QUICC_SPATIALSCHEME_3D_SLFM_HPP
