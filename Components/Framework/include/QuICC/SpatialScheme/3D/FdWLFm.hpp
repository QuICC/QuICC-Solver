/**
 * @file FdWLFm.hpp
 * @brief ID of the Finite Differences sphere + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
 */

#ifndef QUICC_SPATIALSCHEME_3D_FDWLFM_HPP
#define QUICC_SPATIALSCHEME_3D_FDWLFM_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/3D/xLF.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief ID of the Finite Differences sphere + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral m ordering
    */
   class FdWLFm: public xLF
   {
      public:
         /**
          * @brief Constructor
          */
         explicit FdWLFm(const VectorFormulation::Id formulation, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~FdWLFm() = default;

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

#endif // QUICC_SPATIALSCHEME_3D_FDWLFM_HPP
