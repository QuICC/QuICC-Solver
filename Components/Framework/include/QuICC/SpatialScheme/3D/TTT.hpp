/**
 * @file TTT.hpp
 * @brief ID of the Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme
 */

#ifndef QUICC_SPATIALSCHEME_3D_TTT_HPP
#define QUICC_SPATIALSCHEME_3D_TTT_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief ID of the Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme
    */
   class TTT: public ISpatialScheme
   {
      public:
         /**
          * @brief Constructor
          */
         explicit TTT(const VectorFormulation::Id formulation, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~TTT();

         /**
          * @brief Unique id
          */
         static const std::size_t sId;

         /**
          * @brief Create the scheme builder
          */
         virtual std::shared_ptr<IBuilder> createBuilder(ArrayI& dim, const bool needInterpretation) const override;

         /**
          * @brief Create transform for given dimension
          */
         virtual std::shared_ptr<Transform::ITransform> createTransform(const Dimensions::Transform::Id id, std::shared_ptr<Transform::TransformSetup> spSetup) const override;

         /**
          * @brief Create index converter for given dimension
          */
         virtual std::shared_ptr<Parallel::IIndexConv> createIndexConv(const Dimensions::Transform::Id id) const override;

         /**
          * @brief Create spectral decomposition tools
          */
         virtual std::shared_ptr<Equations::Tools::ICoupling> createCouplingTools(const Equations::CouplingIndexType indexType) const override;

         /**
          * @brief Get variant forward transform data type with correct type ininitialized
          */
         virtual VariantTransformDataPointer fwdPtr(const Dimensions::Transform::Id id) const override;

         /**
          * @brief Get variant forward transform data type with correct type ininitialized
          */
         virtual VariantTransformDataPointer bwdPtr(const Dimensions::Transform::Id id) const override;

         /**
          * @brief Create variant scalar variable
          */
         virtual ScalarVariable createSVar(std::shared_ptr<Resolution> spRes) const override;

         /**
          * @brief Create variant vector variable
          */
         virtual VectorVariable createVVar(std::shared_ptr<Resolution> spRes) const override;

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

}
}

#endif // QUICC_SPATIALSCHEME_3D_TTT_HPP
