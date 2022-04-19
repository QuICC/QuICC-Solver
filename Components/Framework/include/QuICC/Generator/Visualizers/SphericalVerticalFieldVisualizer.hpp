/**
 * @file SphericalVerticalFieldVisualizer.hpp
 * @brief Implementation of the spherical vertical component field visualizer
 */

#ifndef QUICC_EQUATIONS_SPHERICALVERTICALFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_SPHERICALVERTICALFIELDVISUALIZER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the spherical vertical component field visualizer
    */
   class SphericalVerticalFieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         SphericalVerticalFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~SphericalVerticalFieldVisualizer();

         /**
          * @brief Set the vertical field name, the base vector field name and requirements
          */
         void setIdentity(const std::size_t vertName, const std::size_t fieldName);

         /**
          * @brief Set which component of the ase vector field to use for vertical component (ie. field, gradient, curl)
          */
         void setFieldType(const FieldType::Id type);

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false);

         /**
          * @brief Copy nonlinear calculation back into physical field
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void useNonlinear(const Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Storage for output field flag
          */
         FieldType::Id mFieldType;

         /**
          * @brief Storage for output curl flag
          */
         std::size_t mFieldName;
   };

   /// Typedef for a shared SphericalVerticalFieldVisualizer
   typedef std::shared_ptr<SphericalVerticalFieldVisualizer> SharedSphericalVerticalFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_SPHERICALVERTICALFIELDVISUALIZER_HPP
