/**
 * @file NonlinearVectorFieldVisualizer.hpp
 * @brief Implementation of a nonlinear vector field visualizer
 */

#ifndef QUICC_EQUATIONS_NONLINEARVECTORFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_NONLINEARVECTORFIELDVISUALIZER_HPP

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
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

   struct NonlinearVectorVisualizerIds
   {
      enum Id {
         CYLINDER_TORPOL_ADVECTION = 0,
      };
   };

   /**
    * @brief Implementation of a nonlinear vector field visualizer
    */
   class NonlinearVectorFieldVisualizer: public IVectorEquation
   {
      public:

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         NonlinearVectorFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NonlinearVectorFieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setNonlinearType(const NonlinearVectorVisualizerIds::Id type);

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const;

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

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents();

      private:
         /**
          * @brief Type of nonlinear term to compute
          */
         NonlinearVectorVisualizerIds::Id mNonlinearType;
   };

   /// Typedef for a shared NonlinearVectorFieldVisualizer
   typedef std::shared_ptr<NonlinearVectorFieldVisualizer> SharedNonlinearVectorFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_NONLINEARVECTORFIELDVISUALIZER_HPP
