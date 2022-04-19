/**
 * @file NonlinearScalarFieldVisualizer.hpp
 * @brief Implementation of a nonlinear scalar field visualizer
 */

#ifndef QUICC_EQUATIONS_NONLINEARSCALARFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_NONLINEARSCALARFIELDVISUALIZER_HPP

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

   struct NonlinearScalarVisualizerIds
   {
      enum Id {
         CYLINDER_HEAT_ADVECTION = 0,
      };
   };

   /**
    * @brief Implementation of a nonlinear scalar field visualizer
    */
   class NonlinearScalarFieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         NonlinearScalarFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NonlinearScalarFieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setNonlinearType(const NonlinearScalarVisualizerIds::Id type);

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

      private:
         /**
          * @brief Type of nonlinear term to compute
          */
         NonlinearScalarVisualizerIds::Id mNonlinearType;
   };

   /// Typedef for a shared NonlinearScalarFieldVisualizer
   typedef std::shared_ptr<NonlinearScalarFieldVisualizer> SharedNonlinearScalarFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_NONLINEARSCALARFIELDVISUALIZER_HPP
