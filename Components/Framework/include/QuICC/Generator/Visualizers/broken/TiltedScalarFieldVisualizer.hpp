/**
 * @file TiltedScalarFieldVisualizer.hpp
 * @brief Implementation of the tilted scalar field visualizer
 */

#ifndef QUICC_EQUATIONS_TILTEDSCALARFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_TILTEDSCALARFIELDVISUALIZER_HPP

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
    * @brief Implementation of the tilted scalar field visualizer
    */
   class TiltedScalarFieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         TiltedScalarFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TiltedScalarFieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set the data field name
          */
         void setDataField(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient);

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
          * @brief Storage for output field flag
          */
         bool mViewField;

         /**
          * @brief Storage for output gradient flag
          */
         bool mViewGradient;

         /**
          * @brief Physical name of data field
          */
         std::size_t mDataField;
   };

   /// Typedef for a shared TiltedScalarFieldVisualizer
   typedef std::shared_ptr<TiltedScalarFieldVisualizer> SharedTiltedScalarFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_TILTEDSCALARFIELDVISUALIZER_HPP
