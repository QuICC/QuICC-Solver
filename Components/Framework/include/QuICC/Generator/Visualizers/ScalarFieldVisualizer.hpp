/**
 * @file ScalarFieldVisualizer.hpp
 * @brief Implementation of the basic scalar field visualizer
 */

#ifndef QUICC_EQUATIONS_SCALARFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_SCALARFIELDVISUALIZER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the basic scalar field visualizer
    */
   class ScalarFieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         ScalarFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ScalarFieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient, const bool viewGradient2 = false);

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
          * @brief Storage for output 2nd order gradient flag
          */
         bool mViewGradient2;
   };

   /// Typedef for a shared ScalarFieldVisualizer
   typedef std::shared_ptr<ScalarFieldVisualizer> SharedScalarFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_SCALARFIELDVISUALIZER_HPP
