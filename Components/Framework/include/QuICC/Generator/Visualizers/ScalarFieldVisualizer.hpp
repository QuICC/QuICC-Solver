/**
 * @file ScalarFieldVisualizer.hpp
 * @brief Implementation of the basic scalar field visualizer
 */

#ifndef QUICC_EQUATIONS_SCALARFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_SCALARFIELDVISUALIZER_HPP

// System includes
//
#include <memory>

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
         virtual ~ScalarFieldVisualizer() = default;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient, const bool viewGradient2 = false);

         /**
          * @brief Set backward path ID
          *
          * @param pathId Path ID
          */
         void setBackwardPath(const std::size_t id);

         /**
          * @brief Set forward path ID
          *
          * @param pathId Path ID
          */
         void setForwardPath(const std::size_t id);

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<Transform::TransformPath> backwardPaths() override;

      protected:
         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents() override;

         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling() override;

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

         /**
          * @brief Backward Transform path ID
          */
         std::size_t mBwdPathId;

         /**
          * @brief Forward Transform path ID
          */
         std::size_t mFwdPathId;
   };

   /// Typedef for a shared ScalarFieldVisualizer
   typedef std::shared_ptr<ScalarFieldVisualizer> SharedScalarFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_SCALARFIELDVISUALIZER_HPP
