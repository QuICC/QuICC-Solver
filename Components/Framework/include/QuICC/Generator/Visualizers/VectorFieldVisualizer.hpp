/**
 * @file VectorFieldVisualizer.hpp
 * @brief Implementation of the basic vector field visualizer
 */

#ifndef QUICC_EQUATIONS_VECTORFIELDVISUALIZER_HPP
#define QUICC_EQUATIONS_VECTORFIELDVISUALIZER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the basic vector field visualizer
    */
   class VectorFieldVisualizer: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VectorFieldVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VectorFieldVisualizer() = default;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient, const bool viewCurl);

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
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling() override;

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents() override;

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
          * @brief Storage for output curl flag
          */
         bool mViewCurl;

         /**
          * @brief Backward Transform path ID
          */
         std::size_t mBwdPathId;

         /**
          * @brief Forward Transform path ID
          */
         std::size_t mFwdPathId;
   };

   /// Typedef for a shared VectorFieldVisualizer
   typedef std::shared_ptr<VectorFieldVisualizer> SharedVectorFieldVisualizer;

}
}

#endif // QUICC_EQUATIONS_VECTORFIELDVISUALIZER_HPP
