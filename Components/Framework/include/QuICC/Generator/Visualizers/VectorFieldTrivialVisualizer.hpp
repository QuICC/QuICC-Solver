/**
 * @file VectorFieldTrivialVisualizer.hpp
 * @brief Implementation of the basic vector field visualizer
 */

#ifndef QUICC_EQUATIONS_VECTORFIELDTRIVIALVISUALIZER_HPP
#define QUICC_EQUATIONS_VECTORFIELDTRIVIALVISUALIZER_HPP

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

   /**
    * @brief Implementation of the basic vector field visualizer
    */
   class VectorFieldTrivialVisualizer: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VectorFieldTrivialVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VectorFieldTrivialVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient, const bool viewCurl);

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
   };

   /// Typedef for a shared VectorFieldTrivialVisualizer
   typedef std::shared_ptr<VectorFieldTrivialVisualizer> SharedVectorFieldTrivialVisualizer;

}
}

#endif // QUICC_EQUATIONS_VECTORFIELDTRIVIALVISUALIZER_HPP
