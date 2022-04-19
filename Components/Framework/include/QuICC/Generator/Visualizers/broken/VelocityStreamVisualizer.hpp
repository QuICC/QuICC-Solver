/**
 * @file VelocityStreamVisualizer.hpp
 * @brief Implementation of the streamfunction + axial velocity to velocity field visualizer
 */

#ifndef QUICC_EQUATIONS_VELOCITYSTREAMVISUALIZER_HPP
#define QUICC_EQUATIONS_VELOCITYSTREAMVISUALIZER_HPP

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
    * @brief Implementation of the streamfunction + axial velocity to velocity field visualizer
    */
   class VelocityStreamVisualizer: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VelocityStreamVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityStreamVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient);

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
   };

   /// Typedef for a shared VelocityStreamVisualizer
   typedef std::shared_ptr<VelocityStreamVisualizer> SharedVelocityStreamVisualizer;

}
}

#endif // QUICC_EQUATIONS_VELOCITYSTREAMVISUALIZER_HPP
