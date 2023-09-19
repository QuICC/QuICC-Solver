/**
 * @file VorticityStreamVisualizer.hpp
 * @brief Implementation of the streamfunction to vorticity visualizer
 */

#ifndef QUICC_EQUATIONS_VORTICITCYSTREAMVISUALIZER_HPP
#define QUICC_EQUATIONS_VORTICITCYSTREAMVISUALIZER_HPP

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
    * @brief Implementation of the streamfunction to vorticity visualizer
    */
   class VorticityStreamVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VorticityStreamVisualizer(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VorticityStreamVisualizer();

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
          * @brief Set the equation coupling information
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
   };

   /// Typedef for a shared VorticityStreamVisualizer
   typedef std::shared_ptr<VorticityStreamVisualizer> SharedVorticityStreamVisualizer;

}
}

#endif // QUICC_EQUATIONS_VORTICITCYSTREAMVISUALIZER_HPP
