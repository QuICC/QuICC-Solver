/**
 * @file IVectorTimeAveragedEquation.hpp
 * @brief Base for the implementation of a vector time averaged equation
 */

#ifndef QUICC_EQUATIONS_IVECTORTIMEAVERAGEDEQUATION_HPP
#define QUICC_EQUATIONS_IVECTORTIMEAVERAGEDEQUATION_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base for the implementation of a vector time averaged equation
    */
   class IVectorTimeAveragedEquation: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorTimeAveragedEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorTimeAveragedEquation();

         /**
          * @brief Current simulation time update
          */
         virtual void setTime(const MHDFloat time, const bool finished) override;

         /**
          * @brief Set the smart pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Framework::Selector::VariantSharedVectorVariable spUnknown) override;

      protected:
         /**
          * @brief Update the stored value with the solver solution (real data)
          */
         virtual MHDVariant updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k) override;

      private:
         /**
          * @brief Multi staged timestep finished?
          */
         bool mTimeFinished;

         /**
          * @brief Timestep since last update
          */
         MHDFloat mTimestep;

         /**
          * @brief Time averaged field
          */
         Framework::Selector::VariantSharedSpectralVectorField mTimeAvg;
   };
}
}

#endif // QUICC_EQUATIONS_IVECTORTIMEAVERAGEDEQUATION_HPP
