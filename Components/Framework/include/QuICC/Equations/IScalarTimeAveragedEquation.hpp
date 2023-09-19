/**
 * @file IScalarTimeAveragedEquation.hpp
 * @brief Base for the implementation of a time averaged scalar equation
 */

#ifndef QUICC_EQUATIONS_ISCALARTIMEAVERAGEDEQUATION_HPP
#define QUICC_EQUATIONS_ISCALARTIMEAVERAGEDEQUATION_HPP

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
    * @brief Base for the implementation of a time averaged scalar equation
    */
   class IScalarTimeAveragedEquation: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarTimeAveragedEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarTimeAveragedEquation();

         /**
          * @brief Current simulation time update
          */
         virtual void setTime(const MHDFloat time, const bool finished) override;

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Framework::Selector::VariantSharedScalarVariable spUnknown) override;

      protected:
         /**
          * @brief Update the stored value with the solver solution (read data)
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
          * @brief Storage of the previous values
          */
         Framework::Selector::VariantSharedScalarField mTimeAvg;
   };
}
}

#endif // QUICC_EQUATIONS_ISCALARTIMEAVERAGEDEQUATION_HPP
