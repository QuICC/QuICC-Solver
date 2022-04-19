/**
 * @file Coordinator.hpp
 * @brief Coordinator for the diagnostics computations
 */

#ifndef QUICC_DIAGNOSTICS_COORDINATOR_HPP
#define QUICC_DIAGNOSTICS_COORDINATOR_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Diagnostics/ICflWrapper.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Coordinator for the diagnostics computations
    */
   class Coordinator
   {
      public:
         /**
          * @brief Constructor
          */
         Coordinator();

         /**
          * @brief Destructor
          */
         virtual ~Coordinator();

         /**
          * @brief Initialise the coordinator
          *
          * @param mesh    Vector of grid values
          * @param scalars Map of shared scalar variables
          * @param vectors Map of shared vector variables
          * @param tstep   Timestep information
          */
         void init(const std::vector<Array>& mesh, const std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>&  scalars, const std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>&  vectors, const Array& tstep, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Compute the initial CFL condition
          */
         void initialCfl();

         /**
          * @brief Compute the current CFL condition
          */
         void updateCfl();

         /**
          * @brief Get CFL condition
          */
         const Matrix& cfl() const;

         /**
          * @brief Get max error goal
          */
         MHDFloat maxError() const;

         /**
          * @brief Get start time
          */
         MHDFloat startTime() const;

         /**
          * @brief Get start timestep
          */
         MHDFloat startTimestep() const;

         /**
          * @brief Use time and timestep from state file
          */
         void useStateTime(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Synchronize diagnostics among CPUs
          */
         void synchronize();

      protected:

      private:
         /**
          * @brief Special CFL location value for max step condition
          */
         const MHDFloat MAXSTEP_LOCATION;

         /**
          * @brief Special CFL location value for min step condition
          */
         const MHDFloat MINSTEP_LOCATION;

         /**
          * @brief Special CFL location value for fixed step condition
          */
         const MHDFloat FIXEDSTEP_LOCATION;

         /**
          * @brief Maximum timestep
          */
         const MHDFloat mcMaxStep;

         /**
          * @brief Minimum timestep
          */
         const MHDFloat mcMinStep;

         /**
          * @brief Fixed timestep
          */
         MHDFloat mFixedStep;

         /**
          * @brief Max error goal
          */
         MHDFloat mMaxError;

         /**
          * @brief Current CFL condition
          */
         Matrix mCfl;

         /**
          * @brief Start simulation time
          */
         MHDFloat mStartTime;

         /**
          * @brief Start simulation timestep
          */
         MHDFloat mStartTimestep;

         /**
          * @brief Shared pointer to a CFL condition wrapper
          */
         SharedICflWrapper  mspCflWrapper;
   };

}
}

#endif // QUICC_DIAGNOSTICS_COORDINATOR_HPP
