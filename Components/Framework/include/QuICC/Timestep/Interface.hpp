/**
 * @file Interface.hpp
 * @brief Interface for timestepper implementation
 */

#ifndef QUICC_TIMESTEP_INTERFACE_HPP
#define QUICC_TIMESTEP_INTERFACE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/BasicTypes.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Io/Ascii/CflWriter.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class Interface
   {
      public:
         /// Typedef for a shared scalar equation iterator
         typedef typename std::vector<Equations::SharedIScalarEquation>::iterator   ScalarEquation_iterator;

         /// Typedef for a shared vector equation iterator
         typedef typename std::vector<Equations::SharedIVectorEquation>::iterator   VectorEquation_iterator;

         /// Typedef for a shared scalar equation range
         typedef std::pair<ScalarEquation_iterator, ScalarEquation_iterator>  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef std::pair<VectorEquation_iterator, VectorEquation_iterator>  VectorEquation_range;

         /// Typedef for a shared scalar variable map
         typedef std::map<std::size_t, Framework::Selector::VariantSharedScalarVariable>  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef std::map<std::size_t, Framework::Selector::VariantSharedVectorVariable>  VectorVariable_map;

         /**
          * @brief Constructor
          *
          * @param time    Initial time value
          * @param cfl     Initial CFL timestep
          * @param error   Max error allowed during timestep
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         Interface(const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Destructor
          */
         virtual ~Interface();

         /**
          * @brief Timestep is finished?
          */
         virtual bool finishedStep() const = 0;

         /**
          * @brief Set solve time
          */
         virtual void setSolveTime(const std::size_t timeId) = 0;

         /**
          * @brief Update equation explicit linear input to solver
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         virtual void getExplicitInput(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar) = 0;

         /**
          * @brief Tune adaptive timestepper
          *
          * \mhdBug Not fully implemented
          */
         virtual void tuneAdaptive(const MHDFloat time) = 0;

         /**
          * @brief Adapt the timestep used
          *
          * @param cfl     CFL conditions
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         virtual void adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq) = 0;

         /**
          * @brief Update control status
          */
         virtual void update();

         /**
          * @brief Compute (partial) forward step
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         virtual void stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar) = 0;

         /**
          * @brief Print timestepper information to stream
          *
          * @param stream  Output stream
          */
         virtual void printInfo(std::ostream& stream) = 0;

         /**
          * @brief Get current simulation time
          */
         virtual MHDFloat time() const;

         /**
          * @brief Get current simulation timestep
          */
         virtual MHDFloat timestep() const;

      protected:
         /**
          * @brief Minimum number of constant timestep before step size increase
          */
         const MHDFloat mcMinCnst;

         /**
          * @brief Maximum timestep jump per step (See Soederlind)
          */
         const MHDFloat mcMaxJump;

         /**
          * @brief No update window for timestep increase
          */
         const MHDFloat mcUpWindow;

         /**
          * @brief Minimal timestep allowed before simulation abort
          */
         const MHDFloat mcMinDt;

         /**
          * @brief Maximum timestep allowed
          */
         const MHDFloat mcMaxDt;

         /**
          * @brief Maximum error allowed
          */
         MHDFloat mMaxError;

         /**
          * @brief Previous timestep length
          */
         MHDFloat mOldDt;

         /**
          * @brief Timestep length
          */
         Matrix mDt;

         /**
          * @brief Current time
          */
         MHDFloat mTime;

         /**
          * @brief Current reference time
          */
         MHDFloat mRefTime;

         /**
          * @brief Constant timestep steps
          */
         MHDFloat mCnstSteps;

         /**
          * @brief Constant timestep steps
          */
         MHDFloat mStepTime;

         /**
          * @brief Shared CFL writer
          */
         Io::Ascii::SharedCflWriter   mspIo;

      private:
   };

} // Timestep
} // QuICC

#endif // QUICC_TIMESTEP_INTERFACE_HPP
