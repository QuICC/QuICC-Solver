/**
 * @file Coordinator.hpp
 * @brief Implementation of a general timestep coordinator structure
 */

#ifndef QUICC_TIMESTEP_COORDINATOR_HPP
#define QUICC_TIMESTEP_COORDINATOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/Timestep/Interface.hpp"

namespace QuICC {

// Forward declaration
namespace Pseudospectral {
   class Coordinator;
}

namespace Timestep {

   /**
    * @brief Implementation of general timestepper structure
    */
   class Coordinator
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
          */
         Coordinator();

         /**
          * @brief Destructor
          */
         ~Coordinator() = default;

         /**
          * @brief Initialise timestepper
          *
          * @param schemeId   ID of timestepping scheme
          * @param time       Initial time value
          * @param cfl        Initial CFL timestep
          * @param error      Max error allowed during timestep
          * @param scalEq     Shared scalar equations
          * @param vectEq     Shared vector equations
          * @param pseudo     Pseudospectral coordinator
          */
         void init(const std::size_t schemeId, const MHDFloat time, const Matrix& cfl, const MHDFloat maxError, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, Pseudospectral::Coordinator& pseudo);

         /**
          * @brief Tune adaptive timestepper
          *
          * \mhdBug Not fully implemented
          */
         void tuneAdaptive(const MHDFloat time);

         /**
          * @brief Adapt the timestep used
          *
          * @param cfl     CFL conditions
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         void adaptTimestep(const Matrix& cfl, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Update control status
          */
         void update();

         /**
          * @brief Compute (partial) forward step
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void stepForward(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

         /**
          * @brief Get current simulation time
          */
         MHDFloat time() const;

         /**
          * @brief Get current simulation timestep
          */
         MHDFloat timestep() const;

         /**
          * @brief Timestep is finished?
          */
         bool finishedStep() const;

         /**
          * @brief Set solve time
          */
         void setSolveTime(const std::size_t timeId);

         /**
          * @brief Update equation explicit linear input to solver
          *
          * @param scalEq Scalar equations
          * @param vectEq Vector equations
          */
         void getExplicitInput(const std::size_t opId, const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

         /**
          * @brief Print timestepper information to stream
          *
          * @param stream  Output stream
          */
         void printInfo(std::ostream& stream);

      protected:

      private:
         /**
          * @brief Pointer to implementation
          */
         std::shared_ptr<Interface> mpImpl;
   };

}
}

#endif // QUICC_TIMESTEP_COORDINATOR_HPP
