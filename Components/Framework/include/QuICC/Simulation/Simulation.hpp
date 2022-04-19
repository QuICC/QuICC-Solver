/**
 * @file Simulation.hpp
 * @brief High level implementation of a simulation
 */

#ifndef QUICC_SIMULATION_HPP
#define QUICC_SIMULATION_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Simulation/SimulationBase.hpp"

namespace QuICC {

   /**
    * @brief High level implementation of a simulation's execution steps.
    */
   class Simulation: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         Simulation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~Simulation();

      protected:

      private:
         /**
          * @brief Initialise the generator specific base components
          */
         virtual void initAdditionalBase();

         /**
          * @brief Do operations required just before starting the time integration
          */
         virtual void preRun();

         /**
          * @brief Do operations required during the main loop
          *
          * \callgraph
          */
         virtual void mainRun();

         /**
          * @brief Do operations required just after finishing the time integration
          *
          * \callgraph
          */
         virtual void postRun();

         /**
          * @brief Write the output if required
          */
         void writeOutput();
   };

   /// Typedef for a shared pointer of a Simulation
   typedef std::shared_ptr<Simulation> SharedSimulation;

}

#endif // QUICC_SIMULATION_HPP
