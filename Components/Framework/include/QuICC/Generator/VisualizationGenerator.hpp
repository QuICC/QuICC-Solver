/**
 * @file VisualizationGenerator.hpp
 * @brief High level implementation of a general visualization file generator
 */

#ifndef QUICC_VISUALIZATIONGENERATOR_HPP
#define QUICC_VISUALIZATIONGENERATOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Simulation/SimulationBase.hpp"

namespace QuICC {

   /**
    * @brief High level implementation of a general visualization file generator
    */
   class VisualizationGenerator: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         VisualizationGenerator() = default;

         /**
          * @brief Simple empty destructor
          */
         virtual ~VisualizationGenerator() = default;

      protected:

      private:
         /**
          * @brief Do operations required just before starting the main loop
          */
         virtual void preRun();

         /**
          * @brief Do operations required during the main loop
          */
         virtual void mainRun();

         /**
          * @brief Do operations required just after finishing the main loop
          */
         virtual void postRun();

         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Allow for implementation specific output tuning
          */
         virtual void tuneOutput();

         /**
          * @brief Allow for additional operators on the initial state input file
          *
          * @param spInitFile Shared initial state file
          */
         virtual void tuneInitialState(Io::Variable::SharedStateFileReader spInitFile);
   };

   /// Typedef for a shared pointer of a VisualizationGenerator
   typedef std::shared_ptr<VisualizationGenerator> SharedVisualizationGenerator;

}

#endif // QUICC_VISUALIZATIONGENERATOR_HPP
