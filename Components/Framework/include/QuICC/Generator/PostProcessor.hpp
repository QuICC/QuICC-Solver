/**
 * @file PostProcessor.hpp
 * @brief High level implementation of a general post processor
 */

#ifndef QUICC_POSTPROCESSOR_HPP
#define QUICC_POSTPROCESSOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Simulation/SimulationBase.hpp"

namespace QuICC {

   /**
    * @brief High level implementation of a post processor state generator
    */
   class PostProcessor: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         PostProcessor() = default;

         /**
          * @brief Simple empty destructor
          */
         virtual ~PostProcessor() = default;

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
          * @brief Allow for implementation specific output tuning
          */
         virtual void tuneOutput();

         /**
          * @brief Write the output if required
          */
         void writeOutput();

         /**
          * @brief Allow for additional operators on the initial state input file
          *
          * @param spInitFile Shared initial state file
          */
         virtual void tuneInitialState(Io::Variable::SharedStateFileReader spInitFile);
   };

   /// Typedef for a shared pointer of a PostProcessor
   typedef std::shared_ptr<PostProcessor> SharedStateGenerator;

}

#endif // QUICC_POSTPROCESSOR_HPP
