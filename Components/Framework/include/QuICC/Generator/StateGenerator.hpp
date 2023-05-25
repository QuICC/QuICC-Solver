/**
 * @file StateGenerator.hpp
 * @brief High level implementation of a general state generator
 */

#ifndef QUICC_STATEGENERATOR_HPP
#define QUICC_STATEGENERATOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Simulation/SimulationBase.hpp"

namespace QuICC {

   /**
    * @brief High level implementation of a general state generator
    */
   class StateGenerator: public SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         StateGenerator() = default;

         /**
          * @brief Simple empty destructor
          */
         virtual ~StateGenerator() = default;

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
   };

   /// Typedef for a shared pointer of a StateGenerator
   typedef std::shared_ptr<StateGenerator> SharedStateGenerator;

}

#endif // QUICC_STATEGENERATOR_HPP
