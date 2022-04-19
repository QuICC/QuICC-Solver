/**
 * @file SimulationRunControl.hpp
 * @brief Implementation of a general simulation control structure
 */

#ifndef QUICC_SIMULATIONRUNCONTROL_HPP
#define QUICC_SIMULATIONRUNCONTROL_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Io/Control/ControlInterface.hpp"

namespace QuICC {

   /**
    * @brief Implementation of simulation control structure
    */
   class SimulationRunControl
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationRunControl();

         /**
          * @brief Destructor
          */
         virtual ~SimulationRunControl();

         /**
          * @brief Update control status with simulation information
          *
          * @param simTime Simulation time
          * @param simDt   Simulation timestep
          */
         void updateSimulation(const MHDFloat simTime, const MHDFloat simDt);

         /**
          * @brief Update control status with cluster information
          *
          * @param wallTime Wall time
          */
         void updateCluster(const MHDFloat wallTime);

         /**
          * @brief Should the simulation keep running?
          */
         std::size_t status() const;

         /**
          * @brief Update the status from control file
          */
         void checkFile();

         /***
          * @brief Set the maximum simulation time
          *
          * @param maxTime New maximum simulation time
          */
         void setMaxSimTime(const MHDFloat maxTime);

         /***
          * @brief Set the maximum wall time
          *
          * @param maxTime New maximum wall time
          */
         void setMaxWallTime(const MHDFloat maxTime);

         /**
          * @brief Print run information
          *
          * @param stream  Output stream
          */
         void printInfo(std::ostream& stream);

         /**
          * @brief Static method to close simulation cleanly through signal
          */
         static void handleSignal(int signum);

      protected:
         /**
          * @brief Status to be switched by signal
          */
         static std::size_t   SIGNAL_STATUS;

         /**
          * @brief Current runtime status
          */
         std::size_t mStatus;

         /**
          * @brief External control file
          */
         Io::Control::ControlInterface  mCtrlFile;

         /**
          * @brief Number of timesteps done
          */
         int   mSteps;

         /**
          * @brief Maximum simulation time
          */
         MHDFloat mMaxSimTime;

         /**
          * @brief Maximum wall time
          */
         MHDFloat mMaxWallTime;

      private:
         /**
          * @brief Initialise the signal handler
          */
         void initSignalHandler();
   };
}

#endif // QUICC_SIMULATIONRUNCONTROL_HPP
