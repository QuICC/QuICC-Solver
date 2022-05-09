/**
 * @file SimulationIoControl.hpp
 * @brief Implementation of the control instance for the simulation IO related operations
 */

#ifndef QUICC_SIMULATIONIOCONTROL_HPP
#define QUICC_SIMULATIONIOCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Io/Ascii/StdOutPipe.hpp"
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"
#include "QuICC/Io/Variable/IVariableHdf5NWriter.hpp"
#include "QuICC/Io/Stats/IStatisticsAsciiWriter.hpp"
#include "QuICC/Io/Config/ConfigurationReader.hpp"
#include "QuICC/Simulation/SimulationConfig.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the control instance for the simulation IO related operations
    */
   class SimulationIoControl
   {
      public:
         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<Io::Variable::SharedIVariableAsciiWriter>::iterator ascii_iterator;

         /// Typedef for an iterator over all the HDF5 writers
         typedef std::vector<Io::Variable::SharedIVariableHdf5NWriter>::iterator hdf5_iterator;

         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<Io::Stats::SharedIStatisticsAsciiWriter>::iterator stats_iterator;

         /**
          * @brief Constructor
          */
         SimulationIoControl();

         /**
          * @brief Destructor
          */
         ~SimulationIoControl();

         /**
          * @brief Set the configuration file
          */
         void setConfigurationFile(Io::Config::SharedConfigurationReader spCfgFile);

         /**
          * @brief Initialise the IO system
          */
         void init();

         /**
          * @brief Cleanup unused memory from IO system
          */
         void cleanup();

         /**
          * @brief Initialise the writer files
          */
         void initWriters();

         /**
          * @brief Finalise the writer files
          */
         void finalizeWriters();

         /**
          * @brief Add an ASCII output file
          *
          * @param spOutFile Shared ASCII writer
          */
         void addAsciiOutputFile(Io::Variable::SharedIVariableAsciiWriter spOutFile);

         /**
          * @brief Add a HDF5 output file
          *
          * @param spOutFile Shared HDF5 writer
          */
         void addHdf5OutputFile(Io::Variable::SharedIVariableHdf5NWriter  spOutFile);

         /**
          * @brief Add an statistics output file
          *
          * @param spOutFile Shared statistics writer
          */
         void addStatsOutputFile(Io::Stats::SharedIStatisticsAsciiWriter spOutFile);

         /**
          * @brief Get interface to configuration file
          */
         const SimulationConfig& config() const;

         /**
          * @brief Get interface to configuration file
          */
         SimulationConfig& rConfig();

         /**
          * @brief Get begin iterator to ASCII files
          */
         ascii_iterator beginAscii();

         /**
          * @brief Get end iterator to ASCII files
          */
         ascii_iterator endAscii();

         /**
          * @brief Get begin iterator to HDF5 files
          */
         hdf5_iterator beginHdf5();

         /**
          * @brief Get end iterator to HDF5 files
          */
         hdf5_iterator endHdf5();

         /**
          * @brief Get begin iterator to statistics files
          */
         stats_iterator beginStats();

         /**
          * @brief Get end iterator to statistics files
          */
         stats_iterator endStats();

         /**
          * @brief Is is time to write ASCII file?
          */
         bool isAsciiTime() const;

         /**
          * @brief Is is time to write HDF5 file?
          */
         bool isHdf5Time() const;

         /**
          * @brief Is is time to update statistics calculations?
          */
         bool isStatsUpdateTime() const;

         /**
          * @brief Update IO status
          */
         void update();

         /**
          * @brief Write output files if conditions are met
          *
          * @param time    Current simulation time
          * @param time    Current simulation timestep
          */
         void writeFiles(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Write ASCII data
          *
          * @param time    Current simulation time
          * @param time    Current simulation timestep
          */
         void writeAscii(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Write HDF5 data
          *
          * @param time    Current simulation time
          * @param time    Current simulation timestep
          */
         void writeHdf5(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Prepare statistics data
          *
          * @param time    Current simulation time
          * @param time    Current simulation timestep
          */
         void prepareStats(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Write statistics data
          */
         void writeStats();

         /**
          * @brief Active stats calculation
          */
         void activateStats();

         /**
          * @brief Disable stats calculation
          */
         void disableStats();

      protected:

      private:
         /**
          * @brief Is is time to write statistics file?
          */
         bool isStatsTime() const;

         /**
          * @brief Vector of ASCII output files
          */
         std::vector<Io::Variable::SharedIVariableAsciiWriter> mAsciiWriters;

         /**
          * @brief Vector of HDF5 output files
          */
         std::vector<Io::Variable::SharedIVariableHdf5NWriter> mHdf5Writers;

         /**
          * @brief Vector of statistics output files
          */
         std::vector<Io::Stats::SharedIStatisticsAsciiWriter> mStatsWriters;

         /**
          * @brief Handle to StdMessage buffer
          */
         Io::Ascii::SharedStdOutPipe   mspStdOut;

         /**
          * @brief Handle to the configuration file
          */
         Io::Config::SharedConfigurationReader    mspCfgFile;

         /**
          * @brief Handle to the configuration file
          */
         std::shared_ptr<SimulationConfig> mspCfg;

         /**
          * @brief Timestep counter
          */
         int mSteps;

         /**
          * @brief ASCII output rate
          */
         int mAsciiRate;

         /**
          * @brief HDF5 output rate
          */
         int mHdf5Rate;

         /**
          * @brief Statistics output rate
          */
         int mStatsRate;

         /**
          * @brief Statistics time average rate
          */
         int mStatsAvgRate;

         /**
          * @brief Activate stats calculation
          */
         bool mActiveStatsUpdate;

         /**
          * @brief Activate stats writing
          */
         bool mActiveStatsWrite;

         /**
          * @brief Initialise the configuration file
          */
         void initCfg();

         /**
          * @brief Initialise the StdMessage file
          */
         void initStdOut();
   };

}

#endif // QUICC_SIMULATIONIOCONTROL_HPP
