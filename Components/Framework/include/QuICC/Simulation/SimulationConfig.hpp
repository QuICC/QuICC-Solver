/**
 * @file SimulationConfig.hpp
 * @brief Interface to simulation configuration file
 */

#ifndef QUICC_SIMULATIONCONFIG_HPP
#define QUICC_SIMULATIONCONFIG_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Io/Config/ConfigurationReader.hpp"

namespace QuICC {

   /**
    * @brief Interface to simulation configuration file
    */
   class SimulationConfig
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationConfig(Io::Config::SharedConfigurationReader spCfgFile);

         /**
          * @brief Destructor
          */
         ~SimulationConfig();

         /**
          * @brief Get the dimension read from the configuration file
          */
         ArrayI dimension() const;

         /**
          * @brief Get the transform implementation
          */
         ArrayI transformSetup() const;

         /**
          * @brief Get the box scale read from the configuration file
          */
         Array boxScale() const;

         /**
          * @brief Get the number of CPUs read from the configuration file
          */
         int nCpu() const;

         /**
          * @brief Get the parallel splitting algorithm read from the configuration file
          */
         Splitting::Algorithms::Id algorithm() const;

         /**
          * @brief Get the parallel communication grouping algorithm read from the configuration file
          */
         Splitting::Groupers::Id grouper() const;

         /**
          * @brief Get the map of physical values read from the configuration file
          */
         const std::map<std::string, MHDFloat>& physical() const;

         /**
          * @brief Update the map of physical values read from the configuration file
          */
         std::map<std::string, MHDFloat>& physical();

         /**
          * @brief Update the map of physical values read from the configuration file
          */
         std::map<std::string, MHDFloat>& rPhysical();

         /**
          * @brief Get the map of boundary conditions read from the configuration file
          */
         std::map<std::string, std::size_t> boundary() const;

         /**
          * @brief Get the map of boundary conditions read from the configuration file
          */
         const std::map<std::string, int>& model(const std::string tag) const;

         /**
          * @brief Get the run options read from the configuration file
          */
         Array run() const;

         /**
          * @brief Get the timestepping options read from the configuration file
          */
         Array timestepping() const;

         /**
          * @brief ID of the timestepper
          */
         std::size_t timestepper() const;

         /**
          * @brief ID of the boundary condition scheme
          */
         std::size_t bcScheme() const;

         /**
          * @brief Use split equation model setup
          */
         bool splitEquation() const;

      protected:

      private:
         /**
          * @brief Handle to the configuration file
          */
         Io::Config::SharedConfigurationReader    mspCfgFile;
   };

}

#endif // QUICC_SIMULATIONCONFIG_HPP
