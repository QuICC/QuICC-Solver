/**
 * @file IConfigurationFile.hpp
 * @brief Implementation of the base for a configuration file
 */

#ifndef QUICC_IO_CONFIG_ICONFIGURATIONFILE_HPP
#define QUICC_IO_CONFIG_ICONFIGURATIONFILE_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/MpiTypes.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Config/IConfigurationNode.hpp"
#include "QuICC/Io/Config/Framework/Block.hpp"
#include "QuICC/Io/Config/Simulation/Block.hpp"
#include "QuICC/Io/Config/Setup/Block.hpp"
#include "QuICC/Io/Config/Model/Block.hpp"

namespace QuICC {

namespace Io {

namespace Config {

   /**
    * @brief Implementation of the base for a configuration file
    */
   template <typename TBase> class IConfigurationFile: public TBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim  Dimensionality of simulation
          * @param name Name of the file
          * @param type Type of the simulation
          */
         IConfigurationFile(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& name, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~IConfigurationFile();

         /**
          * @brief Output run information
          */
         void printInfo() const;

         /**
          * @brief Get framework block
          */
         std::shared_ptr<const Framework::Block> spFramework() const;

         /**
          * @brief Get simulation block
          */
         std::shared_ptr<const Simulation::Block> spSimulation() const;

         /**
          * @brief Get simulation block
          */
         std::shared_ptr<Simulation::Block> rspSimulation();

         /**
          * @brief Get setup block
          */
         std::shared_ptr<const Setup::Block> spSetup() const;

         /**
          * @brief Get model block
          */
         const Model::Block& model() const;

         /**
          * @brief Get model block
          */
         std::shared_ptr<Model::Block> rspModel();

      protected:
         /**
          * @brief Framework part of configuration file
          */
         std::shared_ptr<Framework::Block>  mspFramework;

         /**
          * @brief Framework part of configuration file
          */
         std::shared_ptr<Simulation::Block>  mspSimulation;

         /**
          * @brief Setup block
          */
         std::shared_ptr<Setup::Block>  mspSetup;

         /**
          * @brief Model block
          */
         std::shared_ptr<Model::Block>  mspModel;

         /**
          * @brief Spread parameters over parallel simulation
          */
         void spreadParameters();

      private:
         /**
          * @brief EXTENSION of the configuration file
          */
         static const std::string EXTENSION;

         /**
          * @brief HEADER of the configuration file
          */
         static const std::string HEADER;

         /**
          * @brief VERSION of the configuration file
          */
         static const std::string VERSION;

         /**
          * @brief XML ROOT of the configuration file
          */
         static const std::string ROOT;
   };

   template <typename TBase> inline std::shared_ptr<const Framework::Block> IConfigurationFile<TBase>::spFramework() const
   {
      return this->mspFramework;
   }

   template <typename TBase> inline std::shared_ptr<const Simulation::Block> IConfigurationFile<TBase>::spSimulation() const
   {
      return this->mspSimulation;
   }

   template <typename TBase> inline std::shared_ptr<Simulation::Block> IConfigurationFile<TBase>::rspSimulation()
   {
      return this->mspSimulation;
   }

   template <typename TBase> inline std::shared_ptr<const Setup::Block> IConfigurationFile<TBase>::spSetup() const
   {
      return this->mspSetup;
   }

   template <typename TBase> inline const Model::Block& IConfigurationFile<TBase>::model() const
   {
      return *this->mspModel;
   }

   template <typename TBase> inline std::shared_ptr<Model::Block> IConfigurationFile<TBase>::rspModel()
   {
      return this->mspModel;
   }

   template <typename TBase> const std::string IConfigurationFile<TBase>::EXTENSION = ".cfg";

   template <typename TBase> const std::string IConfigurationFile<TBase>::HEADER = "ConfigurationFile";

   template <typename TBase> const std::string IConfigurationFile<TBase>::VERSION = "1.0";

   template <typename TBase> const std::string IConfigurationFile<TBase>::ROOT = "quicc";

   template <typename TBase> IConfigurationFile<TBase>::IConfigurationFile(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& name, const std::string& type)
      : TBase(name, IConfigurationFile<TBase>::EXTENSION, IConfigurationFile<TBase>::HEADER, type, IConfigurationFile<TBase>::VERSION, IConfigurationFile<TBase>::ROOT)
   {
      // Initialize framework
      this->mspFramework = std::make_shared<Framework::Block>();
      this->mspFramework->init(dim, isPeriodicBox);

      // Initialize setup
      this->mspSetup = std::make_shared<Setup::Block>();
      this->mspSetup->init(dim);

      // Initialize simulation
      this->mspSimulation = std::make_shared<Simulation::Block>();
      this->mspSimulation->init();

      // Initialize model
      this->mspModel = std::make_shared<Model::Block>();
      this->mspModel->init();
   }

   template <typename TBase> IConfigurationFile<TBase>::~IConfigurationFile()
   {
   }

   template <typename TBase> void IConfigurationFile<TBase>::printInfo() const
   {
      // Check if the framework allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         // Create output header
         Tools::Formatter::printNewline(std::cout);
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printCentered(std::cout, "Configuration parameters", '*');
         Tools::Formatter::printLine(std::cout, '-');
         Tools::Formatter::printNewline(std::cout);

         this->mspFramework->printInfo();

         this->mspSimulation->printInfo();

         this->mspSetup->printInfo();

         this->mspModel->printInfo();
      }
   }

   template <typename TBase> void IConfigurationFile<TBase>::spreadParameters()
   {
      //
      // Start of MPI block
      //
      #ifdef QUICC_MPI

      // Create MPI compatible storage for the integer values
      std::vector<int> iData;
      // Create MPI compatible storage for the float values
      std::vector<MHDFloat> fData;

      //
      // Gather data
      //

      // Gather data from framework
      this->mspFramework->gatherParameters(iData, fData);

      // Gather data from simulation
      this->mspSimulation->gatherParameters(iData, fData);

      // Gather data from setup
      this->mspSetup->gatherParameters(iData, fData);

      // Gather data from model
      this->mspModel->gatherParameters(iData, fData);

      // Various MPI broadcast data
      int nBlocks = 2;
      MPI_Aint    displ[nBlocks];
      int         blocks[nBlocks];
      MPI_Datatype   types[nBlocks];
      MPI_Aint    element;

      // Create integer data part
      MPI_Get_address(&iData[0], &element);
      displ[0] = element;
      blocks[0] = iData.size();
      types[0] = MPI_INT;

      // Create float data part
      MPI_Get_address(&fData[0], &element);
      displ[1] = element;
      blocks[1] = fData.size();
      types[1] = Parallel::MpiTypes::type<MHDFloat>();

      // MPI data type for the combined integer and float data
      MPI_Datatype   cfgType;

      // Create MPI datatype
      MPI_Type_create_struct(nBlocks, blocks, displ, types, &cfgType);
      // Commit MPI datatype
      MPI_Type_commit(&cfgType);

      // Broadcast the information
      MPI_Bcast(MPI_BOTTOM, 1, cfgType, QuICCEnv().ioRank(), MPI_COMM_WORLD);
      QuICCEnv().synchronize();

      // Free the datatype
      MPI_Type_free(&cfgType);

      //
      // Distribute data
      //

      // Global integer index
      int iIdx = 0;
      // Global float index
      int fIdx = 0;

      // Scatter data to framework
      this->mspFramework->scatterParameters(iIdx, fIdx, iData, fData);

      // Scatter data to simulation
      this->mspSimulation->scatterParameters(iIdx, fIdx, iData, fData);

      // Scatter data to setup
      this->mspSetup->scatterParameters(iIdx, fIdx, iData, fData);

      // Scatter data to model
      this->mspModel->scatterParameters(iIdx, fIdx, iData, fData);

      //
      // End of MPI block
      //
      #endif // QUICC_MPI
   }

}
}
}

#endif // QUICC_IO_CONFIG_ICONFIGURATIONFILE_HPP
