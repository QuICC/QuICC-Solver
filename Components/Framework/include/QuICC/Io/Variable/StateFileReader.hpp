/** 
 * @file StateFileReader.hpp
 * @brief Implementation of HDF5 state file reader
 */

#ifndef QUICC_IO_VARIABLE_STATEFILEREADER_HPP
#define QUICC_IO_VARIABLE_STATEFILEREADER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/IVariableHdf5Reader.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of HDF5 state file reader
    */
   class StateFileReader: public IVariableHdf5Reader
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       Name of the file
          * @param type       Type of the file (typically scheme name)
          * @param isRegular  Is data regular?
          */
         StateFileReader(std::string name, std::string type, const bool isRegular);

         /**
          * @brief Destructor
          */
         virtual ~StateFileReader();

         /**
          * @brief Read State from file
          */
         virtual void read();

         /**
          * @brief Read resolution and data information from state file
          */
         void readSetup();

         /**
          * @brief Read Partial data from state file
          */
         template <typename TFilter> void readPartial();

         /**
          * @brief Get state file time
          */
         MHDFloat time() const;

         /**
          * @brief Get state file timestep
          */
         MHDFloat timestep() const;
         
      protected:
         /**
          * @brief Time read from file
          */
         MHDFloat mTime;

         /**
          * @brief Timestep read from file
          */
         MHDFloat mTimestep;

         /**
          * @brief Read run information to file
          */
         void readRun();

         /**
          * @brief Read scalar field values from file
          *
          * @param name    Name of the scalar field
          * @param rScalar Storage for the field
          * @param isRequired field is required
          */
         template <typename T> void readSpectralScalar(const std::string& name, typename Framework::Selector::ScalarField<T>& rScalar, const bool isRequired);

         /**
          * @brief Read vector field values from file
          *
          * @param name    Name of the vector field
          * @param rVector Storage for the field
          * @param isRequired field is required
          */
         template <typename T> void readSpectralVector(const std::string& name, std::map<FieldComponents::Spectral::Id,typename Framework::Selector::ScalarField<T> >& rVector, const bool isRequired);

         /**
          * @brief Read vector field component from file
          *
          * @param name    Name of the vector field
          * @param id      ID of the component
          * @param rComp   Storage for field component
          * @param isRequired field is required
          */
         template <typename T> void readSpectralComponent(const std::string& name, FieldComponents::Spectral::Id id, typename Framework::Selector::ScalarField<T>& rComp, const bool isRequired);

         /**
          * @brief Adapt data after reading it in if necessary
          *
          * @param rField   Storage for field
          */
         template <typename T> void adaptData(typename Framework::Selector::ScalarField<T>& rComp);

      private:
   };

   inline MHDFloat StateFileReader::time() const
   {
      return this->mTime;
   }

   inline MHDFloat StateFileReader::timestep() const
   {
      return this->mTimestep;
   }

   template <typename T> void StateFileReader::readSpectralScalar(const std::string& name, typename Framework::Selector::ScalarField<T>& rScalar, const bool isRequired)
   {
      // Open the scalar group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      if(group >= 0)
      {
         // Storage for the field information
         std::vector<std::tuple<int,int, typename Framework::Selector::ScalarField<T>::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(rScalar);

         // Check for data regularity
         if(this->mIsRegular)
         {
            this->readRegularField(group, name, fieldInfo);
         } else
         {
            this->readIrregularField(group, name, fieldInfo);
         }
         
         // close group
         H5Gclose(group);

         // Adapt data if necessary
      //   this->adaptData(rScalar);

      } else if(isRequired)
      {
         throw std::logic_error("Tried to open inexistant HDF5 group");
      }
   }

   template <typename T> void StateFileReader::readSpectralVector(const std::string& name, std::map<FieldComponents::Spectral::Id,typename Framework::Selector::ScalarField<T> >& rVector, const bool isRequired)
   {
      // Open the vector field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      if(group >= 0)
      {
         // Storage for the field information
         std::vector<std::tuple<int,int, typename Framework::Selector::ScalarField<T>::PointType *> > fieldInfo;

         // Check for data regularity
         if(this->mIsRegular)
         {
            for(auto it = rVector.begin(); it != rVector.end(); ++it)
            {
               // create component field information
               fieldInfo = Datatypes::FieldTools::createInfo(it->second);

               // Read component from file 
               this->readRegularField(group,name+"_"+Tools::IdToHuman::toTag(it->first), fieldInfo);

               // Adapt data if necessary
               this->adaptData(it->second);
            }
         } else
         {
            for(auto it = rVector.begin(); it != rVector.end(); ++it)
            {
               // create component field information
               fieldInfo = Datatypes::FieldTools::createInfo(it->second);

               // Read component from file 
               this->readIrregularField(group,name+"_"+Tools::IdToHuman::toTag(it->first), fieldInfo);

               // Adapt data if necessary
               this->adaptData(it->second);
            }
         }
         
         // close group
         H5Gclose(group);

      } else if(isRequired)
      {
         throw std::logic_error("Tried to open inexistant HDF5 group");
      }
   }

   template <typename T> void StateFileReader::readSpectralComponent(const std::string& name, FieldComponents::Spectral::Id id, typename Framework::Selector::ScalarField<T>& rComp, const bool isRequired)
   {
      // Open the vector field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      if(group >= 0)
      {
         // Storage for the field information
         std::vector<std::tuple<int,int, typename Framework::Selector::ScalarField<T>::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(rComp);

         // Check for data regularity
         if(this->mIsRegular)
         {
            // Read the field component
            this->readRegularField(group, name+"_"+Tools::IdToHuman::toTag(id), fieldInfo);
         } else
         {
            // Read the field component
            this->readIrregularField(group, name+"_"+Tools::IdToHuman::toTag(id), fieldInfo);
         }

         // close group
         H5Gclose(group);

         // Adapt data if necessary
         this->adaptData(rComp);

      } else if(isRequired)
      {
         throw std::logic_error("Tried to open inexistant HDF5 group");
      }

   }

   template <typename T> void StateFileReader::adaptData(typename Framework::Selector::ScalarField<T>& rField)
   {
      if(this->res().sim().ss().has(SpatialScheme::Feature::FourierIndex23))
      {
         const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
         {
            if(tRes.idx<Dimensions::Data::DAT2D>(0, k) == 0 && tRes.idx<Dimensions::Data::DAT3D>(k) > this->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
            {
               rField.setProfile(Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)), 0, k);
            }
         }
      }
   }

   /// Typedef for a smart reference counting pointer of a HDF5 state file reader
   typedef std::shared_ptr<StateFileReader>   SharedStateFileReader;

}
}
}

#endif // QUICC_IO_VARIABLE_STATEFILEREADER_HPP
