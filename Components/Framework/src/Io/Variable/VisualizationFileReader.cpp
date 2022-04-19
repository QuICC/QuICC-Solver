/**
 * @file VisualizationFileReader.cpp
 * @brief Source of the implementation of the HDF5 visualisation file reader
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/VisualizationFileReader.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/VisualizationFile.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   VisualizationFileReader::VisualizationFileReader(std::string name, std::string type)
      : IVariableHdf5Reader(Tags::VisualizationFile::BASENAME + name, Tags::VisualizationFile::EXTENSION, Tags::VisualizationFile::HEADER, type, Tags::VisualizationFile::VERSION, Dimensions::Space::PHYSICAL, true)
   {
   }

   VisualizationFileReader::~VisualizationFileReader()
   {
   }

   void VisualizationFileReader::read()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read all the scalars
      VisualizationFileReader::scalar_iterator_range sRange = this->scalarRange();
      VisualizationFileReader::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         // Make sure full field is zero
         std::visit([&](auto&& p){p->setZeros();}, sit->second);

         // Read field values
         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, sit->second))
         {
            std::visit([&](auto&& p){this->readPhysicalScalar(PhysicalNames::Coordinator::tag(sit->first), p->rDom(0).rPhys(), this->isRequired(sit->first));}, sit->second);
         }

         // Read gradient values
         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, sit->second))
         {
            std::visit([&](auto&& p){this->readPhysicalVector(PhysicalNames::Coordinator::tag(sit->first)+"_grad", p->rDom(0).rGrad().rData(), this->isRequired(sit->first));}, sit->second);
         }

         // Read second gradient tensor
         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, sit->second))
         {
            std::visit([&](auto&& p){this->readPhysicalTensor(PhysicalNames::Coordinator::tag(sit->first)+"_grad2", p->rDom(0).rGrad2().rData(), this->isRequired(sit->first));}, sit->second);
         }
      }

      // Read all the vectors
      VisualizationFileReader::vector_iterator_range vRange = this->vectorRange();
      VisualizationFileReader::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         // Make sure full field is zero
         std::visit([&](auto&& p){p->setZeros();}, vit->second);

         // Read field values
         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, vit->second))
         {
            std::visit([&](auto&& p){this->readPhysicalVector(PhysicalNames::Coordinator::tag(vit->first), p->rDom(0).rPhys().rData(), this->isRequired(vit->first));}, vit->second);
         }

         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasCurl());}, vit->second))
         {
            std::visit([&](auto&& p){this->readPhysicalVector(PhysicalNames::Coordinator::tag(vit->first)+"_curl", p->rDom(0).rCurl().rData(), this->isRequired(vit->first));}, vit->second);
         }
      }
   }

   void VisualizationFileReader::readSetup()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();
   }

   void VisualizationFileReader::readMesh()
   {
   }

   void VisualizationFileReader::readPhysicalScalar(const std::string& name, Framework::Selector::PhysicalScalarField& rScalar, const bool isRequired)
   {
      // Open the scalar group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      if(group >= 0)
      {
         // Storage for the field information
         std::vector<std::tuple<int,int, Framework::Selector::PhysicalScalarField::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(rScalar);

         this->readRegularField(group, name, fieldInfo);

         // close group
         H5Gclose(group);

         // Adapt data if necessary
         this->adaptData(rScalar);

      } else if(isRequired)
      {
         throw std::logic_error("Tried to open inexistant HDF5 group");
      }
   }

   void VisualizationFileReader::readPhysicalVector(const std::string& name, std::map<FieldComponents::Physical::Id,Framework::Selector::PhysicalScalarField>& rVector, const bool isRequired)
   {
      // Open the vector field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      if(group >= 0)
      {
         // Storage for the field information
         std::vector<std::tuple<int,int, Framework::Selector::PhysicalScalarField::PointType *> > fieldInfo;

         // Check for data regularity
         for(auto it = rVector.begin(); it != rVector.end(); ++it)
         {
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(it->second);

            // Read component from file
            this->readRegularField(group,name+"_"+Tools::IdToHuman::toTag(it->first), fieldInfo);

            // Adapt data if necessary
            this->adaptData(it->second);
         }

         // close group
         H5Gclose(group);

      } else if(isRequired)
      {
         throw std::logic_error("Tried to open inexistant HDF5 group");
      }
   }

   void VisualizationFileReader::readPhysicalComponent(const std::string& name, FieldComponents::Physical::Id id, Framework::Selector::PhysicalScalarField& rComp, const bool isRequired)
   {
      // Open the vector field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      if(group >= 0)
      {
         // Storage for the field information
         std::vector<std::tuple<int,int, Framework::Selector::PhysicalScalarField::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(rComp);

         // Read the field component
         this->readRegularField(group, name+"_"+Tools::IdToHuman::toTag(id), fieldInfo);

         // close group
         H5Gclose(group);

         // Adapt data if necessary
         this->adaptData(rComp);

      } else if(isRequired)
      {
         throw std::logic_error("Tried to open inexistant HDF5 group");
      }

   }

   void VisualizationFileReader::readPhysicalTensor(const std::string&, std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Framework::Selector::PhysicalScalarField>&, const bool)
   {
   }

   void VisualizationFileReader::adaptData(Framework::Selector::PhysicalScalarField&)
   {
   }

}
}
}
