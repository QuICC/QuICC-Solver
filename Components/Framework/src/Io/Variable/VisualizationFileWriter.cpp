/**
 * @file VisualizationFileWriter.cpp
 * @brief Source of the implementation of the HDF5 visualisation file writer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/VisualizationFile.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   VisualizationFileWriter::VisualizationFileWriter(std::string type)
      : IVariableHdf5NWriter(Tags::VisualizationFile::BASENAME, Tags::VisualizationFile::EXTENSION, Tags::VisualizationFile::HEADER, type, Tags::VisualizationFile::VERSION, Dimensions::Space::PHYSICAL, true)
   {
   }

   VisualizationFileWriter::~VisualizationFileWriter()
   {
   }

   void VisualizationFileWriter::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the Physical parameters
      this->writePhysical();

      // Write the Physical parameters
      this->writeMesh();

      // Write the truncation information
      this->writeTruncation();

      // Write the run information
      this->writeRun();

      // Write all the scalars
      VisualizationFileWriter::scalar_iterator_range sRange = this->scalarRange();
      VisualizationFileWriter::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, sit->second))
         {
            std::visit([&](auto&& p){this->writePhysicalScalar(PhysicalNames::Coordinator::tag(sit->first), p->dom(0).phys());}, sit->second);
         }

         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, sit->second))
         {
            std::visit([&](auto&& p){this->writePhysicalVector(PhysicalNames::Coordinator::tag(sit->first)+"_grad", p->dom(0).grad().data());}, sit->second);
         }

         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad2());}, sit->second))
         {
            std::visit([&](auto&& p){this->writePhysicalTensor(PhysicalNames::Coordinator::tag(sit->first)+"_grad2", p->dom(0).grad2().data());}, sit->second);
         }
      }

      // Write all the vectors
      VisualizationFileWriter::vector_iterator_range vRange = this->vectorRange();
      VisualizationFileWriter::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasPhys());}, vit->second))
         {
            std::visit([&](auto&& p){this->writePhysicalVector(PhysicalNames::Coordinator::tag(vit->first), p->dom(0).phys().data());}, vit->second);
         }

         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasGrad());}, vit->second))
         {
            std::vector<FieldComponents::Spectral::Id> fId;
            fId.push_back(this->res().sim().ss().spectral().ONE());
            fId.push_back(this->res().sim().ss().spectral().TWO());
            fId.push_back(this->res().sim().ss().spectral().THREE());
            for(unsigned int i = 0; i < fId.size(); i ++)
            {
               if(fId.at(i) != FieldComponents::Spectral::NOTUSED)
               {
                  std::visit([&](auto&& p){this->writePhysicalVector(PhysicalNames::Coordinator::tag(vit->first) + "_" + Tools::IdToHuman::toTag(fId.at(i)) + "_grad", p->dom(0).grad(fId.at(i)).data());}, vit->second);
               }
            }
         }

         if(std::visit([&](auto&& p)->bool{return (p->dom(0).hasCurl());}, vit->second))
         {
            std::visit([&](auto&& p){this->writePhysicalVector(PhysicalNames::Coordinator::tag(vit->first)+"_curl", p->dom(0).curl().data());}, vit->second);
         }
      }

      // Close file
      this->postWrite();
   }

   void VisualizationFileWriter::writeMesh()
   {
      // Create the Physical parameters group
      hid_t group = H5Gcreate(this->file(), Tags::VisualizationFile::MESH.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      std::vector<FieldComponents::Physical::Id> gridId;
      gridId.push_back(this->res().sim().ss().physical().ONE());
      gridId.push_back(this->res().sim().ss().physical().TWO());
      gridId.push_back(this->res().sim().ss().physical().THREE());
      for(int i = 0; i < this->res().cpu()->nDim(); i ++)
      {
         if(gridId.at(i) != FieldComponents::Physical::NOTUSED)
         {
            this->writeArray(group, Tags::VisualizationFile::GRID+"_"+Tools::IdToHuman::toTag(gridId.at(i)), this->mMesh.at(i));
         }
      }

      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalScalar(const std::string& name, const Framework::Selector::PhysicalScalarField& scalar)
   {
      // Create the scalar group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tuple<int,int, const Framework::Selector::PhysicalScalarField::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(scalar);

      // Write the scalar field
      this->writeRegularField(group, name, fieldInfo);

      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalVector(const std::string& name, const std::map<FieldComponents::Physical::Id,Framework::Selector::PhysicalScalarField>& vector, const std::string& joint)
   {
      // Create the Field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tuple<int,int, const Framework::Selector::PhysicalScalarField::PointType *> > fieldInfo;

      std::map<FieldComponents::Physical::Id,Framework::Selector::PhysicalScalarField>::const_iterator it;
      for(it = vector.begin(); it != vector.end(); ++it)
      {
         // create component field information
         fieldInfo = Datatypes::FieldTools::createInfo(it->second);

         // Write the vector field
         this->writeRegularField(group, name + joint + Tools::IdToHuman::toTag(it->first), fieldInfo);
      }

      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalTensor(const std::string& name, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Framework::Selector::PhysicalScalarField>& tensor, const std::string& joint)
   {
      // Create the Field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tuple<int,int, const Framework::Selector::PhysicalScalarField::PointType *> > fieldInfo;

      std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Framework::Selector::PhysicalScalarField>::const_iterator it;
      for(it = tensor.begin(); it != tensor.end(); ++it)
      {
         // create component field information
         fieldInfo = Datatypes::FieldTools::createInfo(it->second);

         // Write the tensor field
         this->writeRegularField(group, name + joint + Tools::IdToHuman::toTag(it->first.first) + Tools::IdToHuman::toTag(it->first.second), fieldInfo);
      }

      // close group
      H5Gclose(group);
   }

}
}
}
