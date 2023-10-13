/**
 * @file StateFileReader.cpp
 * @brief Source of the implementation of HDF5 state file reader
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/SpatialScheme/3D/TFF.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/StateFile.hpp"
#include "QuICC/Io/Variable/Tags/VariableHdf5.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   StateFileReader::StateFileReader(std::string name, std::string type, const bool isRegular)
      : IVariableHdf5Reader(Tags::StateFile::BASENAME + name, Tags::StateFile::EXTENSION, Tags::StateFile::HEADER, type, Tags::StateFile::VERSION, Dimensions::Space::SPECTRAL, isRegular), mTime(-1.0), mTimestep(-1.0)
   {
   }

   void StateFileReader::read()
   {
      Profiler::RegionFixture<2> fix("StateFileReader::read");

      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();

      // Read all the scalars
      Profiler::RegionStart<3>("StateFileReader::read-scalars");
      StateFileReader::scalar_iterator_range sRange = this->scalarRange();
      StateFileReader::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         // Make sure full field is zero
         std::visit(
               [&](auto&& p)
               {
                  p->setZeros();
               }, sit->second);

         // Read field values
         std::visit(
               [&](auto&& p)
               {
                  this->readSpectralScalar(PhysicalNames::Coordinator::tag(sit->first), p->rDom(0).rPerturbation(), this->isRequired(sit->first));
               }, sit->second);
      }
      Profiler::RegionStop<3>("StateFileReader::read-scalars");

      // Read all the vectors
      Profiler::RegionStart<3>("StateFileReader::read-vectors");
      StateFileReader::vector_iterator_range vRange = this->vectorRange();
      StateFileReader::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         // Make sure full field is zero
         std::visit(
               [&](auto&& p)
               {
                  p->setZeros();
               }, vit->second);

         // Read field values
         std::visit(
               [&](auto&& p)
               {
                  this->readSpectralVector(PhysicalNames::Coordinator::tag(vit->first), p->rDom(0).rPerturbation().rData(), this->isRequired(vit->first));
               }, vit->second);
      }
      Profiler::RegionStop<3>("StateFileReader::read-vectors");
   }

   void StateFileReader::readSetup()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();
   }

   void StateFileReader::readRun()
   {
      // Open the run paramters group
      hid_t group = H5Gopen(this->file(), Tags::VariableHdf5::RUN.c_str(), H5P_DEFAULT);

      // Read the reached simulation time from file
      this->readScalar(group, Tags::VariableHdf5::RUNTIME, this->mTime);

      // Read the last used timestep from file
      this->readScalar(group, Tags::VariableHdf5::RUNSTEP, this->mTimestep);

      // close group
      H5Gclose(group);
   }

}
}
}
