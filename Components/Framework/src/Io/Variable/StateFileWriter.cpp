/**
 * @file StateFileWriter.cpp
 * @brief Source of the implementation of the HDF5 state file writer
 */

// System includes
//

// Project includes
//
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/Tags/StateFile.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   StateFileWriter::StateFileWriter(std::string type, const bool isRegular)
      : IVariableHdf5NWriter(Tags::StateFile::BASENAME, Tags::StateFile::EXTENSION, Tags::StateFile::HEADER, type, Tags::StateFile::VERSION, Dimensions::Space::SPECTRAL, isRegular)
   {
   }

   StateFileWriter::StateFileWriter(std::string name, std::string type, const bool isRegular)
      : IVariableHdf5NWriter(Tags::StateFile::BASENAME + name, Tags::StateFile::EXTENSION, Tags::StateFile::HEADER, type, Tags::StateFile::VERSION, Dimensions::Space::SPECTRAL, isRegular)
   {
   }

   void StateFileWriter::write()
   {
      Profiler::RegionFixture<2> fix("StateFileWriter::write");

      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the Physical parameters
      this->writePhysical();

      // Write the truncation information
      this->writeTruncation();

      // Write the run information
      this->writeRun();

      // Write all the scalars
      Profiler::RegionStart<3>("StateFileWriter::write-scalars");
      StateFileWriter::scalar_iterator_range sRange = this->scalarRange();
      StateFileWriter::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         std::visit(
               [&](auto&& p)
               {
                  this->writeSpectralScalar(PhysicalNames::Coordinator::tag(sit->first), p->dom(0).perturbation());
               }, sit->second);
      }
      Profiler::RegionStop<3>("StateFileWriter::write-scalars");

      // Write all the vectors
      Profiler::RegionStart<3>("StateFileWriter::write-vectors");
      StateFileWriter::vector_iterator_range vRange = this->vectorRange();
      StateFileWriter::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         std::visit(
               [&](auto&& p)
               {
                  this->writeSpectralVector(PhysicalNames::Coordinator::tag(vit->first), p->dom(0).perturbation().data());
               }, vit->second);
      }
      Profiler::RegionStop<3>("StateFileWriter::write-vectors");

      // Close file
      this->postWrite();
   }

}
}
}
