/**
 * @file IVariableHdf5NWriter.cpp
 * @brief Source of the implementation of the generic variable to HDF5 file writer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/IVariableHdf5NWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Hasher.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/NonDimensional/Coordinator.hpp"
#include "QuICC/Io/Variable/Tags/VariableHdf5.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   IVariableHdf5NWriter::IVariableHdf5NWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const bool isRegular)
      : IHdf5NWriter(name, ext, header, type, version), mTime(-1.0), mTimestep(-1.0), mIsRegular(isRegular), mSpaceId(id)
   {
   }

   IVariableHdf5NWriter::~IVariableHdf5NWriter()
   {
   }

   Dimensions::Space::Id IVariableHdf5NWriter::space() const
   {
      return this->mSpaceId;
   }

   void IVariableHdf5NWriter::setPhysical(const std::map<std::string,MHDFloat>& parameters, const std::map<std::string,int>& boundary)
   {
      // Convert parameters to NonDimensional numbers
      for(auto it = parameters.cbegin(); it != parameters.cend(); ++it)
      {
         size_t nd = Hasher::makeId(it->first);
         this->mPhysical.insert(std::make_pair(nd,NonDimensional::Coordinator::map().find(nd)->second->create(it->second)));
      }

      // Save boundary flags
      this->mBoundary = boundary;
   }

   const Resolution& IVariableHdf5NWriter::res() const
   {
      return *this->mspRes;
   }

   void IVariableHdf5NWriter::setMesh(const std::vector<Array>& mesh)
   {
      this->mMesh = mesh;
   }

   void IVariableHdf5NWriter::setSimTime(const MHDFloat time, const MHDFloat timestep)
   {
      this->mTime = time;

      this->mTimestep = timestep;
   }

   void IVariableHdf5NWriter::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;

      // Set dataset dimensions
      this->setDatasetSize();

      // Set dataset offsets
      this->setDatasetOffsets();

      // This set collective IO operations
      this->setCollIo();
   }

   void IVariableHdf5NWriter::expect(const std::size_t id)
   {
      this->mExpected.insert(id);
   }

   bool IVariableHdf5NWriter::isFull() const
   {
      bool status = true;

      // Check that all expected scalars and vectors are present
      status = status && (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      status = status && this->mspRes;

      if(this->mspRes && (this->mSpaceId == Dimensions::Space::PHYSICAL))
      {
         status = status && (this->mMesh.size() == static_cast<size_t>(this->res().cpu()->nDim()));
      }

      return status;
   }

   void IVariableHdf5NWriter::addScalar(const std::pair<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalar)
   {
      // Only add the variable if it is listed as expected
      if(this->mExpected.count(scalar.first))
      {
         this->mScalars.insert(scalar);

         // Resolution is not yet set extract from scalar
         if(!this->mspRes)
         {
            std::visit([&](auto&& p){this->setResolution(p->dom(0).spRes());}, scalar.second);
         }
      }
   }

   void IVariableHdf5NWriter::addVector(const std::pair<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vector)
   {
      // Only add the variable if it is listed as expected
      if(this->mExpected.count(vector.first))
      {
         this->mVectors.insert(vector);

         // Resolution is not yet set extract from vector
         if(!this->mspRes)
         {
            std::visit([&](auto&& p){this->setResolution(p->dom(0).spRes());}, vector.second);
         }
      }
   }

   IVariableHdf5NWriter::scalar_iterator_range  IVariableHdf5NWriter::scalarRange()
   {
      return std::make_pair(this->mScalars.begin(), this->mScalars.end());
   }

   IVariableHdf5NWriter::vector_iterator_range  IVariableHdf5NWriter::vectorRange()
   {
      return std::make_pair(this->mVectors.begin(), this->mVectors.end());
   }

   void IVariableHdf5NWriter::setDatasetSize()
   {
      // Get dimensions ordered by index access speed (fast -> slow)
      ArrayI oDims;
      if(this->mSpaceId == Dimensions::Space::SPECTRAL)
      {
         oDims = this->res().counter().orderedDimensions(this->mSpaceId);
      }
      else
      {
         oDims = this->res().counter().orderedDimensions(this->mSpaceId);
      }

      int nDims = oDims.size();
      for(int i = 0; i < nDims; ++i)
      {
         // Set the dimension size in reverse order for HDF5
         this->mRegularDims.push_back(oDims(nDims-1-i));
      }
   }

   void IVariableHdf5NWriter::setDatasetOffsets()
   {
      // Compute the offsets
      if(this->mSpaceId == Dimensions::Space::SPECTRAL)
      {
         this->res().counter().computeOffsets(this->mBlock, this->mFileOffsets, this->mSpaceId);
      }
      else
      {
         this->res().counter().computeOffsets(this->mBlock, this->mFileOffsets, this->mSpaceId);
      }
   }

   void IVariableHdf5NWriter::setCollIo()
   {
      this->mCollIoWrite = 0;

      // Select transform dimension depending on dimension space
      Dimensions::Transform::Id  transId;
      if(this->mSpaceId == Dimensions::Space::SPECTRAL || this->mSpaceId == Dimensions::Space::TRANSFORM)
      {
         transId = Dimensions::Transform::TRA1D;
      } else
      {
         transId = static_cast<Dimensions::Transform::Id>(this->res().sim().ss().dimension()-1);
      }

      // Get the maximum number of slowest directions over all CPUs
      for(int i = 0; i < this->res().nCpu(); ++i)
      {
         if(this->mCollIoWrite < this->res().cpu(i)->dim(transId)->dim<Dimensions::Data::DAT3D>())
         {
            this->mCollIoWrite = this->res().cpu(i)->dim(transId)->dim<Dimensions::Data::DAT3D>();
         }
      }
   }

   void IVariableHdf5NWriter::writeTruncation()
   {
      std::ostringstream   oss;

      // Create the truncation parameters group
      hid_t base = H5Gcreate(this->file(), Tags::VariableHdf5::TRUNCATION.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Create the spectral truncation subgroup
      hid_t subGroup = H5Gcreate(base, Tags::VariableHdf5::TRUNCSPECTRAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write spectral resolution information
      for(int i = 0; i < this->res().cpu()->nDim(); i++)
      {
         oss << Tags::VariableHdf5::TRUNCDIM << i+1 << "D";

         // Write truncation value to file
         this->writeScalar(subGroup, oss.str(), this->res().sim().dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::SPECTRAL)-1);

         oss.str("");
      }

      // close group
      H5Gclose(subGroup);

      // Create the transform truncation subgroup
      subGroup = H5Gcreate(base, Tags::VariableHdf5::TRUNCTRANSFORM.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write transform resolution information
      for(int i = 0; i < this->res().cpu()->nDim(); i++)
      {
         oss << Tags::VariableHdf5::TRUNCDIM << i+1 << "D";

         // Write truncation value to file
         this->writeScalar(subGroup, oss.str(), this->res().sim().dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::TRANSFORM)-1);

         oss.str("");
      }

      // close group
      H5Gclose(subGroup);

      // Create the physical truncation subgroup
      subGroup = H5Gcreate(base, Tags::VariableHdf5::TRUNCPHYSICAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write physical resolution information
      for(int i = 0; i < this->res().cpu()->nDim(); i++)
      {
         oss << Tags::VariableHdf5::TRUNCDIM << i+1 << "D";

         // Write truncation value to file
         this->writeScalar(subGroup, oss.str(), this->res().sim().dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::PHYSICAL)-1);

         oss.str("");
      }

      // close group
      H5Gclose(subGroup);

      // close group
      H5Gclose(base);
   }

   void IVariableHdf5NWriter::writePhysical()
   {
      // Create the Physical parameters group
      hid_t group = H5Gcreate(this->file(), Tags::VariableHdf5::PHYSICAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write NonDimensional numbers
      for(auto it = this->mPhysical.cbegin(); it != this->mPhysical.cend(); ++it)
      {
         this->writeScalar(group, it->second->tag(), it->second->value());
      }

      // Write boundary flags
      for(auto it = this->mBoundary.cbegin(); it != this->mBoundary.cend(); ++it)
      {
         // Write reached simulation time to file
         this->writeScalar(group, it->first, it->second);
      }

      // close group
      H5Gclose(group);
   }

   void IVariableHdf5NWriter::writeRun()
   {
      // Create the Run parameters group
      hid_t group = H5Gcreate(this->file(), Tags::VariableHdf5::RUN.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write reached simulation time to file
      this->writeScalar(group, Tags::VariableHdf5::RUNTIME, this->mTime);

      // Write last timestep to file
      this->writeScalar(group, Tags::VariableHdf5::RUNSTEP, this->mTimestep);

      // close group
      H5Gclose(group);
   }

}
}
}
