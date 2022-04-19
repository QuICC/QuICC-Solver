/**
 * @file IVariableHdf5Reader.cpp
 * @brief Source of the implementation of a generic variable data file reader
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
#include "QuICC/Io/Variable/IVariableHdf5Reader.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Io/Variable/Tags/VariableHdf5.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   IVariableHdf5Reader::IVariableHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const bool isRegular)
      : IHdf5Reader(name, ext, header, type, version), mIsRegular(isRegular), mSpaceId(id)
   {
   }

   IVariableHdf5Reader::~IVariableHdf5Reader()
   {
   }

   const Resolution& IVariableHdf5Reader::res() const
   {
      return *this->mspRes;
   }

   void IVariableHdf5Reader::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;

      // This set collective IO operations
      this->setCollIo();
   }

   void IVariableHdf5Reader::expect(const std::size_t id, const bool isRequired)
   {
      this->mExpected.insert(id);

      if(isRequired)
      {
         this->mRequired.insert(id);
      }
   }

   bool IVariableHdf5Reader::isFull() const
   {
      // Check that all expected scalars and vectors are present
      bool sizeStatus = (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      bool resStatus = static_cast<bool>(this->mspRes);

      return (sizeStatus && resStatus);
   }

   void IVariableHdf5Reader::addScalar(const std::pair<std::size_t, Framework::Selector::VariantSharedScalarVariable>& scalar)
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

   void IVariableHdf5Reader::addVector(const std::pair<std::size_t, Framework::Selector::VariantSharedVectorVariable>& vector)
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

   IVariableHdf5Reader::scalar_iterator_range  IVariableHdf5Reader::scalarRange()
   {
      return std::make_pair(this->mScalars.begin(), this->mScalars.end());
   }

   IVariableHdf5Reader::vector_iterator_range  IVariableHdf5Reader::vectorRange()
   {
      return std::make_pair(this->mVectors.begin(), this->mVectors.end());
   }

   void IVariableHdf5Reader::setCollIo()
   {
      this->mCollIoRead = 0;
   }

   bool IVariableHdf5Reader::isRequired(const std::size_t id)
   {
      return this->mRequired.count(id);
   }

   void IVariableHdf5Reader::setReadArguments()
   {
      // Get dimensions ordered by index access speed (fast -> slow)
      ArrayI oDims = this->res().counter().orderedDimensions(this->mSpaceId);

      int nDims = oDims.size();
      for(int i = 0; i < nDims; ++i)
      {
         // Set the dimension size in reverse order for HDF5
         this->mRegularDims.push_back(oDims(nDims-1-i));
      }

      // Compute offsets
      this->res().counter().computeOffsets(this->mBlock, this->mFileOffsets, this->mSpaceId, this->mspFileRes);

      // Get the "global" local minimum for MPI code
      this->mCollIoRead = this->mFileOffsets.size();
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mCollIoRead, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      #endif // QUICC_MPI
   }

   SharedSimulationResolution IVariableHdf5Reader::getFileTruncation() const
   {
      std::ostringstream   oss;

      // Open the truncation parameters group
      hid_t base = H5Gopen(this->file(), Tags::VariableHdf5::TRUNCATION.c_str(), H5P_DEFAULT);

      // Open the spectral truncation subgroup
      hid_t subGroup = H5Gopen(base, Tags::VariableHdf5::TRUNCSPECTRAL.c_str(), H5P_DEFAULT);

      // Get number of dimensions
      hsize_t nDim = 0;
      H5Gget_num_objs(subGroup, &nDim);

      ArrayI spec(nDim);

      // Read spectral resolution information
      for(hsize_t i = 0; i < nDim; i++)
      {
         oss << Tags::VariableHdf5::TRUNCDIM << i+1 << "D";

         // Read dimension from file
         this->readScalar(subGroup, oss.str(), spec(i));

         oss.str("");
      }
      // Increment by one to get size
      spec.array() += 1;

      // close group
      H5Gclose(subGroup);

      // Open the transform truncation parameters group
      subGroup = H5Gopen(base, Tags::VariableHdf5::TRUNCTRANSFORM.c_str(), H5P_DEFAULT);

      ArrayI trans(nDim);

      // Read transform resolution information
      for(hsize_t i = 0; i < nDim; i++)
      {
         oss << Tags::VariableHdf5::TRUNCDIM << i+1 << "D";

         // Read dimension from file
         this->readScalar(subGroup, oss.str(), trans(i));

         oss.str("");
      }
      // Increment by one to get size
      trans.array() += 1;

      // close group
      H5Gclose(subGroup);

      // Open the physical truncation parameters group
      subGroup = H5Gopen(base, Tags::VariableHdf5::TRUNCPHYSICAL.c_str(), H5P_DEFAULT);

      ArrayI phys(nDim);

      // Read physical resolution information
      for(hsize_t i = 0; i < nDim; i++)
      {
         oss << Tags::VariableHdf5::TRUNCDIM << i+1 << "D";

         // Read dimension from file
         this->readScalar(subGroup, oss.str(), phys(i));

         oss.str("");
      }
      // Increment by one to get size
      phys.array() += 1;

      // close group
      H5Gclose(subGroup);

      // close group
      H5Gclose(base);

      auto spSimRes = std::make_shared<SimulationResolution>(phys,spec,trans);
      return spSimRes;
   }

   void IVariableHdf5Reader::getParameters(std::map<std::string, MHDFloat>& rParams) const
   {
      std::ostringstream   oss;

      // Open the truncation parameters group
      hid_t base = H5Gopen(this->file(), Tags::VariableHdf5::PHYSICAL.c_str(), H5P_DEFAULT);

      // Get number of parameters
      hsize_t nParams = 0;
      H5Gget_num_objs(base, &nParams);

      // Iterate over all parameters
      ssize_t size;
      std::string dsname = "";
      MHDFloat value;
      for(hsize_t i = 0; i < nParams; i++)
      {
         // Extract name of dataset
         size = H5Lget_name_by_idx(base, ".", H5_INDEX_NAME, H5_ITER_INC, static_cast<hsize_t>(i), NULL, 0, H5P_DEFAULT);
         dsname.resize(size);
         size = H5Lget_name_by_idx(base, ".", H5_INDEX_NAME, H5_ITER_INC, static_cast<hsize_t>(i), &dsname[0], size+1, H5P_DEFAULT);

         // Read value from dataset
         hid_t dataset = H5Dopen(base, dsname.c_str(), H5P_DEFAULT);
         hid_t type = H5Dget_type(dataset);
         if(Io::Hdf5::Hdf5Types::isSupported<MHDFloat>(type))
         {
            H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
         } else
         {
            throw std::logic_error("Tried to read unimplemented parameter data type");
         }
         H5Tclose(type);
         H5Dclose(dataset);

         // Add parameter to map
         rParams.insert(std::make_pair(dsname, value));
      }

      // close group
      H5Gclose(base);
   }

   void IVariableHdf5Reader::readTruncation()
   {
      this->mspFileRes = this->getFileTruncation();
   }

   void IVariableHdf5Reader::checkTruncation()
   {
      std::ostringstream   oss;

      for(int i = 0; i < this->res().cpu()->nDim(); i++)
      {
         if(this->res().sim().dim(static_cast<Dimensions::Simulation::Id>(i),this->mSpaceId) != this->mspFileRes->dim(static_cast<Dimensions::Simulation::Id>(i),this->mSpaceId))
         {
            oss << "Dimension " << i+1 <<  " doesn't fit";
            Tools::Formatter::printLine(std::cout, '-');
            Tools::Formatter::printCentered(std::cout, oss.str(), '*');

            if(this->res().sim().dim(static_cast<Dimensions::Simulation::Id>(i),this->mSpaceId) > this->mspFileRes->dim(static_cast<Dimensions::Simulation::Id>(i),this->mSpaceId))
            {
               Tools::Formatter::printCentered(std::cout, " ---> Zeros have been added!", ' ');
            } else
            {
               Tools::Formatter::printCentered(std::cout, "---> File data has been truncated!", ' ');
            }
            Tools::Formatter::printLine(std::cout, '-');
            Tools::Formatter::printNewline(std::cout);
            oss.str("");
         }
      }
      Tools::Formatter::printNewline(std::cout);
   }
}
}
}
