/** 
 * @file IHdf5Writer.cpp
 * @brief Source of the general HDF5 writer
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
#include "QuICC/Io/Hdf5/IHdf5Writer.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Hdf5 {

   IHdf5Writer::IHdf5Writer(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : Hdf5File(name, ext, header, type, version), mCollIoWrite(10000)
   {
   }

   IHdf5Writer::~IHdf5Writer()
   {
   }

   void IHdf5Writer::open()
   {
      // Create file
      hid_t fId = H5Fcreate(this->filename().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, this->filePList());

      // Check for successfully opened file
      if(fId < 0)
      {
         throw std::logic_error("Failed to open HDF5 file " + this->filename() + " in write mode!");
      }

      // Store the file handle
      this->setFile(fId);
   }

   void IHdf5Writer::close()
   {
      hid_t fPList = H5Fget_access_plist(this->file());

      // Close file
      H5Fclose(this->file());

      // Free the file property list
      this->freePList(fPList);
   }

   void IHdf5Writer::createFileInfo()
   {
      // Select root of the file
      hid_t loc = H5Gopen(this->file(), "/", H5P_DEFAULT);

      // Define some used variables
      hid_t attr;
      hid_t dspace;

      // Get the string data type
      hid_t type = H5Tcopy(H5T_C_S1);

      // Create header attribute
      dspace = H5Screate(H5S_SCALAR);
      H5Tset_size(type, this->header().size());
      attr = H5Acreate(loc, this->headerTag().c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr, type, this->header().c_str());
      H5Aclose(attr);
      H5Sclose(dspace);

      // Create type attribute
      dspace = H5Screate(H5S_SCALAR);
      H5Tset_size (type, this->type().size());
      attr = H5Acreate(loc, this->typeTag().c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr, type, this->type().c_str());
      H5Aclose(attr);
      H5Sclose(dspace);

      // Create version attribute
      dspace = H5Screate(H5S_SCALAR);
      H5Tset_size (type, this->version().size());
      attr = H5Acreate(loc, this->versionTag().c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr, type, this->version().c_str());
      H5Aclose(attr);
      H5Sclose(dspace);

      H5Gclose(loc);
   }

   void IHdf5Writer::writeString(hid_t loc, const std::string dsname, const std::string str)
   {
      // Get correct data type
      hid_t type = H5Tcopy(Hdf5Types::type<std::string>());
      H5Tset_size(type, str.size());

      // Compute size of the dataspace
      hsize_t dims = 1;

      // Create file dataspace
      hid_t filespace = H5Screate_simple(1, &dims, NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Set memory space
      hid_t memspace = H5Screate_simple(1, &dims, NULL);

      // Write memory
      if(QuICCEnv().allowsIO())
      {
         H5Dwrite(dataset, type, memspace, filespace, H5P_DEFAULT, str.c_str());
      }

      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }
}
}
}
