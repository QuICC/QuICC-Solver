/**
 * @file IHdf5Reader.cpp
 * @brief Source of the HDF5 reader implementation
 */

// System includes
//
#include <sstream>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Hdf5/IHdf5Reader.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Hdf5 {

   IHdf5Reader::IHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : Hdf5File(name, ext, header, type, version), mCollIoRead(0)
   {
   }

   IHdf5Reader::~IHdf5Reader()
   {
   }

   void IHdf5Reader::init()
   {
      // Open the file
      this->open();

      // Check file compatibility
      this->checkCompatibility();
   }

   void IHdf5Reader::open()
   {
      // As this may fail, disable error handler for H5FOpen call
      /* Save old error handler */
      H5E_auto_t err_func;
      void *err_client_data;
      H5Eget_auto(H5E_DEFAULT, &err_func, &err_client_data);
      H5Eset_auto(H5E_DEFAULT, NULL, NULL);

      // Open file
      hid_t fId = H5Fopen(this->filename().c_str(), H5F_ACC_RDONLY, this->filePList());

      // Enable error handling
      H5Eset_auto(H5E_DEFAULT, err_func, err_client_data);

      // Check for successfully opened file
      if(fId < 0)
      {
         throw std::logic_error("Failed to open HDF5 file " + this->filename() + " in read mode!");
      }

      // Store the file handle
      this->setFile(fId);
   }

   void IHdf5Reader::finalize()
   {
      this->close();
   }

   void IHdf5Reader::close()
   {
      hid_t fPList = H5Fget_access_plist(this->file());

      // Close file
      H5Fclose(this->file());

      // Free the property list
      this->freePList(fPList);
   }

   void IHdf5Reader::checkCompatibility()
   {
      // Open the root of the file
      hid_t loc = H5Gopen(this->file(), "/", H5P_DEFAULT);

      // Some useful variables
      hid_t type, ftype;
      hid_t  attr;
      char *cStr;
      size_t size;

      //
      // Read header attribute
      //

      // Open the header attribute
      attr = H5Aopen(loc, this->headerTag().c_str(), H5P_DEFAULT);
      // Get type
      ftype = H5Aget_type(attr);
      // Get string size
      size = H5Tget_size(ftype);
      // Create string type
      type = H5Tcopy (H5T_C_S1);
      H5Tset_strpad(type, H5Tget_strpad(ftype));
      // Set string length
      H5Tset_size(type, size);
      // Alloate memory
      cStr = new char[size];
      // Read attribute
      H5Aread(attr, type, cStr);
      // Create string
      std::string fileHeader(cStr, size);
      // Check compatibility
      if(this->header() != fileHeader)
      {
         throw std::logic_error("Wrong HDF5 file header!");
      }
      delete[] cStr;
      H5Aclose(attr);

      //
      // Read type attribute
      //

      // Open the type attribute
      attr = H5Aopen(loc, this->typeTag().c_str(), H5P_DEFAULT);
      // Get type
      ftype = H5Aget_type(attr);
      // Get string size
      size = H5Tget_size(ftype);
      // Create string type
      type = H5Tcopy (H5T_C_S1);
      H5Tset_strpad(type, H5Tget_strpad(ftype));
      // Set string length
      H5Tset_size(type, size);
      // Alloate memory
      cStr = new char[size];
      // Read attribute
      H5Aread(attr, type, cStr);
      // Create string
      std::string fileType(cStr, size);
      // Check compatibility
      if(this->type() != fileType)
      {
         throw std::logic_error("Wrong HDF5 file type!");
      }
      delete[] cStr;
      H5Aclose(attr);

      //
      // Read version attribute
      //

      // Open the version attribute
      attr = H5Aopen(loc, this->versionTag().c_str(), H5P_DEFAULT);
      // Get type
      ftype = H5Aget_type(attr);
      // Get string size
      size = H5Tget_size(ftype);
      // Create string type
      type = H5Tcopy (H5T_C_S1);
      H5Tset_strpad(type, H5Tget_strpad(ftype));
      // Set string length
      H5Tset_size(type, size);
      // Alloate memory
      cStr = new char[size];
      // Read attribute
      H5Aread(attr, type, cStr);
      // Create string
      std::string fileVersion(cStr, size);
      // Check compatibility
      if(this->version() != fileVersion)
      {
         throw std::logic_error("Wrong HDF5 file version!");
      }
      delete[] cStr;
      H5Aclose(attr);

      H5Gclose(loc);
   }
}
}
}
