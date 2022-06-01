/** 
 * @file Hdf5File.cpp
 * @brief Source of the implementation of a general HDF5 file
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Io/Hdf5/Hdf5File.hpp"

// Project includes
//

namespace QuICC {

namespace Io {

namespace Hdf5 {

   Hdf5File::Hdf5File(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : mName(name), mExt(ext), mHeader(header), mType(type), mVersion(version), mFile(-1)
   {
   }

   Hdf5File::~Hdf5File()
   {
      // Cleanup released resources
      H5garbage_collect();
   }

   hid_t Hdf5File::filePList()
   {
      // 
      // Default options are fine for serial case
      //

      #ifdef QUICC_MPI
         // Create file access property list
         hid_t fPList = H5Pcreate(H5P_FILE_ACCESS);

         #ifdef QUICC_ON_JANUS
            // Deactivate incompatible ROMIO features
            MPI_Info info;
            MPI_Info_create(&info);
            MPI_Info_set(info, (char *)"romio_ds_write", (char *)"disable");
            MPI_Info_set(info, (char *)"romio_ds_read", (char *)"disable");

            // Create the MPI IO access property
            H5Pset_fapl_mpio(fPList, MPI_COMM_WORLD, info);

            // Free info
            MPI_Info_free(&info);
         #else
            // Create the MPI IO access property
            H5Pset_fapl_mpio(fPList, MPI_COMM_WORLD, MPI_INFO_NULL);
         #endif // QUICC_ON_JANUS
      #else
         hid_t fPList(H5P_DEFAULT);
      #endif // QUICC_MPI

      return fPList;
   }

   hid_t Hdf5File::datasetPList(const bool isCollective)
   {
      // 
      // Default options are fine for serial case
      //

      #ifdef QUICC_MPI
         // Create dataset transfer property list
         hid_t dsPList = H5Pcreate(H5P_DATASET_XFER);

         // Set the transfer property to collective IO
         if(isCollective)
         {
            H5Pset_dxpl_mpio(dsPList, H5FD_MPIO_COLLECTIVE);
         } else
         {
            // ... or set the transfer property to independent IO
            H5Pset_dxpl_mpio(dsPList, H5FD_MPIO_INDEPENDENT);
         }
      #else
         hid_t dsPList(H5P_DEFAULT);
      #endif // QUICC_MPI

      return dsPList;
   }

   void Hdf5File::freePList(hid_t dsPList)
   {
      // Let the HDF5 library close the dataset property list
      H5Pclose(dsPList);
   }

   std::string Hdf5File::filename() const
   {
      return this->mName + this->mExt;
   }

   void Hdf5File::resetName(std::string name)
   {
      this->mName = name;
   }

   std::string Hdf5File::name() const
   {
      return this->mName;
   }

   std::string Hdf5File::extension() const
   {
      return this->mExt;
   }

   std::string Hdf5File::headerTag() const
   {
      return Hdf5File::HEADER_TAG;
   }

   std::string Hdf5File::header() const
   {
      return this->mHeader;
   }

   std::string Hdf5File::typeTag() const
   {
      return Hdf5File::TYPE_TAG;
   }

   std::string Hdf5File::type() const
   {
      return this->mType;
   }

   std::string Hdf5File::versionTag() const
   {
      return Hdf5File::VERSION_TAG;
   }

   std::string Hdf5File::version() const 
   {
      return this->mVersion;
   }

   hid_t Hdf5File::file() const
   {
      return this->mFile;
   }

   void Hdf5File::setFile(hid_t fId)
   {
      this->mFile = fId;
   }

   const std::string Hdf5File::HEADER_TAG = "header";

   const std::string Hdf5File::TYPE_TAG = "type";

   const std::string Hdf5File::VERSION_TAG = "version";
}
}
}
