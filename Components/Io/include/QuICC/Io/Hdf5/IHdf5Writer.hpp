/**
 * @file IHdf5Writer.hpp
 * @brief Interface to a general HDF5 writer
 */

#ifndef QUICC_IO_HDF5_IHDF5WRITER_HPP
#define QUICC_IO_HDF5_IHDF5WRITER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <tuple>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/Io/Hdf5/Hdf5Types.hpp"
#include "QuICC/Io/Hdf5/Hdf5File.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Io {

namespace Hdf5 {

   /**
    * @brief Interface to a general HDF5 file writer
    */
   class IHdf5Writer: public Hdf5File
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name     Filename
          * @param ext      File extension
          * @param header   Header string of file
          * @param type     Type string of file
          * @param version  Version string of file
          */
         IHdf5Writer(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
          * @brief Destructor
          */
         virtual ~IHdf5Writer();

         /**
          * @brief Initialise the file
          */
         virtual void init() = 0;

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize() = 0;

      protected:
         /**
          * @brief Number of blocks to read
          */
         std::vector<hsize_t>  mBlock;

         /**
          * @brief Max collective IO Write operations over all CPUs
          */
         int mCollIoWrite;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Create the file info (header, type and version)
          */
         void createFileInfo();

         /**
          * @brief Write std::string dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param str     String to write
          */
         void writeString(hid_t loc, const std::string dsname, const std::string str);

         /**
          * @brief Write scalar dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param val     Scalar to write
          *
          * \tparam T Type of the scalar
          */
         template <typename T> void writeScalar(hid_t loc, const std::string dsname, const T val);

         /**
          * @brief Write Array dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param arr     Array of values to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeArray(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, 1>& arr);

         /**
          * @brief Write Matrix dataset
          *
          * \warning In the MPI case only the IO node is going to write data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param mat     Matrix of values to write
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeMatrix(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat);

         /**
          * @brief Write an irregular field vector dataset
          *
          * For irregular field data the values gets stored in a 2D array
          *
          * \warning This method does not work for 1D data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Storage for the data to read (0 -> rows, 1 -> cols, 2 -> const data pointer)
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeIrregularField(hid_t loc, const std::string dsname, const std::vector<std::tuple<int, int, const T *> >& storage);

         /**
          * @brief Write a regular field vector dataset (on sliced map)
          *
          * For regular field data the values gets stored in a array of the same
          * dimensionality as the data i.e. 3D in 3D array, 2D data in 2D array
          *
          * \warning This method does not work for 1D data
          *
          * @param loc     Base location in HDF5 file
          * @param dsname  HDF5 dataset name
          * @param storage Storage for the data to read (0 -> rows, 1 -> cols, 2 -> const data pointer)
          *
          * \tparam T Type of the scalars
          */
         template <typename T> void writeRegularField(hid_t loc, const std::string dsname, const std::vector<std::tuple<int, int, const T *> >& storage);

      private:
   };

   template <typename T> void IHdf5Writer::writeScalar(hid_t loc, const std::string dsname, const T val)
   {
      // Get correct data type
      hid_t type = Hdf5Types::type<T>();

      // Create scalar dataspace
      hid_t dspace = H5Screate(H5S_SCALAR);

      // Create data set
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write scalar value
      if(QuICCEnv().allowsIO())
      {
         H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
      }

      // close dataspace
      H5Sclose(dspace);

      // close dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Writer::writeArray(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, 1>& arr)
   {
      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Compute size of the dataspace
      hsize_t dims = arr.size();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(1, &dims, NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Set memory space
      hid_t memspace = H5Screate_simple(1, &dims, NULL);

      // Write memory
      if(QuICCEnv().allowsIO())
      {
         H5Dwrite(dataset, type, memspace, filespace, H5P_DEFAULT, arr.data());
      }

      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Writer::writeMatrix(hid_t loc, const std::string dsname, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat)
   {
      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Compute size of the dataspace
      hsize_t dims[2];
      dims[0] = mat.cols();
      dims[1] = mat.rows();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(2, dims, NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Set memory space
      hid_t memspace = H5Screate_simple(2, dims, NULL);

      // Write memory
      if(QuICCEnv().allowsIO())
      {
         H5Dwrite(dataset, type, memspace, filespace, H5P_DEFAULT, mat.data());
      }

      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(dataset);
   }

   template <typename T> void IHdf5Writer::writeRegularField(hid_t loc, const std::string dsname, const std::vector<std::tuple<int, int, const T *> >& storage)
   {
      Profiler::RegionFixture<4> fix("IHdf5Writer::writeRegularField");

      // Get dimensionality of file data
      int nDims = this->mRegularDims.size();

      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(nDims, &this->mRegularDims.front(), NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffsets[nDims];
      // Fix offset for the fastest dimension (always full on each CPU)
      pOffsets[nDims-1] = 0;

      // Compute size of the dataspace
      hsize_t dims[nDims];
      // Always write one block of the slowest dimension
      dims[0] = 1;
      // Fastest dimension is always simulation dimension (even in MPI mode)
      dims[nDims-1] = this->mRegularDims.back();

      // Compute size of the memory dataspace (2D blocks)
      hsize_t iDims[2];
      // Fast dimension is always simulation dimension (even in MPI mode)
      iDims[1] = this->mRegularDims.back();

      // Set PList to parallel access
      hid_t dsPList = this->datasetPList(true);

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space (number of columns of block)
         iDims[0] = std::get<1>(storage.at(i));
         memspace = H5Screate_simple(2, iDims, NULL);

         // Set offsets
         for(int j = 0; j < nDims -1; j++)
         {
            pOffsets[j] = this->mFileOffsets.at(i).at(j);
         }

         // Select corresponding hyperslab (medium dimension is number of block columns)
         dims[nDims-2] = std::get<1>(storage.at(i));
         H5Sselect_hyperslab(filespace, H5S_SELECT_SET, pOffsets, NULL, dims, NULL);

         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, std::get<2>(storage.at(i)));

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Add some ZERO IO calls to allow for collective writes
      //

      // Set zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Select zero file space
      H5Sselect_none(filespace);

      // Zero write IO
      for(unsigned int i = 0; i < (this->mCollIoWrite - storage.size()) ; ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, std::get<2>(storage.at(0)));
      }

      // Close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(dataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freePList(dsPList);
   }

   template <typename T> void IHdf5Writer::writeIrregularField(hid_t loc, const std::string dsname, const std::vector<std::tuple<int, int, const T *> >& storage)
   {
      Profiler::RegionFixture<4> fix("IHdf5Writer::writeIrregularField");

      // Get dimensionality of file data
      int nDims = this->mRegularDims.size();

      // Set data type correctly
      hid_t type = Hdf5Types::type<T>();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(nDims, &(this->mRegularDims.front()), NULL);

      // Create dataset in file
      hid_t dataset = H5Dcreate(loc, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace
      hid_t  memspace;

      // Create offset storage
      hsize_t pOffset[2];
      pOffset[1] = 0;

      // Compute size of the memory dataspace
      hsize_t iDims[2];

      // Set PList to parallel access
      hid_t dsPList = this->datasetPList(false);

      // Storage for the memory offsets
      hsize_t memOffset[2];
      memOffset[0] = 0;
      memOffset[1] = 0;

      // Get data dimensions
      int memRows = 0;
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         memRows = std::max(memRows, std::get<0>(storage.at(i)));
      }

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Create rectangular memory space
         iDims[0] = std::get<1>(storage.at(i));
         iDims[1] = memRows;
         memspace = H5Screate_simple(2, iDims, NULL);

         // Select columns
         H5Sselect_none(memspace);
         iDims[1] = this->mBlock.at(i);
         iDims[0] = 1;
         for(int j = 0; j < std::get<1>(storage.at(i)); ++j)
         {
            memOffset[0] = j;
            H5Sselect_hyperslab(memspace, H5S_SELECT_OR, memOffset, NULL, iDims, NULL);
         }

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Write memory into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, std::get<2>(storage.at(i)));

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      //
      // Add some ZERO IO calls to allow for collective writes
      //

      // Create zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Create zero file selection space
      H5Sselect_none(filespace);

      // Zero write IO
      for(unsigned int i = 0; i < (this->mCollIoWrite - storage.size()); ++i)
      {
         // Write empty data into hyperslab
         H5Dwrite(dataset, type, memspace, filespace, dsPList, NULL);
      }

      // close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(dataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freePList(dsPList);
   }

   /// Typedef for a shared pointer of a IHdf5Writer
   typedef std::shared_ptr<IHdf5Writer> SharedIHdf5Writer;

}
}
}

#endif // QUICC_IO_HDF5_IHDF5WRITER_HPP
