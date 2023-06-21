/**
 * @file Hdf5Typedefs.hpp
 * @brief Definition of simple HDF5 types to avoid including hdf5.h
 */

#ifndef QUICC_IO_HDF5TYPEDEFS_HPP
#define QUICC_IO_HDF5TYPEDEFS_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Io {

#ifdef QUICC_HDF5_V14_HSIZE_T
   /// Typedef for HDF5 hsize_t type. HDF5 version >= 1.14
   typedef long unsigned int QuICC_hsize_t;
#else
   /// Typedef for HDF5 hsize_t type. HDF5 version < 1.14
   typedef long long unsigned int QuICC_hsize_t;
#endif //QUICC_HDF5_V14_HSIZE_T


} // Io
} // QuICC

#endif // QUICC_IO_HDF5TYPEDEFS_HPP
