/** 
 * @file Hdf5Types.hpp
 * @brief Implementation of a simple HDF5 datatype converter
 */

#ifndef QUICC_IO_HDF5_HDF5TYPES_HPP
#define QUICC_IO_HDF5_HDF5TYPES_HPP

// System includes
//
#include <complex>
#include <hdf5.h>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Io {

namespace Hdf5 {

   /**
    * @brief This class implements a few type conversion routines for the HDF5 library
    */
   struct Hdf5Types
   {
      /**
       * @brief Get HDF5 predefined type for a scalar type
       *
       * \tparam T type to convert to HDF5 type
       */
      template <typename T> static hid_t type();

      /**
       * @brief Check if HDF5 datatype is supported 
       *
       * \tparam T type to convert to HDF5 type
       * \param type    HDF5 datatype
       */
      template <typename T> static bool isSupported(hid_t type);

      /**
       * @brief Store complex values as a struct in HDF5 file
       */
      template <typename T> static hid_t complex_asStruct();

      /**
       * @brief Store complex values as an array in HDF5 file
       */
      template <typename T> static hid_t complex_asArray();
   };

   template <> inline hid_t Hdf5Types::complex_asStruct<float>()
   {
      std::complex<float> tmp(0,0);
      float d = 0;
      hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof tmp);
      H5Tinsert (complex_id, "r", 0, H5T_NATIVE_FLOAT);
      H5Tinsert (complex_id, "i", sizeof d, H5T_NATIVE_FLOAT);

      return complex_id;
   }

   template <> inline hid_t Hdf5Types::complex_asStruct<double>()
   {
      std::complex<double> tmp(0,0);
      double d = 0;
      hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof tmp);
      H5Tinsert (complex_id, "r", 0, H5T_NATIVE_DOUBLE);
      H5Tinsert (complex_id, "i", sizeof d, H5T_NATIVE_DOUBLE);

      return complex_id;
   }

   template <> inline hid_t Hdf5Types::complex_asArray<float>()
   {
      hsize_t dims = 2;
      return H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims);
   }

   template <> inline hid_t Hdf5Types::complex_asArray<double>()
   {
      hsize_t dims = 2;
      return H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &dims);
   }

   /**
    * @brief Specialised method for integer type
    */
   template <> inline hid_t Hdf5Types::type<int>()
   {
      return H5T_NATIVE_INT;
   }

   /**
    * @brief Specialised method for std:string type
    */
   template <> inline hid_t Hdf5Types::type<std::string>()
   {
      return H5T_C_S1;
   }

   /**
    * @brief Specialised method for float type
    */
   template <> inline hid_t Hdf5Types::type<float>()
   {
      return H5T_NATIVE_FLOAT;
   }

   /**
    * @brief Specialised method for double type
    */
   template <> inline hid_t Hdf5Types::type<double>()
   {
      return H5T_NATIVE_DOUBLE;
   }

   /**
    * @brief Specialised method for a complex<float>
    */
   template <> inline hid_t Hdf5Types::type<std::complex<float> >()
   {
      // Use simple 2D array for complex data
      #if defined QUICC_HDF5_CMPLX_ARRAY
         return Hdf5Types::complex_asArray<float>();
      // Use struct for complex data
      #elif defined QUICC_HDF5_CMPLX_STRUCT
         return Hdf5Types::complex_asStruct<float>();
      #endif //defined QUICC_HDF5_CMPLX_ARRAY
   }

   /**
    * @brief Specialised method for a complex<double>
    */
   template <> inline hid_t Hdf5Types::type<std::complex<double> >()
   {
      // Use simple 2D array for complex data
      #if defined QUICC_HDF5_CMPLX_ARRAY
         return Hdf5Types::complex_asArray<double>();
      // Use struct for complex data
      #elif defined QUICC_HDF5_CMPLX_STRUCT
         return Hdf5Types::complex_asStruct<double>();
      #endif //defined QUICC_HDF5_CMPLX_ARRAY
   }

   template <typename T> inline bool Hdf5Types::isSupported(hid_t type)
   {
      return true;
   }

   template <> inline bool Hdf5Types::isSupported<std::complex<float> >(hid_t type)
   {
      if(H5Tequal(type, Hdf5Types::complex_asArray<float>()) || H5Tequal(type, Hdf5Types::complex_asStruct<float>()))
      {
         return true;
      } else
      {
         return false;
      }
   }

   template <> inline bool Hdf5Types::isSupported<std::complex<double> >(hid_t type)
   {
      if(H5Tequal(type, Hdf5Types::complex_asArray<double>()) || H5Tequal(type, Hdf5Types::complex_asStruct<double>()))
      {
         return true;
      } else
      {
         return false;
      }
   }

}
}
}

#endif // QUICC_IO_HDF5_HDF5TYPES_HPP
