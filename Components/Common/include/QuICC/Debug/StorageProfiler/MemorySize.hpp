/** 
 * @file MemorySize.hpp
 * @brief Simple template for the memory usage of the different types
 */

#ifndef QUICC_DEBUG_MEMORYSIZE_HPP
#define QUICC_DEBUG_MEMORYSIZE_HPP

// Configuration includes
//

// System includes
//
#include <complex>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Debug {

   /**
    * @brief Simple template for the memory usage of the different types
    */
   template <typename T> struct MemorySize
   {
   };

   /**
    * @brief Template specialization for char
    */
   template <> struct MemorySize<char>
   {
      /// Size of datatype in bytes
      static const int BYTES = 1;
   };


   /**
    * @brief Template specialization for int
    */
   template <> struct MemorySize<int>
   {
      /// Size of datatype in bytes
      static const int BYTES = 4;
   };


   /**
    * @brief Template specialization for int
    */
   template <> struct MemorySize<long>
   {
      /// Size of datatype in bytes
      static const int BYTES = 8;
   };


   /**
    * @brief Template specialization for float
    */
   template <> struct MemorySize<float>
   {
      /// Size of datatype in bytes
      static const int BYTES = 4;
   };

   /**
    * @brief Template specialization for double
    */

   template <> struct MemorySize<double>
   {
      /// Size of datatype in bytes
      static const int BYTES = 8;
   };

   /**
    * @brief Template specialization for complex<float>
    */

   template <> struct MemorySize<std::complex<float> >
   {
      /// Size of datatype in bytes
      static const int BYTES = 8;
   };

   /**
    * @brief Template specialization for complex<double>
    */
   template <> struct MemorySize<std::complex<double> >
   {
      /// Size of datatype in bytes
      static const int BYTES = 16;
   };

}
}

#endif // QUICC_DEBUG_MEMORYSIZE_HPP
