/**
 * @file MpiTypes.hpp
 * @brief Definition of a simple MPI datatype converters
 */

// Only define in MPI case
#ifdef QUICC_MPI

#ifndef QUICC_PARALLEL_MPITYPES_HPP
#define QUICC_PARALLEL_MPITYPES_HPP

// System includes
//
#include <complex>
#include <mpi.h>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Parallel {

   /**
    * @brief This class implements a few type conversion routines for the MPI library
    */
   struct MpiTypes
   {
      /**
       * @brief Get MPI predefined type
       */
      template <typename T> static MPI_Datatype type();
   };

   /**
    * @brief Specialised method for integer type
    */
   template <> inline MPI_Datatype MpiTypes::type<int>()
   {
      return MPI_INT;
   }

   /**
    * @brief Specialised method for long type
    */
   template <> inline MPI_Datatype MpiTypes::type<long>()
   {
      return MPI_LONG;
   }

   /**
    * @brief Specialised method for unsigned integer type
    */
   template <> inline MPI_Datatype MpiTypes::type<unsigned int>()
   {
      return MPI_UNSIGNED;
   }

   /**
    * @brief Specialised method for unsigned long type
    */
   template <> inline MPI_Datatype MpiTypes::type<unsigned long>()
   {
      return MPI_UNSIGNED_LONG;
   }

   /**
    * @brief Specialised method for float type
    */
   template <> inline MPI_Datatype MpiTypes::type<float>()
   {
      return MPI_FLOAT;
   }

   /**
    * @brief Specialised method for double type
    */
   template <> inline MPI_Datatype MpiTypes::type<double>()
   {
      return MPI_DOUBLE;
   }

   /**
    * @brief Specialised method for a complex<double>
    */
   template <> inline MPI_Datatype MpiTypes::type<std::complex<double> >()
   {
      return MPI_DOUBLE_COMPLEX;
   }

}
}

#endif // QUICC_PARALLEL_MPITYPES_HPP

// Only define in MPI case
#endif //QUICC_MPI
