/** 
 * @file Cartesian1DTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII Chebyshev calculation for a Toroidal/Poloidal field in a plane layer
 */

#ifndef QUICC_IO_VARIABLE_CARTESIAN1DTORPOLENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_CARTESIAN1DTORPOLENERGYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/ICartesian1DTorPolEnergyWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the ASCII Chebyshev energy calculation for a Toroidal/Poloidal field in a plane layer
    */
   class Cartesian1DTorPolEnergyWriter: public ICartesian1DTorPolEnergyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DTorPolEnergyWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DTorPolEnergyWriter> SharedCartesian1DTorPolEnergyWriter;

} // Variable
} // Io
} // QuICC

#endif // QUICC_IO_VARIABLE_CARTESIAN1DTORPOLENERGYWRITER_HPP
