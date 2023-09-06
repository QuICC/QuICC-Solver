/** 
 * @file FieldProbeWriter.hpp
 * @brief Implementation of the physical space field probe
 */

#ifndef QUICC_IO_VARIABLE_FIELDPROBEWRITER_HPP
#define QUICC_IO_VARIABLE_FIELDPROBEWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the physical space field proble
    */
   class FieldProbeWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         FieldProbeWriter(const std::string& prefix, const std::string& type, const std::vector<MHDFloat>& position);

         /**
          * @brief Destructor
          */
         virtual ~FieldProbeWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const; 
         
      protected:
         /**
          * @brief Write State to file
          */
         virtual void writeContent();

      private:
         /**
          * @brief Field value
          */
         std::vector<MHDFloat> mValue;

         /*
          * @brief Physical grid position
          */
         std::vector<MHDFloat> mPosition;

         /**
          * @brief Grid indexes corresponding to grid position
          */
         std::vector<int> mIndexes;
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<FieldProbeWriter> SharedFieldProbeWriter;

   inline bool FieldProbeWriter::isHeavy() const
   {
      return false;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_FIELDPROBEWRITER_HPP
