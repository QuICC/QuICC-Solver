/** 
 * @file VisualizationFileWriter.hpp
 * @brief Implementation of the HDF5 visualisation file writer
 */

#ifndef QUICC_IO_VARIABLE_VISUALIZATIONFILEWRITER_HPP
#define QUICC_IO_VARIABLE_VISUALIZATIONFILEWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Variable/IVariableHdf5NWriter.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the HDF5 visualisation file writer
    */
   class VisualizationFileWriter: public IVariableHdf5NWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type    Type of the file (typically scheme name)
          */
         VisualizationFileWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~VisualizationFileWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:
         /**
          * @brief Write the mesh to file
          */
         void writeMesh();

         /**
          * @brief Create group for scalar field
          *
          * @param name    Name of the field
          * @param scalar  Scalar field values
          */
         void writePhysicalScalar(const std::string& name, const Framework::Selector::PhysicalScalarField& scalar);

         /**
          * @brief Create group for vector field
          *
          * @param name    Name of the field
          * @param vector  Vector of components
          */
         void writePhysicalVector(const std::string& name, const std::map<FieldComponents::Physical::Id,Framework::Selector::PhysicalScalarField>& vector, const std::string& joint = "_");

         /**
          * @brief Create group for tensor field
          *
          * @param name    Name of the field
          * @param tensor  Vector of components
          */
         void writePhysicalTensor(const std::string& name, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Framework::Selector::PhysicalScalarField>& tensor, const std::string& joint = "_");

      private:

   };

   /// Typedef for a shared pointer of a Hdf5Writer
   typedef std::shared_ptr<VisualizationFileWriter> SharedVisualizationFileWriter;

}
}
}

#endif // QUICC_IO_VARIABLE_VISUALIZATIONFILEWRITER_HPP
