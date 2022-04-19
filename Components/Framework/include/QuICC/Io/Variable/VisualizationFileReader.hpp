/** 
 * @file VisualizationFileReader.hpp
 * @brief Implementation of the HDF5 visualisation file reader
 */

#ifndef QUICC_IO_VARIABLE_VISUALIZATIONFILEREADER_HPP
#define QUICC_IO_VARIABLE_VISUALIZATIONFILEREADER_HPP

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
#include "QuICC/Io/Variable/IVariableHdf5Reader.hpp"
#include "QuICC/Framework/Selector/ScalarField.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the HDF5 visualisation file reader
    */
   class VisualizationFileReader: public IVariableHdf5Reader
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name    Name of the file
          * @param type    Type of the file (typically scheme name)
          */
         VisualizationFileReader(std::string name, std::string type);

         /**
          * @brief Destructor
          */
         virtual ~VisualizationFileReader();

         /**
          * @brief Reader data to file
          */
         virtual void read();

         /**
          * @brief Read resolution and data information from file
          */
         void readSetup();
         
      protected:
         /**
          * @brief Write the mesh to file
          */
         void readMesh();

         /**
          * @brief Read scalar field from file
          *
          * @param name       Name of the field
          * @param rScalar    Storage for the scalar field
          * @param isRequired Field is required
          */
         void readPhysicalScalar(const std::string& name, Framework::Selector::PhysicalScalarField& rScalar, const bool isRequired);

         /**
          * @brief Read vector field from file
          *
          * @param name       Name of the field
          * @param rVector    Storage for the vector field
          * @param isRequired Field is required
          */
         void readPhysicalVector(const std::string& name, std::map<FieldComponents::Physical::Id,Framework::Selector::PhysicalScalarField>& rVector, const bool isRequired);

         /**
          * @brief Read vector field component from file
          *
          * @param name    Name of the vector field
          * @param id      ID of the component
          * @param rComp   Storage for field component
          * @param isRequired field is required
          */
         void readPhysicalComponent(const std::string& name, FieldComponents::Physical::Id id, Framework::Selector::PhysicalScalarField& rComp, const bool isRequired);

         /**
          * @brief Read tensor field from file
          *
          * @param name       Name of the field
          * @param rTensor    Storage for the tensor field
          * @param isRequired Field is required
          */
         void readPhysicalTensor(const std::string& name, std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Framework::Selector::PhysicalScalarField>& rTensor, const bool isRequired);

         /**
          * @brief Adapt data after reading it in if necessary
          *
          * @param rField   Storage for field
          */
         void adaptData(Framework::Selector::PhysicalScalarField& rComp);

      private:

   };

   /// Typedef for a shared pointer of a Hdf5Writer
   typedef std::shared_ptr<VisualizationFileReader> SharedVisualizationFileReader;

}
}
}

#endif // QUICC_IO_VARIABLE_VISUALIZATIONFILEREADER_HPP
