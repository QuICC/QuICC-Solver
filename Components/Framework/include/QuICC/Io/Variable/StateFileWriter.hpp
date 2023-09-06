/** 
 * @file StateFileWriter.hpp
 * @brief Implementation of the HDF5 state file writer
 */

#ifndef QUICC_IO_VARIABLE_STATEFILEWRITER_HPP
#define QUICC_IO_VARIABLE_STATEFILEWRITER_HPP

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
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the HDF5 state file writer
    */
   class StateFileWriter: public IVariableHdf5NWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type       Type of the file (typically scheme name)
          * @param isRegular  Is data regular?
          */
         StateFileWriter(std::string type, bool isRegular);

         /**
          * @brief Destructor
          */
         virtual ~StateFileWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:
         /**
          * @brief Create group for scalar field
          *
          * @param name    Name of the field
          * @param scalar  Scalar field values
          */
         template <typename T> void writeSpectralScalar(const std::string& name, const typename Framework::Selector::ScalarField<T>& scalar);

         /**
          * @brief Create group for vector field
          *
          * @param name    Name of the field
          * @param vector  Vector of components
          */
         template <typename T> void writeSpectralVector(const std::string& name, const std::map<FieldComponents::Spectral::Id,typename Framework::Selector::ScalarField<T> >& vector);

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<StateFileWriter> SharedStateFileWriter;

   template <typename T> void StateFileWriter::writeSpectralScalar(const std::string& name, const typename Framework::Selector::ScalarField<T>& scalar)
   {
      // Create the scalar group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tuple<int,int,const T *> > fieldInfo = Datatypes::FieldTools::createInfo(scalar);

      // Check for data regularity
      if(this->mIsRegular)
      {
         // Write the scalar values
         this->writeRegularField(group, name, fieldInfo);
      } else
      {
         // Write the scalar values
         this->writeIrregularField(group, name, fieldInfo);
      }
      
      // close group
      H5Gclose(group);
   }

   template <typename T >void StateFileWriter::writeSpectralVector(const std::string& name, const std::map<FieldComponents::Spectral::Id,typename Framework::Selector::ScalarField<T> >& vector)
   {
      // Create the vector field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tuple<int,int, const T *> > fieldInfo;

      // Check for data regularity
      if(this->mIsRegular)
      {
         for(auto it = vector.cbegin(); it != vector.cend(); ++it)
         {
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(it->second);

            // Write vector component
            this->writeRegularField(group, name+"_"+Tools::IdToHuman::toTag(it->first), fieldInfo);
         }
      } else
      {
         for(auto it = vector.cbegin(); it != vector.cend(); ++it)
         {
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(it->second);

            // Write vector component
            this->writeIrregularField(group, name+"_"+Tools::IdToHuman::toTag(it->first), fieldInfo);
         }
      }
      
      // close group
      H5Gclose(group);
   }

}
}
}

#endif // QUICC_IO_VARIABLE_STATEFILEWRITER_HPP
