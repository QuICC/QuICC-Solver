/**
 * @file IVariableHdf5Reader.hpp
 * @brief Implementation of a generic variable data file reader
 */

#ifndef QUICC_IO_VARIABLE_IVARIABLEHDF5READER_HPP
#define QUICC_IO_VARIABLE_IVARIABLEHDF5READER_HPP

// System includes
//
#include <set>
#include <memory>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Hdf5/IHdf5Reader.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of a generic variable data file reader
    */
   class IVariableHdf5Reader: public Io::Hdf5::IHdf5Reader
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       File name
          * @param ext        File extension
          * @param header     File header
          * @param type       Type string of file
          * @param version    Version string of file
          * @param id         ID of the dimension space
          * @param isRegular  Is data regular?
          */
         IVariableHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const bool isRegular);

         /**
          * @brief Destructor
          */
         virtual ~IVariableHdf5Reader();

         /**
          * @brief Add name of expected variable to be added
          */
         void expect(const std::size_t id, const bool isRequired = true);

         /**
          * @brief Make sure all the expected variables have been added
          */
         bool isFull() const;

         /**
          * @brief Add scalar variable to file
          *
          * @param scalar Scalar variable to add
          */
         void addScalar(const std::pair<std::size_t,Framework::Selector::VariantSharedScalarVariable>& scalar);

         /**
          * @brief Add vector variable to file
          *
          * @param vector Vector variable to add
          */
         void addVector(const std::pair<std::size_t,Framework::Selector::VariantSharedVectorVariable>& vector);

         /**
          * @brief Write State to file
          */
         virtual void read() = 0;

         /**
          * @brief Get file truncation information
          */
         SharedSimulationResolution getFileTruncation() const;

         /**
          * @brief Read truncation information
          */
         void readTruncation();

         /**
          * @brief Read physical parameters from file
          */
         void getParameters(std::map<std::string, MHDFloat>& rParams) const;

      protected:
         /// Typedef for the scalar const iterator
         typedef std::map<std::size_t,Framework::Selector::VariantSharedScalarVariable>::iterator  scalar_iterator;

         /// Typedef for the vector const iterator
         typedef std::map<std::size_t,Framework::Selector::VariantSharedVectorVariable>::iterator  vector_iterator;

         /// Typedef for the scalar iterator range
         typedef std::pair<scalar_iterator,scalar_iterator>  scalar_iterator_range;

         /// Typedef for the scalar iterator range
         typedef std::pair<vector_iterator, vector_iterator>  vector_iterator_range;

         /**
          * @brief Get resolution
          */
         const Resolution& res() const;

         /**
          * @brief Set the read arguments
          */
         void setReadArguments();

         /**
          * @brief Check truncation compatibily between data and file
          */
         void checkTruncation();

         /**
          * @brief Get iterator range to scalars
          */
         scalar_iterator_range   scalarRange();

         /**
          * @brief Get iterator range to vectors
          */
         vector_iterator_range   vectorRange();

         /**
          * @brief Is file working on regular data?
          */
         bool mIsRegular;

         /**
          * @brief Is optional variable?
          */
         bool isRequired(const std::size_t id);

      private:
         /**
          * @brief Set the resolution and use it for preliminary initialisation
          *
          * @param spRes      Resolution information
          */
         void setResolution(SharedResolution spRes);

         /**
          * @brief Set the maximum number of IO operations
          */
         void setCollIo();

         /**
          * @brief The dimension space the file is working on
          */
         Dimensions::Space::Id mSpaceId;

         /**
          * @brief Simulation resolution based on file parameters
          */
         SharedSimulationResolution mspFileRes;

         /**
          * @brief Storage for the names of the expected variables
          */
         std::set<std::size_t>  mExpected;

         /**
          * @brief Storage for the names of the required variables
          */
         std::set<std::size_t>  mRequired;

         /**
          * @brief Resolution information
          */
         SharedResolution mspRes;

         /**
          * @brief Storage for the scalars
          */
         std::map<std::size_t,Framework::Selector::VariantSharedScalarVariable>   mScalars;

         /**
          * @brief Storage for the vectors
          */
         std::map<std::size_t,Framework::Selector::VariantSharedVectorVariable>   mVectors;
   };

   /// Typedef for a smart reference counting pointer of a Variable HDF5 numbering reader
   typedef std::shared_ptr<IVariableHdf5Reader>   SharedIVariableHdf5Reader;

}
}
}

#endif // QUICC_IO_VARIABLE_IVARIABLEHDF5READER_HPP
