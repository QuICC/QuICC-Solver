/** 
 * @file IStatisticsAsciiWriter.hpp
 * @brief Implementation of a generic variable to ASCII file writer
 */

#ifndef QUICC_IO_STATS_ISTATISTICSASCIIWRITER_HPP
#define QUICC_IO_STATS_ISTATISTICSASCIIWRITER_HPP

// System includes
//
#include <set>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Io/Ascii/IAsciiWriter.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"

namespace QuICC {

namespace Io {

namespace Stats {

   /**
    * @brief Implementation of a generic variable to ASCII file writer
    */
   class IStatisticsAsciiWriter: public Io::Ascii::IAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       Filename
          * @param ext        File extension
          * @param header     Header string of file
          * @param type       Type string of file
          * @param version    Version string of file
          * @param id         ID of the dimension space
          * @param mode       Write mode of file
          */
         IStatisticsAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode = IAsciiWriter::EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~IStatisticsAsciiWriter();

         /**
          * @brief Add name of expected variable to be added
          *
          * @param id ID of field
          */
         void expect(const std::size_t id);

         /**
          * @brief Get dimension space file is working on
          */
         Dimensions::Space::Id   space() const;

         /**
          * @brief Set the physical parameters of the simulation
          *
          * @param parameters Physical parameters
          * @param boundary Boundary flags
          */
         void setPhysical(const std::map<std::string,MHDFloat>& parameters, const std::map<std::string,std::size_t>& boundary);

         /**
          * @brief Set the mesh grid arrays
          *
          * @param mesh    Grid arrays of the mesh
          */
         void setMesh(const std::vector<Array>& mesh);

         /**
          * @brief Set the simulation time parameters
          *
          * @param time       Reached simulation time
          * @param timestep   Last timestep size
          */
         void setSimTime(const MHDFloat time, const MHDFloat timestep);

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
          * @brief Perform pre calculations
          */
         virtual void preCompute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Perform calculations
          */
         virtual void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Perform post calculations
          */
         virtual void postCompute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const;

      protected:
         /// Typedef for the scalar const iterator
         typedef std::map<std::size_t,Framework::Selector::VariantSharedScalarVariable>::const_iterator  scalar_iterator;

         /// Typedef for the vector const iterator
         typedef std::map<std::size_t,Framework::Selector::VariantSharedVectorVariable>::const_iterator  vector_iterator;

         /// Typedef for the scalar iterator range
         typedef std::pair<scalar_iterator,scalar_iterator>  scalar_iterator_range;

         /// Typedef for the scalar iterator range
         typedef std::pair<vector_iterator, vector_iterator>  vector_iterator_range;

         /**
          * @brief Get iterator range to scalars
          */
         scalar_iterator_range   scalarRange();

         /**
          * @brief Get iterator range to vectors
          */
         vector_iterator_range   vectorRange();

         /**
          * @brief Resolution information
          * 
          * @param spRes      Resolution information
          */
         SharedResolution mspRes;

         /**
          * @brief Physical parameters of the simulation
          */
         std::map<std::size_t,NonDimensional::SharedINumber> mPhysical;

         /**
          * @brief Boundary flags of the simulation
          */
         std::map<std::string,std::size_t> mBoundary;

         /**
          * @brief Storage for the mesh
          */
         std::vector<Array> mMesh;

         /**
          * @brief Time
          */
         MHDFloat mTime;

         /**
          * @brief Timestep
          */
         MHDFloat mTimestep;

      private:
         /**
          * @brief Set the resolution and use it for preliminary initialisation
          */
         void setResolution(SharedResolution spRes);

         /**
          * @brief The dimension space the file is working on
          */
         Dimensions::Space::Id mSpaceId;

         /**
          * @brief Storage for the names of the expecte variables
          */
         std::set<std::size_t>  mExpected;

         /**
          * @brief Storage for the scalars
          */
         std::map<std::size_t,Framework::Selector::VariantSharedScalarVariable>   mScalars;

         /**
          * @brief Storage for the vectors
          */
         std::map<std::size_t,Framework::Selector::VariantSharedVectorVariable>   mVectors;
   };

   /// Typedef for a smart reference counting pointer of a Variable HDF5 numbering writer
   typedef std::shared_ptr<IStatisticsAsciiWriter>   SharedIStatisticsAsciiWriter;

   inline bool IStatisticsAsciiWriter::isHeavy() const
   {
      return false;
   }

}
}
}

#endif // QUICC_IO_STATS_ISTATISTICSASCIIWRITER_HPP
