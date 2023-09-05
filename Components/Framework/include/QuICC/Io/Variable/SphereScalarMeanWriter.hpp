/**
 * @file SphereScalarMeanWriter.hpp
 * @brief Implementation of the mean value of a scalar field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERESCALARMEANWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERESCALARMEANWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   /**
    * @brief Implementation of the mean value of a scalar field in a sphere
    */
   class SphereScalarMeanWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereScalarMeanWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         ~SphereScalarMeanWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         void init() final;

         /**
          * @brief Requires heavy calculation?
          */
         bool isHeavy() const final;

      protected:
         /**
          * @brief Write State to file
          */
         void writeContent() final;

         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

      private:
         /**
          * @brief Mean value
          */
         MHDFloat mMean;

         /*
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Background state
          */
         ArrayZ mBg;

         /**
          * @brief Projector to extract mean
          */
         Matrix mProjector;
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<SphereScalarMeanWriter> SharedSphereScalarMeanWriter;

   inline bool SphereScalarMeanWriter::isHeavy() const
   {
      return false;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_SPHERESCALARMEANWRITER_HPP
