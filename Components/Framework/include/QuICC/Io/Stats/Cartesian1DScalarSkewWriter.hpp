/** 
 * @file Cartesian1DScalarSkewWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 */

#ifndef QUICC_IO_STATS_CARTESIAN1DSCALARSKEWWRITER_HPP
#define QUICC_IO_STATS_CARTESIAN1DSCALARSKEWWRITER_HPP

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
#include "QuICC/Io/Stats/IStatisticsAsciiWriter.hpp"
#include "QuICC/Io/Stats/Cartesian1DScalarAvgWriter.hpp"
#include "QuICC/Io/Stats/Cartesian1DScalarRMSWriter.hpp"

namespace QuICC {

namespace Io {

namespace Stats {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) skew calculation for a scalar field
    */
   class Cartesian1DScalarSkewWriter: public IStatisticsAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarSkewWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarRMSWriter& RMS, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarSkewWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute skew for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Post Compute skew for scalar field
          */
         void postCompute(Transform::TransformCoordinatorType& coord);

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

          MHDFloat mArea;
         /**
          * @brief Storage for the scalar energy
          */
          SharedCartesian1DScalarAvgWriter mAvg;
          SharedCartesian1DScalarRMSWriter mRMS;
         Array mSkew;
   };

   inline bool Cartesian1DScalarSkewWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DScalarSkewWriter> SharedCartesian1DScalarSkewWriter;

}
}
}

#endif // QUICC_IO_STATS_CARTESIAN1DSCALARSKEWWRITER_HPP
