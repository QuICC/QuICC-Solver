/** 
 * @file Cartesian1DScalarAvgWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 */

#ifndef QUICC_IO_STATS_CARTESIAN1DSCALARAVGWRITER_HPP
#define QUICC_IO_STATS_CARTESIAN1DSCALARAVGWRITER_HPP

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

namespace QuICC {

namespace Io {

namespace Stats {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) horizontal average calculation for a scalar field
    */
   class Cartesian1DScalarAvgWriter: public IStatisticsAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarAvgWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarAvgWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Pre computation stage for average of scalar field
          */
         void preCompute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Compute average of scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Post computation stage for average of scalar field
          */
         void postCompute(Transform::TransformCoordinatorType& coord);

         /**
          * To share average with other stats
          */
         const Array& average() const;

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
          * @brief Storage for the scalar energy
          */
         Array mAvg;
   };

   inline bool Cartesian1DScalarAvgWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<Cartesian1DScalarAvgWriter> SharedCartesian1DScalarAvgWriter;

}
}
}
   
#endif // QUICC_IO_STATS_CARTESIAN1DSCALARAVGWRITER_HPP
