/**
 * @file SphereMaxAbsoluteFieldValueWriter.hpp
 * @brief Implementation of the MaxAbsoluteFieldValue  in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHEREMAXABSOLUTEFIELDVALUEWRITER_HPP
#define QUICC_IO_VARIABLE_SPHEREMAXABSOLUTEFIELDVALUEWRITER_HPP

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
    * @brief Implementation of the absolute Velocity number in a sphere
    */
   class SphereMaxAbsoluteFieldValueWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereMaxAbsoluteFieldValueWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         ~SphereMaxAbsoluteFieldValueWriter() = default;

         /**
          * @brief Initialise the operator, transform and file
          */
         void init() final;

	      /**
          * @brief Compute max for vector field
          */
         void compute(Transform::TransformCoordinatorType& coord);

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
          * @brief MaxAbsoluteFieldValue number
          */
         MHDFloat mMaxAbsoluteFieldValue;

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef std::shared_ptr<SphereMaxAbsoluteFieldValueWriter> SharedSphereMaxAbsoluteFieldValueWriter;

   inline bool SphereMaxAbsoluteFieldValueWriter::isHeavy() const
   {
      return true;
   }

}
}
}

#endif // QUICC_IO_VARIABLE_SPHEREMAXABSOLUTEFIELDVALUEWRITER_HPP
