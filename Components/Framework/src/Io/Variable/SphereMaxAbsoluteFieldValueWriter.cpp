/**
 * @file SphereMaxAbsoluteFieldValueWriter.cpp
 * @brief Source of the implementation of the ASCII absolute field value in a sphere
 */

// System includes
//
#include <ctime>
#include <iomanip>
#include <stdexcept>

// Project includes
//
#include "QuICC/Io/Variable/SphereMaxAbsoluteFieldValueWriter.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "Types/Math.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Io/Variable/Tags/MaxAbsoluteFieldValue.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   SphereMaxAbsoluteFieldValueWriter::SphereMaxAbsoluteFieldValueWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + Tags::MaxAbsoluteFieldValue::BASENAME, Tags::MaxAbsoluteFieldValue::EXTENSION, prefix + Tags::MaxAbsoluteFieldValue::HEADER, type, Tags::MaxAbsoluteFieldValue::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mMaxAbsoluteFieldValue(std::numeric_limits<MHDFloat>::quiet_NaN())
   {
   }

   void SphereMaxAbsoluteFieldValueWriter::init()
   {
      this->mHasMOrdering = this->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123);

      IVariableAsciiWriter::init();
   }

   void SphereMaxAbsoluteFieldValueWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);

      this->mMaxAbsoluteFieldValue = std::visit(
            [&](auto&& p)
            {
               auto val = (
                  p->dom(0).phys().comp(FieldComponents::Physical::R).data().array().pow(2) +
                  p->dom(0).phys().comp(FieldComponents::Physical::THETA).data().array().pow(2) +
                  p->dom(0).phys().comp(FieldComponents::Physical::PHI).data().array().pow(2)
                  ).maxCoeff();
               return val;
            } ,vRange.first->second);

      this->mMaxAbsoluteFieldValue = std::sqrt(this->mMaxAbsoluteFieldValue);
   }

   void SphereMaxAbsoluteFieldValueWriter::writeContent()
   {
      // Create file
      this->preWrite();

      // Get the "global" max of the absolute field value from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mMaxAbsoluteFieldValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      using Tools::Formatter::ioFW;
      int ioPrec = 14;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::scientific;
         this->mFile << std::setprecision(ioPrec) << ioFW(ioPrec) << this->mTime << "\t" << ioFW(ioPrec) << this->mMaxAbsoluteFieldValue << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if absolute velocity is NaN
      if(std::isnan(this->mMaxAbsoluteFieldValue))
      {
         QuICCEnv().abort("Sphere absolute field value is NaN!");
      }
   }

}
}
}

