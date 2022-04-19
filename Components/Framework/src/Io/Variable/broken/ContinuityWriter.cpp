/**
 * @file ContinuityWriter.cpp
 * @brief Source of the implementation of the ASCII maximal continuity value writer
 */

// Configuration includes
//

// System includes
//
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/ContinuityWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/ContinuityTags.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ContinuityWriter::ContinuityWriter(std::string type)
      : IVariableAsciiWriter(ContinuityTags::BASENAME, ContinuityTags::EXTENSION, ContinuityTags::HEADER, type, ContinuityTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   ContinuityWriter::~ContinuityWriter()
   {
   }

   void ContinuityWriter::write()
   {
      // Create file
      this->preWrite();

      ContinuityWriter::vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      ContinuityWriter::vector_iterator  vIt = vRange.first;

      MHDFloat continuity;
      FieldComponents::Spectral::Id specONE = vIt->second->dom(0).res().sim().ss().spectral().ONE();
      FieldComponents::Spectral::Id specTWO = vIt->second->dom(0).res().sim().ss().spectral().TWO();
      FieldComponents::Physical::Id physONE = vIt->second->dom(0).res().sim().ss().physical().ONE();
      FieldComponents::Physical::Id physTWO = vIt->second->dom(0).res().sim().ss().physical().TWO();
      if(vIt->second->dom(0).res().sim().ss().dimension() == 3)
      {
         FieldComponents::Spectral::Id specTHREE = vIt->second->dom(0).res().sim().ss().spectral().THREE();
         FieldComponents::Physical::Id physTHREE = vIt->second->dom(0).res().sim().ss().physical().THREE();
         continuity = (vIt->second->dom(0).grad(specONE).comp(physONE).data() + vIt->second->dom(0).grad(specTWO).comp(physTWO).data() + vIt->second->dom(0).grad(specTHREE).comp(physTHREE).data()).array().abs().maxCoeff();
      } else
      {
         continuity = (vIt->second->dom(0).grad(specONE).comp(physONE).data() + vIt->second->dom(0).grad(specTWO).comp(physTWO).data()).array().abs().maxCoeff();
      }

      // Get the "global" velocity divergence from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &continuity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << continuity << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if continuity is NaN
      if(std::isnan(continuity))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Continuity is NaN!");
      }
   }

}
}
}
