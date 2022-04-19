/**
 * @file Cartesian1DNusseltZWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer through the Z boundary
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
#include "QuICC/Io/Variable/Cartesian1DNusseltZWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/NusseltTags.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DNusseltZWriter::Cartesian1DNusseltZWriter(std::string type)
      : IVariableAsciiWriter(Tags::Nusselt::BASENAME, Tags::Nusselt::EXTENSION, Tags::Nusselt::HEADER, type, Tags::Nusselt::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   Cartesian1DNusseltZWriter::~Cartesian1DNusseltZWriter()
   {
   }

   void Cartesian1DNusseltZWriter::write()
   {
      // Create file
      this->preWrite();

      Cartesian1DNusseltZWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      Cartesian1DNusseltZWriter::scalar_iterator sit = sRange.first;

      ArrayI mode = sit->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->mode(0);
      MHDFloat nusselt = 0.0;
      if(mode(2) == 0 && mode(3) == 0)
      {
         // Create boundary operator
         Array bc = 2.0*Array::Ones(this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
         bc(0) = bc(0)/2.0;

         // Compute Nusselt number
         nusselt = -bc.dot(sit->second->dom(0).perturbation().profile(0,0).real());
      }

      // Get the "global" Nusselt number from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &nusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << nusselt << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if Nusselt number is NaN
      if(std::isnan(nusselt))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Nusselt number is NaN!");
      }
   }

}
}
}
