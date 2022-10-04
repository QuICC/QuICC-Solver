/**
 * @file Cartesian1DScalarAvgWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) statistics RMS calculation for scalar field
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
#include "QuICC/Io/Stats/Cartesian1DScalarAvgWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Stats/AvgTags.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Transform/Backward/P.hpp"

namespace QuICC {

   namespace Io {

   namespace Stats {

      Cartesian1DScalarAvgWriter::Cartesian1DScalarAvgWriter(const std::string& prefix, const std::string& type)
         : IStatisticsAsciiWriter(prefix + AvgTags::BASENAME, AvgTags::EXTENSION, prefix + AvgTags::HEADER, type, AvgTags::VERSION, Dimensions::Space::SPECTRAL), mAvg(-Array::Ones(2))
      {
      }

      Cartesian1DScalarAvgWriter::~Cartesian1DScalarAvgWriter()
      {
      }

      void Cartesian1DScalarAvgWriter::init()
      {
         IStatisticsAsciiWriter::init();

         if(QuICCEnv().allowsIO())
         {
            this->mFile << "# " << std::setprecision(14) <<  (1.0 + this->mMesh.at(0).transpose().reverse().array())/2.0 << std::endl;
         }
      }

      void Cartesian1DScalarAvgWriter::preCompute(Transform::TransformCoordinatorType& coord)
      {
         // Dealias variable data
         scalar_iterator_range sRange = this->scalarRange();
         assert(std::distance(sRange.first, sRange.second) == 1);

         if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0 && this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0) == 0)
         {
            std::visit(
                  [&](auto&& p)
                  {
                     coord.communicator().transferForward(Dimensions::Transform::SPECTRAL, p->rDom(0).rTotal(), false);
                  }, sRange.first->second);

            // Recover dealiased BWD data
            auto pInVar = coord.ss().bwdPtr(Dimensions::Transform::TRA1D);
            coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd(pInVar);

            // Get FWD storage
            auto pOutVar = coord.ss().fwdPtr(Dimensions::Transform::TRA1D);
            coord.communicator().storage(Dimensions::Transform::TRA1D).provideFwd(pOutVar);

            // Compute projection transform for first dimension
            std::visit([&](auto&& pOut, auto&& pIn){coord.transform1D().backward(pOut->rData(), pIn->data(), Transform::Backward::P::id());}, pOutVar, pInVar);

            // set mAvg to be the 0th mode of the field
            this->mAvg = std::visit([](auto&& p)->Array{return p->slice(0).col(0).real();}, pOutVar);

            // Free BWD storage
            coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(pInVar);

            // Free FWD storage
            coord.communicator().storage(Dimensions::Transform::TRA1D).freeFwd(pOutVar);
         } else
         {
            this->mAvg = Array::Zero(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>());
         }
      }

      const Array& Cartesian1DScalarAvgWriter::average() const
      {
         return this->mAvg;
      }

      void Cartesian1DScalarAvgWriter::compute(Transform::TransformCoordinatorType& coord)
      {
         // nothing here for this one
      }

      void Cartesian1DScalarAvgWriter::postCompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

         #ifdef QUICC_MPI
            MPI_Allreduce(MPI_IN_PLACE, this->mAvg.data(), this->mAvg.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif //QUICC_MPI
      }

      void Cartesian1DScalarAvgWriter::writeContent()
      {
         // Create file
         this->preWrite();

         // Check if the workflow allows IO to be performed
         if(QuICCEnv().allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mAvg.transpose() << std::endl;
         }

         // Close file
         this->postWrite();

         // Abort if kinetic energy is NaN
         if(std::isnan(this->mAvg.sum()))
         {
            QuICCEnv().abort("Horizontal Avg is NaN!");
         }
      }

   }
   }
}
