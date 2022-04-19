/**
 * @file Cartesian1DScalarKurtWriter.cpp
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
#include "QuICC/Io/Stats/Cartesian1DScalarKurtWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Stats/KurtTags.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

   namespace Io {

   namespace Stats {

      Cartesian1DScalarKurtWriter::Cartesian1DScalarKurtWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarRMSWriter& RMS, const std::string& type)
         : IStatisticsAsciiWriter(prefix + KurtTags::BASENAME, KurtTags::EXTENSION, prefix + KurtTags::HEADER, type, KurtTags::VERSION, Dimensions::Space::SPECTRAL),mArea(-1), mAvg(Avg), mRMS(RMS), mKurt(-Array::Ones(2))
      {
      }

      Cartesian1DScalarKurtWriter::~Cartesian1DScalarKurtWriter()
      {
      }

      void Cartesian1DScalarKurtWriter::init()
      {
         IStatisticsAsciiWriter::init();

         if(QuICCEnv().allowsIO())
         {
            this->mFile << "# " << std::setprecision(14) << (1.0 + this->mMesh.at(0).transpose().reverse().array())/2.0 <<std::endl;
         }
         int dimZ = this->mspRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL);
         int dimX = this->mspRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::PHYSICAL);
         int dimY = this->mspRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::PHYSICAL);

         this->mArea = dimX*dimY;
         this->mKurt = Array::Zero(dimZ);

      }


      void Cartesian1DScalarKurtWriter::compute(Transform::TransformCoordinatorType& coord)
      {
         // calculate the transforms and calculate the RMS

         // Initialize the Skewness
         this->mKurt.setConstant(0.0);

         scalar_iterator_range sRange = this->scalarRange();
         assert(std::distance(sRange.first, sRange.second) == 1);

         // go through each vertical level and find the 0th mode (horizontal average)
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
            // Mean
            //this->mKurt(k_) = (this->rInVar.phys().slice(k) - this->mAvg->average()(k_)).array().pow(4).sum()/(this->mRMS->RMS()(k_)).pow(4);
            //this->mKurt(k_) = (sRange.first->second->dom(0).phys().slice(k).array() - mAvg->average()(k_)).array().pow(3).sum()/(mRMS->RMS    ()(k_)).pow(3);
            this->mKurt(k_) = (std::visit([](auto&& p)->auto&& {return p->dom(0).phys();}, sRange.first->second).slice(k).array() - mAvg->average()(k_)).array().pow(4).sum();
            this->mKurt(k_) = this->mKurt(k_)/this->mArea;

         }

      }

      void Cartesian1DScalarKurtWriter::postCompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

#ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mKurt.data(), this->mKurt.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif //QUICC_MPI

            this->mKurt = this->mKurt.array()/(mRMS->RMS().array()*mRMS->RMS().array()*mRMS->RMS().array()*mRMS->RMS().array());
      }

      void Cartesian1DScalarKurtWriter::writeContent()
      {
         // Create file
         this->preWrite();

         // Check if the workflow allows IO to be performed
         if(QuICCEnv().allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mTime << " \t" << this->mKurt.transpose() << std::endl;
         }

         // Close file
         this->postWrite();

         // Abort if kinetic energy is NaN
         if(std::isnan(this->mKurt.sum()) && this->mRMS->RMS().sum() != 0)
         {
            QuICCEnv().abort(99);

            throw std::logic_error("Kurtosis is NaN!");
         }
      }

   }
   }
}
