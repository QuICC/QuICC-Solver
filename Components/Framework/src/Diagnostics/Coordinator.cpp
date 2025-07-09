/**
 * @file Coordinator.cpp
 * @brief Source of the diagnostic coordinator
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Diagnostics/Coordinator.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "Environment/MpiTypes.hpp"
#include "QuICC/Diagnostics/SphericalTorPolWrapper.hpp"
#include "QuICC/Diagnostics/CartesianTorPolWrapper.hpp"
#include "QuICC/Diagnostics/StreamVerticalWrapper.hpp"
#include "QuICC/PhysicalNames/Streamfunction.hpp"
#include "QuICC/PhysicalNames/VelocityZ.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Timestep/Constants.hpp"

namespace QuICC {

namespace Diagnostics {

#ifdef QUICC_MPI
// declaration for computing min with MPI
void mpi_cfl_min(void* a, void* b, int* len, MPI_Datatype* type);
#endif // QUICC_MPI

Coordinator::Coordinator() :
    mcMaxStep(Timestep::LIMIT_MAXSTEP),
    mcMinStep(Timestep::LIMIT_MINSTEP),
    mFixedStep(-1),
    mMaxError(-1.0),
    mCfl(2, 1),
    mStartTime(0.0),
    mStartTimestep(0.0)
{
   this->mCfl.setZero();
}

void Coordinator::addCfl(SharedICflWrapper spCfl)
{
   if(spCfl->isActive())
   {
      this->mCflOps.push_back(spCfl);
   }
}

void Coordinator::init(const std::vector<Array>& mesh,
   const std::map<std::size_t,
      Framework::Selector::VariantSharedScalarVariable>& scalars,
   const std::map<std::size_t,
      Framework::Selector::VariantSharedVectorVariable>& vectors,
   const Array& tstep,
   const std::map<std::size_t, NonDimensional::SharedINumber>& params)
{
   // Check for constant timestep setup
   if (tstep(1) > 0)
   {
      this->mFixedStep = std::max(this->mcMinStep, tstep(1));
      this->mFixedStep = std::min(this->mFixedStep, this->mcMaxStep);

      // Clear CFL operators
      this->mCflOps.clear();
   }
   else if(this->mCflOps.size() > 0)
   {
      for(auto cfl: this->mCflOps)
      {
         auto ids = cfl->fieldIds();
         for(auto id: ids)
         {
            if(vectors.count(id) > 0)
            {
               namespace S = SpatialScheme;
               std::shared_ptr<const S::ISpatialScheme> spScheme = std::visit(
                  [](auto&& p) { return p->dom(0).res().sim().spSpatialScheme(); },
                  vectors.find(id)->second);

               // Create a toroidal/poloidal spherical wrapper
               if (spScheme->formulation() == VectorFormulation::TORPOL)
               {
                  if (spScheme->has(SpatialScheme::Feature::ShellGeometry) || 
                        spScheme->has(SpatialScheme::Feature::SphereGeometry))
                  {
                     auto spField = std::make_shared<SphericalTorPolWrapper>(
                        vectors.find(id)->second);
                     cfl->setField(id, spField);
                  }
                  else if (spScheme->has(SpatialScheme::Feature::CartesianGeometry))
                  {
                     auto spField = std::make_shared<CartesianTorPolWrapper>(
                        vectors.find(id)->second);
                     cfl->setField(id, spField);
                  }
               }
            }
            else if(scalars.count(id) > 0)
            {
               namespace S = SpatialScheme;
               std::shared_ptr<const S::ISpatialScheme> spScheme = std::visit(
                  [](auto&& p) { return p->dom(0).res().sim().spSpatialScheme(); },
                  scalars.find(id)->second);

               if (spScheme->has(SpatialScheme::Feature::CartesianGeometry))
               {
                  if(id == PhysicalNames::Streamfunction::id() &&
                        scalars.count(PhysicalNames::VelocityZ::id() > 0))
                  {
                     auto spVelocity = std::make_shared<StreamVerticalWrapper>(
                        scalars.find(id)->second,
                        scalars.find(PhysicalNames::VelocityZ::id())->second);
                     cfl->setField(id, spVelocity);
                  }
               }
            }
         }
      }

      this->mFixedStep = tstep(1);
   }
   // Required wrapper is not implemented
   else
   {
      this->mFixedStep = tstep(1);

      // Clear CFL operators
      this->mCflOps.clear();
   }

   if (this->mCflOps.size() > 0)
   {
      int cols = 1;
      for(auto op: this->mCflOps)
      {
         assert(op);
         op->init(mesh);
         cols += op->initialCfl().cols();
      }

      this->mCfl.resize(2, cols);
      this->mCfl.setZero();
   }

   // Store configuration file start time
   this->mStartTime = tstep(0);

   // Store configuration file start step
   this->mStartTimestep = tstep(1);

   // Store error goal from configuration file (not enabled if fixed timestep is
   // used)
   if (tstep(2) > 0)
   {
      this->mMaxError = tstep(2);
   }
}

void Coordinator::initialCfl()
{
   // Used fixed timestep
   if (this->mFixedStep > 0 && this->mMaxError > 0)
   {
      this->mCfl(0, 0) = this->mcMinStep;
      this->mCfl(1, 0) = Timestep::MINSTEP_LOCATION;
   }
   else if (this->mFixedStep > 0)
   {
      this->mCfl(0, 0) = this->mFixedStep;
      this->mCfl(1, 0) = Timestep::FIXEDSTEP_LOCATION;
   }
   // Compute initial CFL condition
   else if (this->mCflOps.size() > 0)
   {
      // Compute CFL for initial state
      int col = 1;
      for(auto op: this->mCflOps)
      {
         assert(op);

         auto tmp = op->initialCfl();
         this->mCfl.block(0, col, 2, tmp.cols()) = tmp;
         col += tmp.cols();
      }
      this->updateCflMatrix(this->mCfl);

      if (this->mcMinStep < this->mCfl(0, 0))
      {
         this->mCfl(0, 0) = this->mcMinStep;
         this->mCfl(1, 0) = Timestep::MINSTEP_LOCATION;
      }
   }
}

void Coordinator::updateCfl()
{
   // Used fixed timestep
   if (this->mFixedStep > 0)
   {
      this->mCfl(0, 0) = this->mFixedStep;
      this->mCfl(1, 0) = Timestep::FIXEDSTEP_LOCATION;
   }
   // Compute CFL condition
   else if (this->mCflOps.size() > 0)
   {
      // Safety assert
      assert(this->mCflOps.size() > 0);

      int col = 1;
      for(auto op: this->mCflOps)
      {
         assert(op);

         auto tmp = op->cfl();
         this->mCfl.block(0, col, 2, tmp.cols()) = tmp;
         col += tmp.cols();
      }
      this->updateCflMatrix(this->mCfl);

      // Check for maximum timestep
      if (this->mcMaxStep < this->mCfl(0, 0))
      {
         this->mCfl(0, 0) = this->mcMaxStep;
         this->mCfl(1, 0) = Timestep::MAXSTEP_LOCATION;
      }

      if (-this->mFixedStep < this->mCfl(0, 0))
      {
         this->mCfl(0, 0) = -this->mFixedStep;
         this->mCfl(1, 0) = Timestep::FIXEDSTEP_LOCATION;
      }
   }
}

void Coordinator::synchronize()
{
//
// Start of MPI block
//
#ifdef QUICC_MPI

   if (this->mFixedStep <= 0 && this->mCflOps.size() > 0)
   {
      // Create MPI operation
      MPI_Op op;
      MPI_Op_create(mpi_cfl_min, true, &op);

      // Create MPI datatype
      MPI_Datatype ctype;
      MPI_Type_contiguous(2, Environment::MpiTypes::type<MHDFloat>(), &ctype);
      MPI_Type_commit(&ctype);

      // Reduce CFL on all CPUs to the global minimum
      MPI_Allreduce(MPI_IN_PLACE, this->mCfl.data(), this->mCfl.cols(), ctype,
         op, MPI_COMM_WORLD);
   }

//
// End of MPI block
//
#endif // QUICC_MPI
}

MHDFloat Coordinator::maxError() const
{
   return this->mMaxError;
}

const Matrix& Coordinator::cfl() const
{
   return this->mCfl;
}

MHDFloat Coordinator::startTime() const
{
   return this->mStartTime;
}

MHDFloat Coordinator::startTimestep() const
{
   return this->mStartTimestep;
}

void Coordinator::useStateTime(const MHDFloat time, const MHDFloat timestep)
{
   // Configuration requests use of state time
   if (this->mStartTime < 0)
   {
      this->mStartTime = time;
   }

   // Configuration requests use of state timestep
   if (this->mStartTimestep < 0)
   {
      this->mStartTimestep = timestep;
   }
}

void Coordinator::updateCflMatrix(Matrix& cfl) const
{
   int idx;
   cfl(0,0) = cfl.row(0).tail(cfl.cols()-1).minCoeff(&idx);

   if(cfl.rows() > 1)
   {
      cfl.col(0).tail(cfl.rows()-1) = cfl.col(idx+1).tail(cfl.rows()-1);
   }
}

#ifdef QUICC_MPI
void mpi_cfl_min(void* a, void* b, int* len, MPI_Datatype* type)
{
   MHDFloat* in = static_cast<MHDFloat*>(a);
   MHDFloat* inout = static_cast<MHDFloat*>(b);

   int cols = *len;
   int rows;
   MPI_Type_size(*type, &rows);
   rows /= sizeof(MHDFloat);

   for (int i = 0; i < cols; i++)
   {
      if (in[i * rows] < inout[i * rows])
      {
         inout[i * rows] = in[i * rows];
         for (int j = 1; j < rows; j++)
         {
            inout[i * rows + j] = in[i * rows + j];
         }
      }
   }
}
#endif // QUICC_MPI

} // namespace Diagnostics
} // namespace QuICC
