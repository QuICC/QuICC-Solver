/**
 * @file StorageProfiler.cpp
 * @brief Source of the implementation of a storage profiler
 */

// Configuration includes
//

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "QuICC/Debug/StorageProfiler/StorageProfiler.hpp"

// Project includes
//

namespace QuICC {

namespace Debug {

   const unsigned int StorageProfiler::mActiveLvl = StorageProfiler::LVL1;

   std::map<StorageProfiler::StoragePoint, MHDFloat> StorageProfiler::mRequirements = std::map<StorageProfiler::StoragePoint, MHDFloat>();

   void StorageProfiler::init()
   {
   }

   bool StorageProfiler::isActive(StorageProfiler::StoragePoint point)
   {
      return ((static_cast<int>(point) % mActiveLvl) == 0);
   }

   MHDFloat StorageProfiler::requirement(StorageProfiler::StoragePoint point)
   {
      return mRequirements.at(point);
   }

   void StorageProfiler::update(StorageProfiler::StoragePoint point, MHDFloat memory)
   {
      if(isActive(point))
      {
         if(mRequirements.count(point) == 0)
         {
            mRequirements.insert(std::make_pair(point, 0.0));
         }

         // Increment required storage
         mRequirements.at(point) += memory;
      }
   }

   size_t StorageProfiler::size()
   {
      return mRequirements.size();
   }

   void StorageProfiler::getRequirements(Array& reqs)
   {
      reqs.resize(StorageProfiler::size());
      int i = 0;
      for(auto it = mRequirements.cbegin(); it != mRequirements.cend(); ++it)
      {
         reqs(i) = StorageProfiler::requirement(it->first);
         i++;
      }
   }

   void StorageProfiler::analyze(Array& reqs, Array& min, Array& max)
   {
      StorageProfiler::getRequirements(reqs);

      // Get the "global" requirements from MPI code
      #ifdef QUICC_MPI
         min.resize(reqs.size());
         max.resize(reqs.size());

         // Get the max values
         MPI_Allreduce(reqs.data(), max.data(), reqs.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         // Get the min values
         MPI_Allreduce(reqs.data(), min.data(), reqs.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

         // Get the mean values
         MPI_Allreduce(MPI_IN_PLACE, reqs.data(), reqs.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         int size;
         // Get the number of CPUs involved
         MPI_Comm_size(MPI_COMM_WORLD, &size);
         // Compute mean requirements per CPU
         reqs /= static_cast<MHDFloat>(size);
      #else
         min.resize(0);
         max.resize(0);
      #endif // QUICC_MPI
   }

   std::map<StorageProfiler::StoragePoint,std::string>  StorageProfiler::namePoints()
   {
      std::map<StorageProfiler::StoragePoint,std::string>  names;
      for(auto it = mRequirements.cbegin(); it != mRequirements.cend(); ++it)
      {
         names.insert(std::make_pair(it->first, StorageProfiler::pointName(it->first)));
      }

      return names;
   }

   std::string StorageProfiler::pointName(StorageProfiler::StoragePoint point)
   {
      switch(point)
      {
         //////////////////////////
         case VARIABLES:
            return "Variables";

         case VARIABLESSPEC:
            return "Spectral";
         case VARIABLESPHYS:
            return "Physical";

         //////////////////////////
         case TRANSFORMS:
            return "Transforms";

         case TRAALEGENDRE:
            return "Assoc. Legendre";
         case TRAMIXEDFFT:
            return "Mixed FFT";
         case TRACOMPLEXFFT:
            return "ComplexFFT 3D";
         case TRASHELLCHEBYSHEV:
            return "ComplexFFT 3D";
         case TRAANNULUSCHEBYSHEV:
            return "ComplexFFT 3D";
         case TRACARTESIANCHEBYSHEV:
            return "ComplexFFT 3D";
         case TRASPHEREWORLAND:
            return "Sphere Worland";
         case TRACYLINDERWORLAND:
            return "ComplexFFT 3D";

         //////////////////////////
         case TEMPORARIES:
            return "Temporaries";

         case TEMPTRA1D:
            return "Transform 1D";
         case TEMPTRA2D:
            return "Transform 2D";
         case TEMPTRA3D:
            return "Transform 3D";
         case TEMPSPECTRAL:
            return "Transform Spectral";

         //////////////////////////
         case EQUATIONS:
            return "Equations";

         case TRIVIALEQUATION:
            return "Trivial";
         case DIAGNOSTICEQUATION:
            return "Diagnostic";
         case PROGNOSTICEQUATION:
            return "Prognostic";

         //////////////////////////
         case MPI:
            return "MPI";

         case MPIBUFFERS:
            return "Buffers";
         case MPITYPES:
            return "Types";
         case MPICOMM:
            return "Comm";

         //////////////////////////
         case IO:
            return "IO";

         // Default output for unknown measure point
         default:
            return "Unknown measure point";
      }

   }

   void StorageProfiler::reset()
   {
      mRequirements.clear();
   }

}
}
