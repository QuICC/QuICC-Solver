// External includes
//
#include <iostream>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "ViewOps/Reduction/Cuda/Reduction.hpp"


namespace QuICC {
namespace Graph {

void MapOps::setCflMagVel(mlir::quiccir::CflMagVelOp op)
{
   // Get index from MLIR source
   std::uint64_t index = op.getImplptr().value();
   if (index >= _thisArr.size())
   {
      _thisArr.resize(index + 1, nullptr);
   }
   if (_thisArr[index] == nullptr)
   {
      using Tout = R_DCCSC3D_t;
      using Tin = R_DCCSC3D_t;
      double scaling = 1.0;

      if (_isCpu)
      {
      }
#ifdef QUICC_HAS_CUDA_BACKEND
      else
      {
         using namespace QuICC::Reduction::Cuda;
         using op_t = QuICC::Reduction::Cuda::OpCfl<QuICC::Reduction::Cuda::MagVelFunctor<double>, Tout, Tin, Tin, Tin, Tin, Tin, Tin>;
         _ops.push_back(std::make_unique<op_t>(QuICC::Reduction::Cuda::MagVelFunctor<double>(1.0)));
        auto* ptr =
            std::get<std::shared_ptr<NaryOp<Tout, Tin, Tin, Tin, Tin, Tin, Tin>>>(_ops.back()).get();
         assert(ptr != nullptr);
         _thisArr[index] = ptr;
      }
#endif
   }
   else
   {
#ifndef NDEBUG
      std::cout << "operator already allocated\n";
#endif
   }
}

} // namespace Graph
} // namespace QuICC
