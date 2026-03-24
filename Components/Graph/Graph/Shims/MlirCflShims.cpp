#include <cassert>
#include <complex>
#include <iostream>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/Types.hpp"
#include "ViewOps/Reduction/Cuda/Reduction.hpp"
#include "Profiler/Interface.hpp"


using namespace QuICC::Graph;



#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a dot operator
/// gpu backend
/// @param obj pointer to operator implementation
/// @param pRet dot product
/// @param pU0 lhs vector
/// @param pU1 lhs vector
/// @param pU2 lhs vector
/// @param pV0 rhs vector
/// @param pV1 rhs vector
/// @param pV2 rhs vector
extern "C" void
_ciface_quiccir_cfl_magvel_f64_DCCSC3D_f64_DCCSC3D_gpu(
   void* obj, view3_t* ret, view3_t* pUval0, view3_t* pUval1, view3_t* pUval2, view3_t* pMval0, view3_t* pMval1, view3_t* pMval2)
{
    QuICC::Profiler::RegionFixture<3> fixTotal("MlirCflShims");

#ifndef NDEBUG
   std::cout << "_ciface_quiccir_cfl_magvel_f64_DCCSC3D_f64_DCCSC3D_gpu\n";
#endif
   assert(obj != nullptr);
   assert(pUval0 != nullptr);
   assert(pMval0 != nullptr);
   assert(QuICC::Cuda::isDeviceMemory(pUval0->data));
   assert(QuICC::Cuda::isDeviceMemory(pMval0->data));
   // op
   using namespace QuICC::Reduction::Cuda;
   using namespace QuICC::Reduction;
   using Tout = R_DCCSC3D_t;
   using Tin = R_DCCSC3D_t;
   using op_t = OpCfl<MagVelFunctor<double>, Tout, Tin, Tin, Tin, Tin, Tin, Tin>;
   ;
   // views
   using namespace QuICC::View;
   constexpr std::uint32_t rank = 3;
   ViewBase<std::uint32_t> pointers[rank];
   pointers[1] = ViewBase<std::uint32_t>(pUval0->pos, pUval0->posSize);
   ViewBase<std::uint32_t> indices[rank];
   indices[1] = ViewBase<std::uint32_t>(pUval0->coo, pUval0->cooSize);
   Tout viewRet(ret->data, ret->dataSize, ret->dims, pointers, indices);
   Tin viewUval0(pUval0->data, pUval0->dataSize, pUval0->dims, pointers, indices);
   Tin viewUval1(pUval1->data, pUval1->dataSize, pUval1->dims, pointers, indices);
   Tin viewUval2(pUval2->data, pUval2->dataSize, pUval2->dims, pointers, indices);
   Tin viewMval0(pMval0->data, pMval0->dataSize, pMval0->dims, pointers, indices);
   Tin viewMval1(pMval1->data, pMval1->dataSize, pMval1->dims, pointers, indices);
   Tin viewMval2(pMval2->data, pMval2->dataSize, pMval2->dims, pointers, indices);
   // call
   auto cl = reinterpret_cast<op_t*>(obj);
   cl->apply(viewRet, viewUval0, viewUval1, viewUval2, viewMval0, viewMval1, viewMval2);
};
#endif

/// @brief C Interface to MLIR for a dot operator
/// @param obj pointer to operator implementation
/// @param pRet dot product
/// @param pU0 lhs vector
/// @param pU1 lhs vector
/// @param pU2 lhs vector
/// @param pV0 rhs vector
/// @param pV1 rhs vector
/// @param pV2 rhs vector
extern "C" void
_ciface_quiccir_cfl_magvel_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D(
   void* obj, view3_t* ret, view3_t* pUval0, view3_t* pUval1, view3_t* pUval2, view3_t* pMval0, view3_t* pMval1, view3_t* pMval2)
{
#ifdef QUICC_HAS_CUDA_BACKEND
   if (QuICC::Cuda::isDeviceMemory(pUval0->data))
   {
      _ciface_quiccir_cfl_magvel_f64_DCCSC3D_f64_DCCSC3D_gpu(
         obj, ret, pUval0, pUval1, pUval2, pMval0, pMval1, pMval2);
   }
   else
   {
      
   }
#else

#endif
}
