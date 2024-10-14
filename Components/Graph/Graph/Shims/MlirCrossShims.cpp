#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "ViewOps/Pointwise/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"
#include "Graph/Types.hpp"


using namespace QuICC::Graph;

/// @brief C Interface to MLIR for a cross operator
/// @param obj pointer to operator implementation
/// @param pRet0 cross product
/// @param pRet0 cross product
/// @param pRet0 cross product
/// @param pU0 lhs vector
/// @param pU1 lhs vector
/// @param pU2 lhs vector
/// @param pV0 rhs vector
/// @param pV1 rhs vector
/// @param pV2 rhs vector
extern "C" void _ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_cpu(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet0,
    ViewDescriptor<double, std::uint32_t, 3>* pRet1,
    ViewDescriptor<double, std::uint32_t, 3>* pRet2,
    ViewDescriptor<double, std::uint32_t, 3>* pU0,
    ViewDescriptor<double, std::uint32_t, 3>* pU1,
    ViewDescriptor<double, std::uint32_t, 3>* pU2,
    ViewDescriptor<double, std::uint32_t, 3>* pV0,
    ViewDescriptor<double, std::uint32_t, 3>* pV1,
    ViewDescriptor<double, std::uint32_t, 3>* pV2
    )
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_cpu\n";
    #endif
    assert(obj != nullptr);
    assert(pRet0 != nullptr);
    assert(pRet1 != nullptr);
    assert(pRet2 != nullptr);
    assert(pU0 != nullptr);
    assert(pU1 != nullptr);
    assert(pU2 != nullptr);
    assert(pV0 != nullptr);
    assert(pV1 != nullptr);
    assert(pV2 != nullptr);
    // op
    using namespace QuICC::Pointwise::Cpu;
    using namespace QuICC::Pointwise;
    using T = QuICC::View::View<double, QuICC::View::DCCSC3D>;
    using op_t = Op<CrossCompFunctor<double>, T, T, T, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pU0->pos, pU0->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pU0->coo, pU0->cooSize);
    T viewU0(pU0->data, pU0->dataSize, pU0->dims, pointers, indices);
    T viewU1(pU1->data, pU1->dataSize, pU1->dims, pointers, indices);
    T viewU2(pU2->data, pU2->dataSize, pU2->dims, pointers, indices);
    T viewV0(pV0->data, pV0->dataSize, pV0->dims, pointers, indices);
    T viewV1(pV1->data, pV1->dataSize, pV1->dims, pointers, indices);
    T viewV2(pV2->data, pV2->dataSize, pV2->dims, pointers, indices);
    T viewRet0(pRet0->data, pRet0->dataSize, pRet0->dims, pointers, indices);
    T viewRet1(pRet1->data, pRet1->dataSize, pRet1->dims, pointers, indices);
    T viewRet2(pRet2->data, pRet2->dataSize, pRet2->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet0, viewU1, viewU2, viewV1, viewV2);
    cl->apply(viewRet1, viewU2, viewU0, viewV2, viewV0);
    cl->apply(viewRet2, viewU0, viewU1, viewV0, viewV1);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a cross operator
/// @param obj pointer to operator implementation
/// @param pRet0 cross product
/// @param pRet0 cross product
/// @param pRet0 cross product
/// @param pU0 lhs vector
/// @param pU1 lhs vector
/// @param pU2 lhs vector
/// @param pV0 rhs vector
/// @param pV1 rhs vector
/// @param pV2 rhs vector
extern "C" void _ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_gpu(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet0,
    ViewDescriptor<double, std::uint32_t, 3>* pRet1,
    ViewDescriptor<double, std::uint32_t, 3>* pRet2,
    ViewDescriptor<double, std::uint32_t, 3>* pU0,
    ViewDescriptor<double, std::uint32_t, 3>* pU1,
    ViewDescriptor<double, std::uint32_t, 3>* pU2,
    ViewDescriptor<double, std::uint32_t, 3>* pV0,
    ViewDescriptor<double, std::uint32_t, 3>* pV1,
    ViewDescriptor<double, std::uint32_t, 3>* pV2
    )
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_gpu\n";
    #endif
    assert(obj != nullptr);
    assert(pRet0 != nullptr);
    assert(pRet1 != nullptr);
    assert(pRet2 != nullptr);
    assert(pU0 != nullptr);
    assert(pU1 != nullptr);
    assert(pU2 != nullptr);
    assert(pV0 != nullptr);
    assert(pV1 != nullptr);
    assert(pV2 != nullptr);
    // op
    using namespace QuICC::Pointwise::Cuda;
    using namespace QuICC::Pointwise;
    using T = QuICC::View::View<double, QuICC::View::DCCSC3D>;
    using op_t = Op<CrossCompFunctor<double>, T, T, T, T, T>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    ViewBase<std::uint32_t> pointers[rank];
    pointers[1] = ViewBase<std::uint32_t>(pU0->pos, pU0->posSize);
    ViewBase<std::uint32_t> indices[rank];
    indices[1] = ViewBase<std::uint32_t>(pU0->coo, pU0->cooSize);
    T viewU0(pU0->data, pU0->dataSize, pU0->dims, pointers, indices);
    T viewU1(pU1->data, pU1->dataSize, pU1->dims, pointers, indices);
    T viewU2(pU2->data, pU2->dataSize, pU2->dims, pointers, indices);
    T viewV0(pV0->data, pV0->dataSize, pV0->dims, pointers, indices);
    T viewV1(pV1->data, pV1->dataSize, pV1->dims, pointers, indices);
    T viewV2(pV2->data, pV2->dataSize, pV2->dims, pointers, indices);
    T viewRet0(pRet0->data, pRet0->dataSize, pRet0->dims, pointers, indices);
    T viewRet1(pRet1->data, pRet1->dataSize, pRet1->dims, pointers, indices);
    T viewRet2(pRet2->data, pRet2->dataSize, pRet2->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewRet0, viewU1, viewU2, viewV1, viewV2);
    cl->apply(viewRet1, viewU2, viewU0, viewV2, viewV0);
    cl->apply(viewRet2, viewU0, viewU1, viewV0, viewV1);
};
#endif


/// @brief C Interface to MLIR for a cross operator
/// @param obj pointer to operator implementation
/// @param pRet0 cross product
/// @param pRet0 cross product
/// @param pRet0 cross product
/// @param pU0 lhs vector
/// @param pU1 lhs vector
/// @param pU2 lhs vector
/// @param pV0 rhs vector
/// @param pV1 rhs vector
/// @param pV2 rhs vector
extern "C" void _ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D(void* obj,
    ViewDescriptor<double, std::uint32_t, 3>* pRet0,
    ViewDescriptor<double, std::uint32_t, 3>* pRet1,
    ViewDescriptor<double, std::uint32_t, 3>* pRet2,
    ViewDescriptor<double, std::uint32_t, 3>* pU0,
    ViewDescriptor<double, std::uint32_t, 3>* pU1,
    ViewDescriptor<double, std::uint32_t, 3>* pU2,
    ViewDescriptor<double, std::uint32_t, 3>* pV0,
    ViewDescriptor<double, std::uint32_t, 3>* pV1,
    ViewDescriptor<double, std::uint32_t, 3>* pV2
    )
{
    #ifdef QUICC_HAS_CUDA_BACKEND
    assert(pRet0 != nullptr);
    assert(pRet1 != nullptr);
    assert(pRet2 != nullptr);
    assert(pU0 != nullptr);
    assert(pU1 != nullptr);
    assert(pU2 != nullptr);
    assert(pV0 != nullptr);
    assert(pV1 != nullptr);
    assert(pV2 != nullptr);
    if (!QuICC::Cuda::isDeviceMemory(pRet0->data))
    {
        assert(!QuICC::Cuda::isDeviceMemory(pRet1->data));
        assert(!QuICC::Cuda::isDeviceMemory(pRet2->data));
        assert(!QuICC::Cuda::isDeviceMemory(pU0->data));
        assert(!QuICC::Cuda::isDeviceMemory(pU1->data));
        assert(!QuICC::Cuda::isDeviceMemory(pU2->data));
        assert(!QuICC::Cuda::isDeviceMemory(pV0->data));
        assert(!QuICC::Cuda::isDeviceMemory(pV1->data));
        assert(!QuICC::Cuda::isDeviceMemory(pV2->data));
        _ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_cpu(obj, pRet0, pRet1, pRet2, pU0, pU1, pU2, pV0, pV1, pV2);
    }
    else
    {
        _ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_gpu(obj, pRet0, pRet1, pRet2, pU0, pU1, pU2, pV0, pV1, pV2);
    }
    #else
    _ciface_quiccir_cross_inertia_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_f64_DCCSC3D_cpu(obj, pRet0, pRet1, pRet2, pU0, pU1, pU2, pV0, pV1, pV2);
    #endif
}