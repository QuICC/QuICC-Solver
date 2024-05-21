#include <iostream>
#include <complex>
#include <cassert>

#include "Graph/Shims/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/Types.hpp"

using namespace QuICC::Graph;

namespace QuICC::Graph::details
{
    /// @brief map meta producer pointer to meta pointer
    std::unordered_map<void*, void*> metaStore;

    /// @brief map meta pointer to size for deallocation
    std::unordered_map<void*, std::uint32_t> metaSize;

    void cleanupMeta()
    {
        for (auto& m : metaSize)
        {
            if (m.second > 0)
            {
                // Check memory space
                bool isCpuMem = true;
                #ifdef QUICC_HAS_CUDA_BACKEND
                isCpuMem = !QuICC::Cuda::isDeviceMemory(m.first);
                #endif
                if (isCpuMem)
                {
                    ::operator delete(m.first, m.second,
                        static_cast<std::align_val_t>(sizeof(std::uint32_t)));
                }
                #ifdef QUICC_HAS_CUDA_BACKEND
                else
                {
                    cudaErrChk(cudaFree(m.first));
                }
                #endif
                #ifndef NDEBUG
                std::cout << "_ciface_quiccir_dealloc_meta, bytes: " << m.second << '\n';
                #endif
                m.second = 0;
            }
        }
    }

} // namespace details

/// @brief C Interface to MLIR for a transpose operator
/// cpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_C_DCCSC3DJIK_t_C_DCCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_C_DCCSC3DJIK_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_C_DCCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_C_DCCSC3D_t_C_DCCSC3DJIK_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_C_DCCSC3D_t_C_DCCSC3DJIK_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tout = C_DCCSC3D_t;
    using Tin = C_DCCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // Op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_S1CLCSC3D_t;
    using Tout = C_DCCSC3D_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_201_C_DCCSC3DJIK_t_C_S1CLCSC3DJIK_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_201_C_DCCSC3DJIK_t_C_S1CLCSC3DJIK_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // Op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tin = C_S1CLCSC3DJIK_t;
    using Tout = C_DCCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p201_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

/// @brief C Interface to MLIR for a transpose operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_C_S1CLCSC3D_t_C_DCCSC3D_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_C_S1CLCSC3D_t_C_DCCSC3D_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // Op
    using namespace QuICC::Transpose::Cpu;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3D_t;
    using Tout = C_S1CLCSC3D_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};

#ifdef QUICC_HAS_CUDA_BACKEND
/// @brief C Interface to MLIR for a transpose operator
/// gpu backend
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_transpose_120_C_S1CLCSC3DJIK_t_C_DCCSC3DJIK_t(void* obj, view3_cd_t* pOut, const view3_cd_t* pIn)
{
    #ifndef NDEBUG
    std::cout <<
        "_ciface_quiccir_transpose_120_C_S1CLCSC3DJIK_t_C_DCCSC3DJIK_t\n";
    #endif
    assert(obj != nullptr);
    assert(pIn != nullptr);
    assert(pOut != nullptr);
    assert(pIn->dataSize == pOut->dataSize);
    // Op
    using namespace QuICC::Transpose::Cuda;
    using namespace QuICC::Transpose;
    using Tin = C_DCCSC3DJIK_t;
    using Tout = C_S1CLCSC3DJIK_t;
    using op_t = Op<Tout, Tin, p120_t>;
    // views
    using namespace QuICC::View;
    constexpr std::uint32_t rank = 3;
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> pointers[rank];
    // not used for dense transpose, not setting up
    ViewBase<std::uint32_t> indices[rank];
    Tin viewIn(pIn->data, pIn->dataSize, pIn->dims, pointers, indices);
    Tout viewOut(pOut->data, pOut->dataSize, pOut->dims, pointers, indices);
    // call
    auto cl = reinterpret_cast<op_t*>(obj);
    cl->apply(viewOut, viewIn);
};
#endif

namespace QuICC::Graph::details
{
    void alloc_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
    {
        // Alloc meta for fully populated tensor
        // we will need to add a hashmap with counter
        // to know when we can deallocate the meta data
        pNewBuffer->posSize = pNewBuffer->dims[2]+1;
        std::size_t sizeByte = sizeof(std::uint32_t) * pNewBuffer->posSize;
        // Check if this metadata was already allocated;
        assert(pProdBuffer->pos != nullptr);
        if (details::metaStore[pProdBuffer->pos])
        {
            #ifndef NDEBUG
            std::cout << "meta pos buffer already allocated\n";
            #endif
            pNewBuffer->pos = reinterpret_cast<std::uint32_t*>(details::metaStore[pProdBuffer->pos]);
        }
        else
        {
            details::alloc_ptr(&pNewBuffer->pos, pNewBuffer->posSize, pProdBuffer->data);
            // Store meta
            details::metaStore[pProdBuffer->pos] = pNewBuffer->pos;
            // Register meta for dealloc
            details::metaSize[pNewBuffer->pos] = sizeByte;
            std::atexit(details::cleanupMeta);
            // Temporary auto move to host if needed
            QuICC::View::ViewBase<std::uint32_t> viewPos(pNewBuffer->pos, pNewBuffer->posSize);
            using namespace QuICC::Memory;
            tempOnHostMemorySpace converter(viewPos, TransferMode::write | TransferMode::block);
            // Populate meta for fully populated tensor
            viewPos[0] = 0;
            for (std::size_t i = 1; i < pNewBuffer->posSize; ++i) {
                viewPos[i] = viewPos[i-1]+pNewBuffer->dims[1];
            }
        }
        pNewBuffer->cooSize = pNewBuffer->dims[1]*pNewBuffer->dims[2];
        sizeByte = sizeof(std::uint32_t) * pNewBuffer->cooSize;
        assert(pProdBuffer->coo != nullptr);
        if (details::metaStore[pProdBuffer->coo])
        {
            #ifndef NDEBUG
            std::cout << "meta coo buffer already allocated\n";
            #endif
            pNewBuffer->coo = reinterpret_cast<std::uint32_t*>(details::metaStore[pProdBuffer->coo]);
        }
        else
        {
            details::alloc_ptr(&pNewBuffer->coo, pNewBuffer->cooSize, pProdBuffer->data);
            // Store meta
            details::metaStore[pProdBuffer->coo] = pNewBuffer->coo;
            // Register meta for dealloc
            details::metaSize[pNewBuffer->coo] = sizeByte;
            // Temporary auto move to host if needed
            QuICC::View::ViewBase<std::uint32_t> viewCoo(pNewBuffer->coo, pNewBuffer->cooSize);
            using namespace QuICC::Memory;
            tempOnHostMemorySpace converter(viewCoo, TransferMode::write | TransferMode::block);
            // Populate meta for fully populated tensor
            for (std::size_t i = 0; i < pNewBuffer->cooSize; ++i) {
                viewCoo[i] = i % pNewBuffer->dims[1];
            }
        }

        // Alloc buffer
        pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
        details::alloc_ptr(&pNewBuffer->data, pNewBuffer->dataSize, pProdBuffer->data);
    }



    void alloc_C_DCCSC3D_t_C_S1CLCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
    {
        // Alloc meta for fully populated tensor
        // we will need to add a hashmap with counter
        // to know when this was already allocated
        // and when we can deallocate the meta data
        auto meta = denseTransposePtrAndIdx<C_DCCSC3D_t, C_S1CLCSC3D_t>(
            {pNewBuffer->dims[0], pNewBuffer->dims[1], pNewBuffer->dims[2]});

        pNewBuffer->posSize = meta.ptr.size();
        std::size_t sizeByte = sizeof(std::uint32_t) * pNewBuffer->posSize;

        // Check if this metadata was already allocated;
        assert(pProdBuffer->pos != nullptr);
        if (details::metaStore[pProdBuffer->pos])
        {
            #ifndef NDEBUG
            std::cout << "meta pos buffer already allocated\n";
            #endif
            pNewBuffer->pos = reinterpret_cast<std::uint32_t*>(details::metaStore[pProdBuffer->pos]);
        }
        else
        {
            details::alloc_ptr(&pNewBuffer->pos, pNewBuffer->posSize, pProdBuffer->data);
            // Store meta
            details::metaStore[pProdBuffer->pos] = pNewBuffer->pos;
            // Register metaSize for dealloc
            details::metaSize[pNewBuffer->pos] = sizeByte;
            std::atexit(details::cleanupMeta);
            // Temporary auto move to host if needed
            QuICC::View::ViewBase<std::uint32_t> viewPos(pNewBuffer->pos, pNewBuffer->posSize);
            using namespace QuICC::Memory;
            tempOnHostMemorySpace converter(viewPos, TransferMode::write | TransferMode::block);
            // Populate meta for fully populated tensor
            for (std::size_t i = 0; i < pNewBuffer->posSize; ++i) {
                viewPos[i] = meta.ptr[i];
            }
        }
        pNewBuffer->cooSize = meta.idx.size();
        sizeByte = sizeof(std::uint32_t) * pNewBuffer->cooSize;
        assert(pProdBuffer->coo != nullptr);
        if (details::metaStore[pProdBuffer->coo])
        {
            #ifndef NDEBUG
            std::cout << "meta coo buffer already allocated\n";
            #endif
            pNewBuffer->coo = reinterpret_cast<std::uint32_t*>(details::metaStore[pProdBuffer->coo]);
        }
        else
        {
            details::alloc_ptr(&pNewBuffer->coo, pNewBuffer->cooSize, pProdBuffer->data);
            // Store meta
            details::metaStore[pProdBuffer->coo] = pNewBuffer->coo;
            // Register meta for dealloc
            details::metaSize[pNewBuffer->coo] = sizeByte;
            // Temporary auto move to host if needed
            QuICC::View::ViewBase<std::uint32_t> viewCoo(pNewBuffer->coo, pNewBuffer->cooSize);
            using namespace QuICC::Memory;
            tempOnHostMemorySpace converter(viewCoo, TransferMode::write | TransferMode::block);
            // Populate meta for fully populated tensor
            for (std::size_t i = 0; i < pNewBuffer->cooSize; ++i) {
                viewCoo[i] = meta.idx[i];
            }
        }

        // Alloc buffer
        pNewBuffer->dataSize = pNewBuffer->dims[0] * pNewBuffer->cooSize;
        details::alloc_ptr(&pNewBuffer->data, pNewBuffer->dataSize, pProdBuffer->data);
    }

    void alloc_C_S1CLCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
    {
        // Alloc meta for fully populated tensor
        // we will need to add a hashmap with counter
        // to know when this was already allocated
        // and when we can deallocate the meta data
        auto meta = denseTransposePtrAndIdx<C_S1CLCSC3D_t, C_DCCSC3D_t>(
            {pNewBuffer->dims[0], pNewBuffer->dims[1], pNewBuffer->dims[2]});

        pNewBuffer->posSize = meta.ptr.size();
        std::size_t sizeByte = sizeof(std::uint32_t) * pNewBuffer->posSize;

        // Check if this metadata was already allocated;
        assert(pProdBuffer->pos != nullptr);
        if (details::metaStore[pProdBuffer->pos])
        {
            #ifndef NDEBUG
            std::cout << "meta pos buffer already allocated\n";
            #endif
            pNewBuffer->pos = reinterpret_cast<std::uint32_t*>(details::metaStore[pProdBuffer->pos]);
        }
        else
        {
            details::alloc_ptr(&pNewBuffer->pos, pNewBuffer->posSize, pProdBuffer->data);
            // Store meta
            details::metaStore[pProdBuffer->pos] = pNewBuffer->pos;
            // Register metaSize for dealloc
            details::metaSize[pNewBuffer->pos] = sizeByte;
            std::atexit(details::cleanupMeta);
            // Temporary auto move to host if needed
            QuICC::View::ViewBase<std::uint32_t> viewPos(pNewBuffer->pos, pNewBuffer->posSize);
            using namespace QuICC::Memory;
            tempOnHostMemorySpace converter(viewPos, TransferMode::write | TransferMode::block);
            // Populate meta for fully populated tensor
            for (std::size_t i = 0; i < pNewBuffer->posSize; ++i) {
                viewPos[i] = meta.ptr[i];
            }
        }
        pNewBuffer->cooSize = meta.idx.size();
        sizeByte = sizeof(std::uint32_t) * pNewBuffer->cooSize;
        assert(pProdBuffer->coo != nullptr);
        if (details::metaStore[pProdBuffer->coo])
        {
            #ifndef NDEBUG
            std::cout << "meta coo buffer already allocated\n";
            #endif
            pNewBuffer->coo = reinterpret_cast<std::uint32_t*>(details::metaStore[pProdBuffer->coo]);
        }
        else
        {
            details::alloc_ptr(&pNewBuffer->coo, pNewBuffer->cooSize, pProdBuffer->data);
            // Store meta
            details::metaStore[pProdBuffer->coo] = pNewBuffer->coo;
            // Register meta for dealloc
            details::metaSize[pNewBuffer->coo] = sizeByte;
            // Temporary auto move to host if needed
            QuICC::View::ViewBase<std::uint32_t> viewCoo(pNewBuffer->coo, pNewBuffer->cooSize);
            using namespace QuICC::Memory;
            tempOnHostMemorySpace converter(viewCoo, TransferMode::write | TransferMode::block);
            // Populate meta for fully populated tensor
            for (std::size_t i = 0; i < pNewBuffer->cooSize; ++i) {
                viewCoo[i] = meta.idx[i];
            }
        }

        // Temporary auto move to host if needed
        QuICC::View::ViewBase<std::uint32_t> viewPos(pNewBuffer->pos, pNewBuffer->posSize);
        using namespace QuICC::Memory;
        tempOnHostMemorySpace converter(viewPos, TransferMode::read | TransferMode::block);
        // Alloc buffer
        std::size_t cumSliceSize = 0;
        for (std::size_t i = 0; i < pNewBuffer->posSize - 1; ++i) {
            auto width = viewPos[i+1] - viewPos[i];
            auto height = pNewBuffer->dims[0] - i;
            cumSliceSize += height * width;
        }
        pNewBuffer->dataSize = cumSliceSize;
        details::alloc_ptr(&pNewBuffer->data, pNewBuffer->dataSize, pProdBuffer->data);
    }

} // namespace details


/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[2]);

    details::alloc_C_DCCSC3D_t_C_DCCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_201_C_DCCSC3DJIK_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_201_C_DCCSC3DJIK_t_C_DCCSC3D_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[2]);

    details::alloc_C_DCCSC3D_t_C_DCCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_120_C_DCCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_120_C_DCCSC3D_t_C_DCCSC3D_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[2]);

    details::alloc_C_DCCSC3D_t_C_DCCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_120_C_DCCSC3D_t_C_DCCSC3DJIK_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_120_C_DCCSC3D_t_C_DCCSC3DJIK_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[2]);

    details::alloc_C_DCCSC3D_t_C_DCCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_201_C_DCCSC3D_t_C_S1CLCSC3D_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[2]);

    details::alloc_C_DCCSC3D_t_C_S1CLCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_201_C_DCCSC3DJIK_t_C_S1CLCSC3DJIK_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_201_C_DCCSC3DJIK_t_C_S1CLCSC3DJIK_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[2]);

    details::alloc_C_DCCSC3D_t_C_S1CLCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_120_C_S1CLCSC3D_t_C_DCCSC3D_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_120_C_S1CLCSC3D_t_C_DCCSC3D_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[2]);

    details::alloc_C_S1CLCSC3D_t_C_DCCSC3D_t(pNewBuffer, pProdBuffer);
};

/// @brief Temporary allocator
/// @param pNewBuffer
/// @param pProdBuffer
extern "C" void _ciface_quiccir_alloc_transpose_120_C_S1CLCSC3DJIK_t_C_DCCSC3DJIK_t(view3_cd_t* pNewBuffer, const view3_cd_t* pProdBuffer)
{
    #ifndef NDEBUG
    std::cout << "_ciface_quiccir_alloc_transpose_120_C_S1CLCSC3DJIK_t_C_DCCSC3DJIK_t\n";
    #endif
    // This operation allocates for the serial transpose operator
    // therefore it assumes dense 3D tensors

    // Check slice dimensions according to permutation
    // note that the mangling comes from the MLIR
    // convention where the layers are leftmost
    assert(pNewBuffer->dims[1] == pProdBuffer->dims[0]);
    assert(pNewBuffer->dims[2] == pProdBuffer->dims[1]);
    assert(pNewBuffer->dims[0] == pProdBuffer->dims[2]);

    details::alloc_C_S1CLCSC3D_t_C_DCCSC3D_t(pNewBuffer, pProdBuffer);
};
