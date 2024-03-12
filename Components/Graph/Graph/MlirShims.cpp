#include <iostream>
#include <Quiccir-c/Utils.hpp>

/// @brief C Interface to MLIR for a jw prj operator
/// @param op
/// @param uval
/// @param umod
extern "C" void _ciface_quiccir_jw_prj_layoutUval_layoutUmod(void* obj, view3_t* pUval, view3_t* pUmod)
{
//   auto cl = reinterpret_cast<JWOp*>(obj);
//   cl->apply(pUval, pUmod);
    std::cout << "_ciface_quiccir_jw_prj_layoutUval_layoutUmod\n";
};

/// @brief C Interface to MLIR for a jw int operator
/// @param op
/// @param umod
/// @param uval
extern "C" void _ciface_quiccir_jw_int_layoutUmod_layoutUval(void* obj, view3_t* pUmod, view3_t* pUval)
{
//   auto cl = reinterpret_cast<JWOp*>(obj);
//   cl->apply(pUmod, pUval);
    std::cout << "_ciface_quiccir_jw_int_layoutUmod_layoutUval\n";
};

/// @brief C Interface to MLIR for an allocator
/// @param uval
/// @param op
/// @param umod
extern "C" void _ciface_quiccir_alloc_jw_prj_layoutUval_layoutUmod(view3_t* pNewBuffer, view3_t* pProdBuffer)
{
  /// proto is dense
//   pNewBuffer->coo = nullptr;
//   pNewBuffer->cooSize = 0;
//   pNewBuffer->pos = nullptr;
//   pNewBuffer->posSize = 0;
//   std::size_t size = pNewBuffer->dims[0] * pNewBuffer->dims[1] * pNewBuffer->dims[2];
//   pNewBuffer->dataSize = size;
//   std::size_t sizeByte = sizeof(double) * pNewBuffer->dataSize;
//   pNewBuffer->data = reinterpret_cast<double*>(::operator new(sizeByte, static_cast<std::align_val_t>(sizeof(double))));
    std::cout << "_ciface_quiccir_alloc_jw_prj_layoutUval_layoutUmod\n";
};

/// @brief C Interface to MLIR for a  deallocator
/// @param umod
extern "C" void _ciface_quiccir_dealloc_layoutUval(view3_t* pBuffer)
{
//   /// proto is dense
//   // pNewBuffer->coo = nullptr;
//   // pNewBuffer->cooSize = 0;
//   // pNewBuffer->pos = nullptr;
//   // pNewBuffer->posSize = 0;
//   std::size_t sizeByte = sizeof(double) * pBuffer->dataSize;
//   ::operator delete(pBuffer->data, sizeByte, static_cast<std::align_val_t>(sizeof(double)));
//   pBuffer->dataSize = 0;
    std::cout << "_ciface_quiccir_dealloc_layoutUval\n";
};