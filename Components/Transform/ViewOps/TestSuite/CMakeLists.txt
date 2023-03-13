include(BundleCatch2)

add_executable(ViewOpsTest "")

set(Tags "Cpu;Cuda")

set(CpuExtension "cpp")
set(CudaExtension "cu")

set(CudaHeader "#include <cuComplex.h>
#include \"Cuda/CudaUtil.hpp\"
#include \"Memory/Memory.hpp\"
#include \"Memory/Cuda/Malloc.hpp\"
")


Set(CudaView3D "// Get offsets of packed pointers
    auto psz = pointers[1].size();
    auto isz = indices[1].size();

    // device mem
    QuICC::Memory::Cuda::Malloc mem_res;
    QuICC::Memory::MemBlock<typename mods_t::ScalarType> modsIn_d(S, &mem_res);
    QuICC::Memory::MemBlock<typename mods_t::ScalarType> modsOut_d(S, &mem_res);
    QuICC::Memory::MemBlock<std::uint32_t> ptrAndIdx_d(psz + isz, &mem_res);

    {
    // Set device pointers (pointers and indices within the device block)
    constexpr size_t rank = 3;
    ViewBase<std::uint32_t> pointers_d[rank];
    pointers_d[1] = ViewBase<std::uint32_t>(ptrAndIdx_d.data(), psz);
    ViewBase<std::uint32_t> indices_d[rank];
    indices_d[1] = ViewBase<std::uint32_t>(ptrAndIdx_d.data() + psz, isz);

    // device view
    mods_t modsInV(modsIn_d.data(), modsIn.size(), dimensions.data(), pointers_d, indices_d);
    mods_t modsOutV(modsOut_d.data(), modsOut.size(), dimensions.data(), pointers_d, indices_d);

    // Copy to gpu
    cudaErrChk(cudaMemcpy(modsIn_d.data(), modsIn.data(),
        S*sizeof(typename mods_t::ScalarType), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpy(modsInV.pointers()[1].data(), pointers[1].data(),
        psz*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
    cudaErrChk(cudaMemcpy(modsInV.indices()[1].data(), indices[1].data(),
        isz*sizeof(std::uint32_t), cudaMemcpyHostToDevice));
"
)

set(CudaCopyBack3D "// copy back
    cudaErrChk(cudaMemcpy(modsOut.data(), modsOut_d.data(),
        S*sizeof(typename mods_t::ScalarType), cudaMemcpyDeviceToHost));
    }
   "
)

foreach(Tag IN LISTS Tags)
    set(Extension "${${Tag}Extension}")
    set(Header "${${Tag}Header}")
    set(View3D "${${Tag}View3D}")
    set(CopyBack3D "${${Tag}CopyBack3D}")
    configure_file(
        "MixedDiffTests.cpp.in"
        "${CMAKE_CURRENT_BINARY_DIR}/MixedDiffTest.${Extension}"
    )
    target_sources(ViewOpsTest PRIVATE
        "${CMAKE_CURRENT_BINARY_DIR}/MixedDiffTest.${Extension}"
    )
    configure_file(
        "ComplexDiffTests.cpp.in"
        "${CMAKE_CURRENT_BINARY_DIR}/ComplexDiffTest.${Extension}"
    )
    target_sources(ViewOpsTest PRIVATE
        "${CMAKE_CURRENT_BINARY_DIR}/ComplexDiffTest.${Extension}"
    )
endforeach()

target_sources(ViewOpsTest PRIVATE
"MainDiffTests.cpp"
)

target_link_libraries(ViewOpsTest
    ${QUICC_NAMESPACE}${QUICC_CURRENT_COMPONENT_LIB}
    Catch2::Catch2
    QuICC::Memory
)

add_test(
    NAME ViewOpsTest
    COMMAND ViewOpsTest
    WORKING_DIRECTORY
    "${QUICC_WORK_DIR}"
)
