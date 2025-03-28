set(Tags "viewCpu_t")
if(QUICC_USE_CUFFT)
    list(APPEND Tags "viewGpu_t")
endif()
if(QUICC_USE_VKFFT)
    list(APPEND Tags "viewGpuVkFFT_t")
endif()

# Memres
set(viewCpu_tMemRes "Cpu::NewDelete")
set(viewGpu_tMemRes "Cuda::Malloc")
set(viewGpuVkFFT_tMemRes "Cuda::Malloc")
set(IntegratorMem "mem")

set(IntegratorOps "P")
set(IntegratorOpDirection "fwd_t")
set(IntegratorOpInOut "mods_t, phys_t")

# View Cpu - Integrator
set(IntegratorviewCpu_tViewStorage "// views
    mods_t mods({rOut.data(),
        static_cast<std::size_t>(rOut.size())},
        {{static_cast<std::uint32_t>(rOut.rows()),
            columns, layers}}, pointers, indices);
    phys_t phys({const_cast<std::remove_reference_t<decltype(in)>::Scalar*>(in.data()),
        static_cast<std::size_t>(in.size())},
            {{static_cast<std::uint32_t>(in.rows()),
            columns, layers}}, pointers, indices);")

# Indices/Pointers
set(IndicesPointers "// get 3D dimensions
    auto idBlocks = mspSetup->idBlocks();
    std::uint32_t activeLayers = idBlocks.rows();

    // here we guestimate based on highest index
    // since we don't how how the problem is setup globally
    std::uint32_t layers = idBlocks(activeLayers-1,0)+1;

    std::uint32_t columns = 0;
    for(std::uint32_t i = 0; i < activeLayers; ++i)
    {
        std::uint32_t maxCol = idBlocks(i, 1);
        columns = std::max(maxCol, columns);
    }

    std::uint32_t activeColumns = mPadded.cols();

    std::array<std::vector<std::uint32_t>, 3> pointers{};
    std::array<std::vector<std::uint32_t>, 3> indices{};
    auto psz = layers+1;
    auto isz = activeColumns;
    pointers[1].resize(psz);
    indices[1].resize(isz);

    // set pointers and indices
    pointers[1][0] = 0;
    std::uint32_t idk = 0;
    for(std::uint32_t ii = 0; ii < layers; ++ii)
    {
        std::uint32_t k = idBlocks(idk, 0);
        std::uint32_t numCol = 0;
        if (k == ii)
        {
            numCol = idBlocks(idk, 1);
            idk += 1;
        }

        pointers[1][ii+1] = pointers[1][ii] + numCol;

        std::uint32_t j = 0;
        for(std::uint32_t idj = pointers[1][ii]; idj < pointers[1][ii+1]; ++idj)
        {
            // making j index up
            indices[1][idj] = j;
            j += 1;
        }
    }"
)

# PreOp
set(IntegratorPreOp "void applyPreOperator(Matrix& out, const MatrixZ& in, const bool useReal) const final {};"
)

# PostOp
set(IntegratorPostOp 
  "void applyPostOperator(Matrix& rOut) const final {};
  void applyPostOperator(MatrixZ& rOut, const Matrix& in, const bool useReal) const final {};"
)

set(viewCpu_tSetDataViews
"// set pointers
    ViewBase<std::uint32_t> dataPointers[rank];
    dataPointers[1] = ViewBase<std::uint32_t>(mDataPointers.data(), mDataPointers.size());
    ViewBase<std::uint32_t> dataIndices[rank];
    dataIndices[1] = ViewBase<std::uint32_t>(mDataIndices.data(), mDataIndices.size());

    // set views
    mModsView = mods_t(mModsFlat.data(), mModsFlat.size(), modsDims.data(), dataPointers, dataIndices);
    mPhysView = phys_t(mPhysFlat.data(), mPhysFlat.size(), physDims.data(), dataPointers, dataIndices);"
)

# loop
set(Kind "Integrator")
set(OpDirection "${${Kind}OpDirection}")
set(Ops "${${Kind}Ops}")
set(Mem "${${Kind}Mem}")
set(OpInOut "${${Kind}OpInOut}")

set(PreOp "${${Kind}PreOp}")
set(PostOp "${${Kind}PostOp}")

foreach(Op IN LISTS Ops)
  foreach(Tag IN LISTS Tags)
    set(CuhHeader "${${Tag}CuhHeader}")
    set(MemBlock "${${Tag}MemBlock}")
    set(MemRes "${${Tag}MemRes}")
    set(IndicesPointers "${IndicesPointers}")
    set(SetDataViews "${${Tag}SetDataViews}")
    configure_file(
        "../OpV.hpp.in"
        "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Chebyshev/LinearMap/${Kind}/${Op}${Tag}.hpp.inc"
    )
    set_target_properties(${QUICC_CURRENT_COMPONENT_LIB}
        PROPERTIES
        PUBLIC_HEADER "${CMAKE_BINARY_DIR}/include/QuICC/Transform/Wrappers/Chebyshev/LinearMap/${Kind}/${Op}${Tag}.hpp.inc"
    )
    set(CppHeader "${CppHeader}")
    set(CuHeader "${${Tag}CuHeader}")
    set(Alloc "${${Tag}Alloc}")
    set(ViewStorage "${${Kind}${Tag}ViewStorage}")
    set(CopyBack "${${Kind}${Tag}CopyBack}")

    message(DEBUG "CppHeader: ${CppHeader}")
    message(DEBUG "PtrParam: ${PtrParam}")
    configure_file(
        "OpV.cpp.in"
        "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Chebyshev/LinearMap/${Kind}/${Op}${Tag}.cpp"
    )
    target_sources(${QUICC_CURRENT_COMPONENT_LIB} PRIVATE
      "${CMAKE_BINARY_DIR}/${QUICC_CURRENT_COMPONENT_DIR}/src/Wrappers/Chebyshev/LinearMap/${Kind}/${Op}${Tag}.cpp"
    )
  endforeach()
endforeach()
