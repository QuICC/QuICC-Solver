/**
 * @file Builder.hpp
 * @brief Generic Worland parallALT operator builder
 */
#pragma once

// External includes
//
#include <memory>
#include <cstdint>

// Project includes
//
#include "Operator/Binary.hpp"
#include "Profiler/Interface.hpp"
#include "Types/Internal/Typedefs.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "ViewOps/Worland_parallALT/TypeTraits.hpp"

#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <nvrtc.h>
#include "parallALT.hpp"

namespace QuICC {
namespace Transform {
namespace Worland_parallALT {
using namespace QuICC::Operator;
using type = QuICC::Memory::Cuda::Malloc;
/// @brief Derived classes implement the differentiation in modal space
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Order of differentiation
/// @tparam Direction Fft direction tag
/// @tparam Treatment special treatment mask, typically of mode zero or dealiasing
template<class Tout, class Tin, std::int64_t Direction, std::int64_t Type>
class ParallaltOp : public UnaryBaseOp<ParallaltOp<Tout, Tin, Direction, Type>, Tout, Tin> {
public:
    /// @brief Default constructor
   ParallaltOp(std::shared_ptr<Memory::memory_resource> mem) : _mem(mem) {
   };
    /// @brief dtor
    ~ParallaltOp()
    {
       deleteParallALT(&VkGPU, &appContainer);
    };

    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void initImpl(Tout& out, const Tin& in){
        std::uint32_t Ntheta = (Direction) ? out.dims()[0] : in.dims()[0];//igrid.size();
        
        //std::uint32_t nLayers = static_cast<std::uint32_t>(this->mspSetup->slowSize());

        ///\todo this should be the full matrix size
        //std::uint32_t M = out.pointers()->size() - 1;

        PfSolve::PfSolveResult resPfSolve = PfSolve::PFSOLVE_SUCCESS;
        config = {};
	    config.Ntheta = Ntheta;//this->mspSetup->bwdSize();
        config.radialTransform = Type;
        config.useGraphs = 1;
	    config.testMerge = 0;
	    config.testAccuracy = 1;
	    config.doALTOnly = 1;
	    config.useMatMulConnection = 1;
	    config.use_tc = 2;
	    config.mergeType = 1;
	    config.numMergedIterMatMul = 8;
	    config.numMergedIterMax = 1;
	    config.numMergedIterMin = 1;
	    config.disableCaching = 1;
	    config.fixAccuracy = 1;
	    config.numRadialBatches = 1;
	    config.profile_iter = 1;
	    config.profile_iter_combined = 1;
        config.specifyBuffersAtLaunch = 1;
	    config.WMMA_M = 8;
	    config.WMMA_N = 8;
	    config.WMMA_K = 4;
	    appContainer = {};
        config.projector = Direction;

         std::uint32_t* temp_pointers =
           (std::uint32_t*)calloc(out.pointers()[1].size(), sizeof(std::uint32_t));
        cudaErrChk(cudaMemcpyAsync(temp_pointers, out.pointers()[1].data(),
        out.pointers()[1].size() * sizeof(std::uint32_t), cudaMemcpyDeviceToHost));
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
           int numRHS = temp_pointers[i] - temp_pointers[i-1];
           if ((((i - 1) % 2) == 0) && (numRHS!=0))
            {
                config.num_m_even++;
            }
            if ((((i-1)%2) == 1 ) && (numRHS!=0))
            {
                config.num_m_odd++;
            }
        }

        int start = 0;
        int iter = 0;
        int* m_even = (int*)calloc(config.num_m_even, sizeof(int));
        int* m_even_endBatch = (int*)calloc(config.num_m_even, sizeof(int));
        int* m_odd = (int*)calloc(config.num_m_odd, sizeof(int));
        int* m_odd_endBatch = (int*)calloc(config.num_m_odd, sizeof(int));
        config.m_even_list = m_even;
        config.m_even_endBatch = m_even_endBatch;
        config.m_odd_list = m_odd;
        config.m_odd_endBatch = m_odd_endBatch;
        
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
            int numRHS = temp_pointers[i] - temp_pointers[i-1];
            if ((((i - 1) % 2) == 0) && (numRHS!=0))
            {
                start += numRHS;
                m_even[iter] = i-1;
                m_even_endBatch[iter] = 2*start;
                printf("%d %d %d \n", m_even[iter], m_even_endBatch[iter], iter);
                iter++;
            }
        }
    
        //current_nCols = 0;
        start = 0;
        iter = 0;
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
            int numRHS = temp_pointers[i] - temp_pointers[i-1];
            if ((((i - 1) % 2) == 1) && (numRHS!=0))
            {
                start += numRHS;
                m_odd[iter] = i-1;
                m_odd_endBatch[iter] = 2*start;
                printf("%d %d %d \n", m_odd[iter], m_odd_endBatch[iter], iter);
                iter++;
            }

        }
        free(temp_pointers);
        config.M = ((m_even[config.num_m_even - 1]) / 2 + 1) * 2;// M;
        config.L = (Direction) ? in.dims()[0] : out.dims()[0];// 3 * M / 2;

	    //appContainer.input_buffer_S = (double*)in.data();
        //appContainer.buffer_S = (double*)out.data();
       initializeParallALT(&VkGPU, config, &appContainer);
    };
    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void applyImpl(Tout& out, const Tin& in){
        if (appContainer.config.Ntheta == 0)
        {
            initImpl(out, in);
        }
      Profiler::RegionFixture<5> fix("ParallaltOp::applyImpl");

      assert(out.size() == in.size());
      //assert(out.dims()[0] == in.dims()[0]);
      assert(out.dims()[1] == in.dims()[1]);
      assert(out.dims()[2] == in.dims()[2]);
      assert(QuICC::Cuda::isDeviceMemory(out.data()));
      assert(QuICC::Cuda::isDeviceMemory(in.data()));
      parallALT_launchParams launchParams;
      launchParams.input_buffer_S = (double*)in.data();
      launchParams.buffer_S = (double*)out.data();

      launchApp_parallALT(&appContainer, &launchParams);
      };

 private:

    /**
    * @brief parallALT configurations
    */
    mutable parallALT_configuration config = {};

    /**
    * @brief parallALT app pointers
    */
    mutable parallALT_app appContainer = {};

    /**
    * @brief parallALT app pointers
    */
    mutable PfSolve::VkGPU VkGPU = {};

    /// @brief memory resource
    /// needs shared ptr for memory pools
    /// note, this must call the dtor last
    /// otherwise we cannot dealloc data
    /// \todo consider removing shared ptr and using singleton
    std::shared_ptr<Memory::memory_resource> _mem;

    /// @brief Give access to base class
    //friend BinaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin, ScaleType>;

};

} // namespace Worland_parallALT
} // namespace Transform
} // namespace QuICC
