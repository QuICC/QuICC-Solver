/**
 * @file Builder.hpp
 * @brief Generic ALegendre parallALT operator builder
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
#include "ViewOps/ALegendre_parallALT/TypeTraits.hpp"

#include <cuComplex.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <nvrtc.h>
#include "parallALT.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre_parallALT {
using namespace QuICC::Operator;
using type = QuICC::Memory::Cuda::Malloc;
/// @brief Derived classes implement the differentiation in modal space
/// @tparam Tout differentiated modes type
/// @tparam Tin input modes type
/// @tparam Order of differentiation
/// @tparam Direction Fft direction tag
/// @tparam Treatment special treatment mask, typically of mode zero or dealiasing
template<class Tout, class Tin>
class ParallaltOp : public UnaryBaseOp<ParallaltOp<Tout, Tin>, Tout, Tin> {
public:
    /// @brief Default constructor
   ParallaltOp(std::shared_ptr<Memory::memory_resource> mem) : _mem(mem) {
   };
    /// @brief dtor
    ~ParallaltOp()
    {
        if (appContainer.config.Ntheta != 0)
            deleteParallALT(&VkGPU, &appContainer);

       if (temp_buffer != 0)
        {
           cudaFree(temp_buffer);
           temp_buffer = 0;
        };
        if (temp_buffer2 != 0)
        {
           cudaFree(temp_buffer2);
           temp_buffer2 = 0;
        };
        if (temp_buffer3 != 0)
        {
           cudaFree(temp_buffer3);
           temp_buffer3 = 0;
        };
        if (temp_buffer4 != 0)
        {
           cudaFree(temp_buffer4);
           temp_buffer4 = 0;
        };
    };

    /// @brief Action implementation
    /// @param out differentiatied modes
    /// @param in input modes
    void initImpl(Tout& out, const Tin& in){
        std::uint32_t Ntheta = (transformDirection) ? out.dims()[0] : in.dims()[0];//igrid.size();
        
        //std::uint32_t nLayers = static_cast<std::uint32_t>(this->mspSetup->slowSize());

        ///\todo this should be the full matrix size
        //std::uint32_t M = out.pointers()->size() - 1;

        PfSolve::PfSolveResult resPfSolve = PfSolve::PFSOLVE_SUCCESS;
        config = {};
	    config.Ntheta = Ntheta;//this->mspSetup->bwdSize();
        config.radialTransform = 0;
        config.vectorSHT = (transformType/10)%10;
        if ((config.vectorSHT == 1) && ((transformType/100) == 0))
        {
           config.vectorSHT = 3;
        }
        config.shTransform = (transformType % 10) ? BLOCK_MUL_LL : 0;
        if ((!transformDirection) && (config.vectorSHT)) config.shTransform = 0;//BLOCK_DIV_LL;
        config.useGraphs = 0;
        config.initializeGraph = 0;
        config.useUberKernel = 1;
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
        config.convertPackedStrided = 2;
	    config.WMMA_M = 8;
	    config.WMMA_N = 8;
	    config.WMMA_K = 4;
	    appContainer = {};
        config.projector = transformDirection;
        
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

        int* m_list = (int*)calloc(config.num_m_even + config.num_m_odd, sizeof(int));
        int* m_endBatch = (int*)calloc(config.num_m_even+config.num_m_odd, sizeof(int));
	   
        for (std::uint32_t i = 1; i < out.pointers()[1].size(); ++i)
        {
            int numRHS = temp_pointers[i] - temp_pointers[i-1];
           if (numRHS != 0)
           {
              start += numRHS;
              m_endBatch[iter] = 2 * start;
              m_list[iter] = i - 1;
              printf("%d %d %d \n", m_list[iter], m_endBatch[iter], iter);
              iter++;
           }
        }
        config.m_list = m_list;
	    config.m_endBatch = m_endBatch;

        config.m_even_list = m_even;
        config.m_even_endBatch = m_even_endBatch;
        config.m_odd_list = m_odd;
        config.m_odd_endBatch = m_odd_endBatch;
        start = 0;
        iter = 0;
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
        config.M = (m_even[config.num_m_even - 1] > m_odd[config.num_m_odd - 1]) ? ((m_even[config.num_m_even - 1]) / 2 + 1) * 2 : ((m_odd[config.num_m_odd - 1]) / 2 + 1) * 2;// M;
        if (config.projector) {
            config.L = in.dims()[0];
		    config.inputBufferStride = in.dims()[0];
		    config.outputBufferStride = Ntheta;
        }
        else
        {
            config.L = Ntheta;
            config.inputBufferStride = Ntheta;
		    config.outputBufferStride = out.dims()[0];
        }
	    //appContainer.input_buffer_S = (double*)in.data();
        //appContainer.buffer_S = (double*)out.data();
       initializeParallALT(&VkGPU, config, &appContainer);
        free(m_even);
        free(m_even_endBatch);
        free(m_odd);
        free(m_odd_endBatch);
        free(m_list);
        free(m_endBatch);
    };

    void applyImpl(Tout& out, const Tin& in){
       //printf("%d\n", appContainer.config.Ntheta);
        if (appContainer.config.Ntheta == 0)
        {
            initImpl(out, in);
        }
      Profiler::RegionFixture<5> fix("ParallaltOp::applyImpl");

      //assert(out.size() == in.size());
      //assert(out.dims()[0] == in.dims()[0]);
      assert(out.dims()[1] == in.dims()[1]);
      assert(out.dims()[2] == in.dims()[2]);

      assert(QuICC::Cuda::isDeviceMemory(out.data()));
      assert(QuICC::Cuda::isDeviceMemory(in.data()));

      if (temp_buffer == 0)
      {
         cudaMalloc((void**)&temp_buffer,
            appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
               sizeof(double));
      }
      if (temp_buffer2 == 0)
      {
         cudaMalloc((void**)&temp_buffer2,
            appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
               sizeof(double));
      }
      
      
        launchParams.input_buffer_S = (double*)in.data();
        launchParams.temp_buffer_S = temp_buffer;
        launchParams.buffer_S = temp_buffer2;
        launchParams.output_buffer_S = (double*)out.data();
        launchApp_parallALT(&appContainer, &launchParams);
     
     /* double* xx =
                (double*)calloc(2 * appContainer.config.Ntheta *
                                         (appContainer.config.sizeEvenBlock +
                                            appContainer.config.sizeOddBlock),
            sizeof(double));
        
      cudaDeviceSynchronize();
          cudaMemcpy(xx, launchParams.input_buffer_S,
            68 *
               2*300 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (appContainer.config.projector)
         {
            for (int j = 0; j < 100; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
         printf("in_scalar\n");
        
      cudaMemcpy(xx, launchParams.output_buffer_S,
               68 *
               2*300 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (appContainer.config.projector)
         {
            for (int j = 0; j < 100; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
         printf("out_scalar\n");

         free(xx);*/
				
      
      };
    
      void applyImpl(Tout& out, Tout& out2, const Tin& in, const Tin& in2){
       //printf("%d\n", appContainer.config.Ntheta);
        if (appContainer.config.Ntheta == 0)
        {
            initImpl(out, in);
        }
      Profiler::RegionFixture<5> fix("ParallaltOp::applyImpl");

      //assert(out.size() == in.size());
      //assert(out.dims()[0] == in.dims()[0]);
      assert(out.dims()[1] == in.dims()[1]);
      assert(out.dims()[2] == in.dims()[2]);

      assert(QuICC::Cuda::isDeviceMemory(out.data()));
      assert(QuICC::Cuda::isDeviceMemory(in.data()));

      if (temp_buffer == 0)
      {
         cudaMalloc((void**)&temp_buffer,
            appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
               sizeof(double));
      }
      if (temp_buffer2 == 0)
      {
         cudaMalloc((void**)&temp_buffer2,
            appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
               sizeof(double));
      }
      
        if (temp_buffer3 == 0)
          {
             cudaMalloc((void**)&temp_buffer3,
                appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
                   sizeof(double));
          }
          if (temp_buffer4 == 0)
          {
             cudaMalloc((void**)&temp_buffer4,
                appContainer.config.Ntheta * (appContainer.config.sizeEvenBlock + appContainer.config.sizeOddBlock) *
                   sizeof(double));
          }

            launchParams.input_buffer_S = (double*)in.data();
            launchParams.temp_buffer_S = temp_buffer;
            launchParams.buffer_S = temp_buffer2;
            launchParams.output_buffer_S = (double*)out.data();
            
             launchParams.input_buffer_T = (double*)in2.data();
             launchParams.temp_buffer_T = temp_buffer3;
             launchParams.buffer_T = temp_buffer4;
             launchParams.output_buffer_T = (double*)out2.data();
             launchApp_parallALT(&appContainer, &launchParams);
             /* double* xx =
                (double*)calloc(2 * appContainer.config.Ntheta *
                                         (appContainer.config.sizeEvenBlock +
                                            appContainer.config.sizeOddBlock),
            sizeof(double));
        
       cudaDeviceSynchronize();
          cudaMemcpy(xx, launchParams.input_buffer_S,
            68 *
               2*300 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (!appContainer.config.projector)
         {
            for (int j = 0; j < 128; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
         printf("in_S_store\n");
        
          cudaMemcpy(xx, launchParams.input_buffer_T,
            68 *
               2*300 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (!appContainer.config.projector)
         {
            for (int j = 0; j < 128; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
         printf("in_T\n");
      cudaMemcpy(xx, launchParams.output_buffer_S,
               68 *
               2*300 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (!appContainer.config.projector)
         {
            for (int j = 0; j < 128; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
         printf("out_S\n");
      cudaMemcpy(xx, launchParams.output_buffer_T,
               68 *
               2*300 *
               sizeof(double),
            cudaMemcpyDeviceToHost);
         if (!appContainer.config.projector)
         {
            for (int j = 0; j < 128; j++)
            {
               for (int i = 0; i < 4; i++)
               {
                  printf("%.17e %.17e | ", xx[2 * i+ 2*j * 68], xx[2 * i + 1+ 2*j * 68]);
               }
               printf("\n");
            }
         }
          printf("out_T\n");
         free(xx);*/

     
     // if (testAccuracy) {
				
      
      };
    void setType(int inputType, int inputDirection)
      {
         transformType = inputType;
         transformDirection = inputDirection;
      };
    int getType()
      {
         return transformType;
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

    mutable double* temp_buffer = 0;

    mutable double* temp_buffer2 = 0;

    mutable double* temp_buffer3 = 0;

    mutable double* temp_buffer4 = 0;
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

    mutable int transformType = 0;
    mutable int transformDirection = 0;
    mutable parallALT_launchParams launchParams = {};
    /// @brief Give access to base class
    //friend BinaryBaseOp<DiffOp<Tout, Tin, Order, Direction, Treatment>, Tout, Tin, ScaleType>;

};

} // namespace ALegendre_parallALT
} // namespace Transform
} // namespace QuICC
