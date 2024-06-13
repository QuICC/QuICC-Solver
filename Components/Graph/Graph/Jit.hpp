/**
 * @file Jit.hpp
 * @brief Mlir JIT wrapper
 */
#pragma once

// External includes
//

// Passes
#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/Transforms/QuiccirPasses.h>
#include <Quiccir/Pipelines/Passes.h>
// #include <Quiccir-c/Utils.h>

#include <mlir/InitAllDialects.h>
#include <mlir/Dialect/Func/Extensions/AllExtensions.h>
#include <mlir/IR/MLIRContext.h>
#include <mlir/Parser/Parser.h>
#include <mlir/Transforms/Passes.h>

#include <mlir/Pass/PassManager.h>
// #include <llvm/Support/CommandLine.h>

// JIT
#include <mlir/ExecutionEngine/ExecutionEngine.h>
#include <mlir/ExecutionEngine/OptUtils.h>
#include <llvm/Support/SourceMgr.h>
#include <llvm/Support/TargetSelect.h>
#include <mlir/Target/LLVMIR/Dialect/Builtin/BuiltinToLLVMIRTranslation.h>
#include <mlir/Target/LLVMIR/Dialect/LLVMIR/LLVMToLLVMIRTranslation.h>

// Project includes
//
#include "Graph/Shims/MlirShims.hpp"

namespace QuICC {
/// @brief This namespace groups the code needed to connect QuICC and quiccir, a MLIR dialect
namespace Graph {

template<class T, class ATT>
ViewDescriptor<T, std::uint32_t, 3> getViewDescriptor(View::View<T, ATT> view)
{
    return ViewDescriptor<T, std::uint32_t, 3>{{view.dims()[0], view.dims()[1], view.dims()[2]},
    view.pointers()[1].data(), (std::uint32_t)view.pointers()[1].size(),
    view.indices()[1].data(), (std::uint32_t)view.indices()[1].size(),
    view.data(), (std::uint32_t)view.size()};
}

/// @brief classe to setup and JIT the mlir graph
/// @tparam RANK spatial dimensions
template <std::uint32_t RANK = 3u>
class Jit {
private:

    /// @brief do we need to add a wrapper
    bool needViewWrapper = false;
    /// @brief storage for MLIR registry
    mlir::DialectRegistry _registry;
    /// @brief storage for MLIR context
    mlir::MLIRContext _ctx;
    /// @brief storage for MLIR module to be JITted
    mlir::OwningOpRef<mlir::ModuleOp> _module;
    /// @brief storage for MLIR execution engine
    std::unique_ptr<mlir::ExecutionEngine> _engine;
    /// @brief  storage for JITted graph function
    void* _funSym;
    /// @brief mat to QuICC operators
    QuICC::Graph::MapOps _storeOp;

    /// @brief register needed MLIR dialects
    void setDialects();
    /// @brief set grid and spectral dimensions
    /// @param physDims
    /// @param modsDims
    /// @param outStage
    /// @param inStage
    void insertWrapper(const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage);
    /// @brief set grid and spectral dimensions
    /// @param physDims
    /// @param modsDims
    void setDimensions(const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims);
    /// @brief set layout for each stage
    /// @param layOpt
    void setLayout(const std::array<std::array<std::string, 2>, 3> layOpt);
    /// @brief map operators
    /// @param physDims
    /// @param modsDims
    void setMap(const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::shared_ptr<Memory::memory_resource> mem);
    /// @brief lower graph to LLVM IR
    void lower();
    /// @brief JIT and store graph function
    void setEngineAndJit();

public:
    /// @brief Setup Graph from string
    /// @param modStr
    /// @param physDims
    /// @param modsDims
    /// @param layOpt
    Jit(const std::string modStr,
        const std::shared_ptr<Memory::memory_resource> mem,
        const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage)
    {
        setDialects();

        // Load module
        _module = mlir::parseSourceString<mlir::ModuleOp>(modStr, &_ctx);
        insertWrapper(physDims, modsDims, layOpt, outStage, inStage);
        setDimensions(physDims, modsDims);
        setLayout(layOpt);
        setMap(physDims, modsDims, mem);
        lower();
        setEngineAndJit();
    }

    /// @brief Setup Graph from source file
    /// @param sourceMgr
    /// @param physDims
    /// @param modsDims
    /// @param layOpt
    Jit(const llvm::SourceMgr sourceMgr,
        const std::shared_ptr<Memory::memory_resource> mem,
        const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage)
    {
        setDialects();

        // Load module
        _module = mlir::parseSourceFile<mlir::ModuleOp>(sourceMgr, &_ctx);
        insertWrapper(physDims, modsDims, layOpt, outStage, inStage);
        setDimensions(physDims, modsDims);
        setLayout(layOpt);
        setMap(physDims, modsDims, mem);
        lower();
        setEngineAndJit();
    }

    template <class Trets, class Targs>
    void apply(Trets ret0V, Targs arg0V)
    {
        // Map view to MLIR view descritor
        auto arg0 = getViewDescriptor(arg0V);
        auto ret0 = getViewDescriptor(ret0V);

        // Get operators map
        auto thisArr = _storeOp.getThisArr();

        // Apply graph
        auto fun = (void (*)(void*, decltype(ret0)*, decltype(arg0)*))_funSym;
        fun(thisArr.data(), &ret0, &arg0);
    }

    template <class Trets, class Targs>
    void apply(Trets ret0V, Targs arg0V, Targs arg1V, Targs arg2V)
    {
        // Map view to MLIR view descritor
        auto arg0 = getViewDescriptor(arg0V);
        auto arg1 = getViewDescriptor(arg1V);
        auto arg2 = getViewDescriptor(arg2V);
        auto ret0 = getViewDescriptor(ret0V);

        // Get operators map
        auto thisArr = _storeOp.getThisArr();

        // Apply graph
        auto fun = (void (*)(void*, decltype(ret0)*,
            decltype(arg0)*,
            decltype(arg1)*,
            decltype(arg2)*
            ))_funSym;
        fun(thisArr.data(), &ret0, &arg0, &arg1, &arg2);
    }
};

template <std::uint32_t RANK>
void Jit<RANK>::setDialects()
{
    mlir::func::registerAllExtensions(_registry);
    // Add the following to include *all* MLIR Core dialects, or selectively
    // include what you need like above. You only need to register dialects that
    // will be *parsed* by the tool, not the one generated
    registerAllDialects(_registry);
    _ctx.appendDialectRegistry(_registry);
    // Load our Dialect in this MLIR Context.
    _ctx.loadDialect<mlir::quiccir::QuiccirDialect>();
}

namespace details
{
template <std::uint32_t RANK>
void setOpt(const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::array<std::array<std::string, 2>, RANK> layout,
    const Stage stage, std::array<std::int64_t, RANK>& dim, std::string& lay)
{

    switch(stage)
    {
        /// MLIR convention has layer first!
        case Stage::PPP :
            dim[1] = physDims[0];
            dim[2] = physDims[1];
            dim[0] = physDims[2];
            lay = layout[0][0];
            break;;
        case Stage::MPP :
            dim[1] = modsDims[0];
            dim[2] = physDims[1];
            dim[0] = physDims[2];
            lay = layout[0][1];
            break;
        case Stage::PPM :
            dim[1] = physDims[1];
            dim[2] = physDims[2];
            dim[0] = modsDims[0];
            lay = layout[1][0];
            break;
        case Stage::MPM :
            dim[1] = modsDims[1];
            dim[2] = physDims[2];
            dim[0] = modsDims[0];
            lay = layout[1][1];
            break;
        case Stage::PMM :
            dim[1] = physDims[2];
            dim[2] = modsDims[0];
            dim[0] = modsDims[1];
            lay = layout[2][0];
            break;
        case Stage::MMM :
            dim[1] = modsDims[2];
            dim[2] = modsDims[0];
            dim[0] = modsDims[1];
            lay = layout[2][1];
            break;
    }
}
}

template <std::uint32_t RANK>
void Jit<RANK>::insertWrapper(const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::array<std::array<std::string, 2>, RANK> lay,
    const Stage outStage, const Stage inStage)
{
    // Inline and insert wrapper
    mlir::PassManager pmPre(&_ctx);
    pmPre.addPass(mlir::createInlinerPass());
    mlir::quiccir::QuiccirViewWrapperOptions opt;

    std::array<std::int64_t, RANK> dimArgs;
    details::setOpt<RANK>(physDims, modsDims, lay, inStage, dimArgs, opt.layArgs);
    opt.dimArgs = dimArgs;
    std::array<std::int64_t, RANK> dimRets;
    details::setOpt<RANK>(physDims, modsDims, lay, outStage, dimRets, opt.layRets);
    opt.dimRets = dimRets;

    pmPre.addPass(mlir::quiccir::createViewWrapperPass(opt));

    if (mlir::failed(pmPre.run(*_module))) {
        throw std::runtime_error("Failed to insert wrapper.");
    }
}


template <std::uint32_t RANK>
void Jit<RANK>::setDimensions(const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims)
{
    // Inline and set dimensions
    mlir::PassManager pmPre(&_ctx);
    pmPre.addPass(mlir::createInlinerPass());
    mlir::OpPassManager &nestedFuncPmPre = pmPre.nest<mlir::func::FuncOp>();
    // Reorder dimensions based on mlir convention
    std::array<std::int64_t, RANK> phys{physDims[2], physDims[0], physDims[1]};
    std::array<std::int64_t, RANK> mods{modsDims[2], modsDims[0], modsDims[1]};
    nestedFuncPmPre.addPass(mlir::quiccir::createSetDimensionsPass(phys, mods));

    if (mlir::failed(pmPre.run(*_module))) {
        throw std::runtime_error("Failed to set dimensions.");
    }
}

template <std::uint32_t RANK>
void Jit<RANK>::setLayout(const std::array<std::array<std::string, 2>, 3> layOpt)
{
    // Inline and set layout attributes
    mlir::PassManager pmPre(&_ctx);
    mlir::OpPassManager &nestedFuncPmPre = pmPre.nest<mlir::func::FuncOp>();
    nestedFuncPmPre.addPass(mlir::quiccir::createSetViewLayoutPass(layOpt));
    pmPre.addPass(mlir::createCanonicalizerPass());

    if (mlir::failed(pmPre.run(*_module))) {
        throw std::runtime_error("Failed to set layout.");
    }
}

template <std::uint32_t RANK>
void Jit<RANK>::setMap(const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::shared_ptr<Memory::memory_resource> mem)
{
    // setup ops map and store
    _storeOp = std::move(QuICC::Graph::MapOps(*_module, mem, physDims, modsDims));
}

template <std::uint32_t RANK>
void Jit<RANK>::lower()
{
    // Top level (module) pass manager
    mlir::PassManager pm(&_ctx);
    mlir::quiccir::quiccLibCallPipelineBuilder(pm);

    // Lower
    if (mlir::failed(pm.run(*_module))) {
        throw std::runtime_error("Failed to lower module.");
    }
}

template <std::uint32_t RANK>
void Jit<RANK>::setEngineAndJit()
{
    _module->dump();

    // Initialize LLVM targets.
    llvm::InitializeNativeTarget();
    llvm::InitializeNativeTargetAsmPrinter();

    // Register the translation from MLIR to LLVM IR, which must happen before we
    // can JIT-compile.
    mlir::registerBuiltinDialectTranslation(*_module->getContext());
    mlir::registerLLVMDialectTranslation(*_module->getContext());

    // An optimization pipeline to use within the execution engine.
    auto optPipeline = mlir::makeOptimizingTransformer(
        /*optLevel=*/0, /*sizeLevel=*/0,
        /*targetMachine=*/nullptr);

    // Create an MLIR execution engine. The execution engine eagerly JIT-compiles
    // the module.
    mlir::ExecutionEngineOptions engineOptions;
    engineOptions.transformer = optPipeline;
    // engineOptions.sharedLibPaths = executionEngineLibs;
    auto maybeEngine = mlir::ExecutionEngine::create(*_module, engineOptions);
    if (!maybeEngine) {
        throw std::runtime_error("failed to construct an execution engine");
    }
    auto& engine = maybeEngine.get();
    _engine = std::move(engine);

    // Invoke the JIT-compiled function.
    std::string symbol = "_view_entry";
    auto funSym = _engine->lookup(symbol);
    if (auto E = funSym.takeError()) {
        throw std::runtime_error("JIT invocation failed");
    }
    _funSym = funSym.get();
}


} // namespace Graph
} // namespace QuICC
