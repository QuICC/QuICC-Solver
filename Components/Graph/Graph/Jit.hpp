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
namespace Graph {

template<class T, class ATT>
ViewDescriptor<T, std::uint32_t, 3> getViewDescriptor(View::View<T, ATT> view)
{
    return ViewDescriptor<T, std::uint32_t, 3>{{view.dims()[0], view.dims()[1], view.dims()[2]},
    view.pointers()[1].data(), (std::uint32_t)view.pointers()[1].size(),
    view.indices()[1].data(), (std::uint32_t)view.indices()[1].size(),
    view.data(), (std::uint32_t)view.size()};
}

template <std::uint32_t RANK = 3u>
class Jit {
private:
    /// @brief data
    mlir::DialectRegistry _registry;
    mlir::MLIRContext _ctx;
    mlir::OwningOpRef<mlir::ModuleOp> _module;
    void* _funSym;
    std::unique_ptr<mlir::ExecutionEngine> _engine;

    QuICC::Graph::MapOps _storeOp;
    /// @brief memory resource for internal op allocation
    std::shared_ptr<Memory::memory_resource> _mem;


    void setDialects()
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

    void setLayout(const std::array<std::array<std::string, 2>, 3> layOpt)
    {
        // Inline and set view attributes
        mlir::PassManager pmPre(&_ctx);
        pmPre.addPass(mlir::createInlinerPass());
        mlir::OpPassManager &nestedFuncPmPre = pmPre.nest<mlir::func::FuncOp>();
        nestedFuncPmPre.addPass(mlir::quiccir::createSetViewLayoutPass(layOpt));

        if (mlir::failed(pmPre.run(*_module))) {
            throw std::runtime_error("Failed to set layout.");
        }
    }

    void setMap(const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims)
    {
        // setup ops map and store
        _storeOp = std::move(QuICC::Graph::MapOps(*_module, _mem, physDims, modsDims));
    }

    void lower()
    {
        // Top level (module) pass manager
        mlir::PassManager pm(&_ctx);
        mlir::quiccir::quiccLibCallPipelineBuilder(pm);

        // Lower
        if (mlir::failed(pm.run(*_module))) {
            throw std::runtime_error("Failed to lower module.");
        }
    }

    void setEngineAndJit()
    {
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
        std::string symbol = "entry";
        auto funSym = _engine->lookup(symbol);
        if (auto E = funSym.takeError()) {
           throw std::runtime_error("JIT invocation failed");
        }
        _funSym = funSym.get();
    }

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
        const std::array<std::array<std::string, 2>, RANK> layOpt = {{}}) : _mem(mem)
    {
        setDialects();

        // Load module
        _module = mlir::parseSourceString<mlir::ModuleOp>(modStr, &_ctx);

        if (layOpt[0][0].size() > 0)
        {
            setLayout(layOpt);
        }

        setMap(physDims, modsDims);
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
        const std::array<std::array<std::string, 2>, RANK> layOpt = {{}}) : _mem(mem)
    {
        setDialects();

        // Load module
        _module = mlir::parseSourceFile<mlir::ModuleOp>(sourceMgr, &_ctx);

        if (layOpt[0][0].size() > 0)
        {
            setLayout(layOpt);
        }

        setMap(physDims, modsDims);
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


} // namespace Graph
} // namespace QuICC
