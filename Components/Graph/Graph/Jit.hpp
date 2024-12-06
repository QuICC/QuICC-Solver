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
#include <mlir/Debug/CLOptionsSetup.h>

// Project includes
//
#include "Graph/OpsMap.hpp"
#include "Graph/Tags.hpp"
#include "Graph/Types.hpp"
#include "Graph/Shims/MlirShims.hpp"

namespace QuICC {
/// @brief This namespace groups the code needed to connect QuICC and quiccir, a MLIR dialect
namespace Graph {

/// @brief Return a ViewDescritort struct
/// (which is a 1:1 map to view type in quiccir)
/// from a QuICC::View::View
/// @tparam T scalar type
/// @tparam ATT attribute encoding (not used at the moment)
/// @param view
/// @return
template<class T, class ATT>
ViewDescriptor<T, std::uint32_t, 3> getViewDescriptor(View::View<T, ATT> view)
{
    return ViewDescriptor<T, std::uint32_t, 3>{{view.dims()[0], view.dims()[1], view.dims()[2]},
    view.pointers()[1].data(), (std::uint32_t)view.pointers()[1].size(),
    view.indices()[1].data(), (std::uint32_t)view.indices()[1].size(),
    view.data(), (std::uint32_t)view.size()};
}

/// @brief collection of pass options used in the
/// lowering pipeline
struct PipelineOptions
{
    /// @brief Wrapper pass options
    mlir::quiccir::QuiccirViewWrapperOptions wrap;
};

/// @brief classe to setup and JIT the mlir graph
/// @tparam RANK spatial dimensions
template <std::uint32_t RANK = 3u>
class Jit {
private:

    /// @brief storage for MLIR registry
    mlir::DialectRegistry _registry;
    /// @brief storage for MLIR context
    mlir::MLIRContext _ctx;

    // mlir::tracing::DebugConfig _debugConfig;
    /// @brief storage for MLIR module to be JITted
    mlir::OwningOpRef<mlir::ModuleOp> _module;
    /// @brief storage for MLIR execution engine
    std::unique_ptr<mlir::ExecutionEngine> _engine;
    /// @brief  storage for JITted graph function
    void* _funSym;
    /// @brief mat to QuICC operators
    QuICC::Graph::MapOps _storeOp;

    /// @brief storage for intermediate stages
    /// metadata idx/ptr pairs
    /// (0, 1) -> stage 0
    /// (2, 3) -> stage 1
    /// (4, 5) -> stage 2
    std::array<MemRefDescriptor<std::uint32_t, 1>*, 6> _metaMap;
    std::array<MemRefDescriptor<std::uint32_t, 1>, 6> _metaStore;

    /// @brief storage for scaling parameters
    PhysicalParameters<double> _physParams;

    /// @brief storage for pipeline pass options
    PipelineOptions _opt;

    /// @brief register needed MLIR dialects
    void setDialects();
    /// @brief set grid and spectral dimensions
    /// @param physDims physical dimensions
    /// order v012, RThetaPhi
    /// @param modsDims spectral sapce dimensions
    /// order v012, NLM
    /// @param outStage
    /// @param inStage
    void insertWrapper(const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage);

    /// @brief set grid and spectral dimensions
    /// @param physDims physical dimensions
    /// order v012, RThetaPhi
    /// @param modsDims spectral sapce dimensions
    /// order v012, NLM
    void setDimensions(const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims);

    /// @brief set layout for each stage
    /// @param layOpt
    void setLayout(const std::array<std::array<std::string, 2>, 3> layOpt);

    /// @brief set meta data for intermediate stages
    void setMeta(const std::vector<View::ViewBase<std::uint32_t>>& meta);

    /// @brief map operators
    /// @param mem memory resource to be passed to operators
    void setMap(const std::shared_ptr<Memory::memory_resource> mem);

    /// @brief lower graph to LLVM IR
    void lower();

    /// @brief JIT and store graph function
    void setEngineAndJit();

    /// @brief Passes to set the wrappers and dimensions
    /// @param mem memory resource
    /// @param physDims physical dimensions
    /// order v012, RThetaPhi
    /// @param modsDims spectral sapce dimensions
    /// order v012, NLM
    /// @param layOpt
    /// @param outStage
    /// @param inStage
    /// @param meta
    void setWrappers(const std::shared_ptr<Memory::memory_resource> mem,
        const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage,
        const std::vector<View::ViewBase<std::uint32_t>>& meta);

public:
    ~Jit()
    {
        /// \todo clearup meta
    }

    /// @brief Setup Graph from string
    /// @param modStr string describing mlir module
    /// @param mem memory resource
    /// @param physDims physical dimensions
    /// order v012, RThetaPhi
    /// @param modsDims spectral sapce dimensions
    /// order v012, NLM
    /// @param layOpt
    /// @param outStage
    /// @param inStage
    /// @param meta
    Jit(const std::string modStr,
        const std::shared_ptr<Memory::memory_resource> mem,
        const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage,
        const std::vector<View::ViewBase<std::uint32_t>>& meta,
        const PhysicalParameters<double> physParams = PhysicalParameters(),
        const PipelineOptions options = PipelineOptions());

    /// @brief Setup Graph from source file
    /// @param sourceMgr source manager containing mlir module
    /// @param mem memory resource
    /// @param physDims physical dimensions
    /// order v012, RThetaPhi
    /// @param modsDims spectral sapce dimensions
    /// order v012, NLM
    /// @param layOpt
    /// @param outStage
    /// @param inStage
    /// @param meta
    Jit(const llvm::SourceMgr sourceMgr,
        const std::shared_ptr<Memory::memory_resource> mem,
        const std::array<std::uint32_t, RANK> physDims,
        const std::array<std::uint32_t, RANK> modsDims,
        const std::array<std::array<std::string, 2>, RANK> layOpt,
        const Stage outStage, const Stage inStage,
        const std::vector<View::ViewBase<std::uint32_t>>& meta,
        const PhysicalParameters<double> physParams = PhysicalParameters(),
        const PipelineOptions options = PipelineOptions());

    /// @brief Apply graph
    /// Catch all, supposed to fail
    /// @tparam Targs...
    /// @param args
    template <class... Targs>
    void apply(Targs...);

    /// @brief Apply graph
    /// 1 input, 1 output
    /// @tparam Trets
    /// @tparam Targs
    /// @param ret0V
    /// @param arg0V
    template <class Trets, class Targs>
    void apply(Trets ret0V, Targs arg0V);

    /// @brief Apply graph
    /// 3 inputs, 1 output
    /// @tparam Trets
    /// @tparam Targs
    /// @param ret0V
    /// @param arg0V
    /// @param arg1V
    /// @param arg2V
    template <class Trets, class Targs>
    void apply(Trets ret0V, Targs arg0V, Targs arg1V, Targs arg2V);

    /// @brief Apply graph
    /// 3 inputs, 3 outputs
    /// @tparam Trets
    /// @tparam Targs
    /// @param ret0V
    /// @param ret1V
    /// @param ret2V
    /// @param arg0V
    /// @param arg1V
    /// @param arg2V
    template <class Trets, class Targs>
    void apply(Trets ret0V, Trets ret1V, Trets ret2V, Targs arg0V, Targs arg1V, Targs arg2V);

    /// @brief Apply graph
    /// 3 inputs, 3+3 outputs
    /// Pure hydro models returning also physical space velocity
    /// @tparam TyOnerets
    /// @tparam TyTworets
    /// @tparam Targs
    /// @param retTyOne0V
    /// @param retTyOne1V
    /// @param retTyOne2V
    /// @param arg0V
    /// @param arg1V
    /// @param arg2V
    template <class TyOnerets, class TyTworets, class Targs>
    void apply(
        TyOnerets retTyOne0V, TyOnerets retTyOne1V, TyOnerets retTyOne2V,
        TyTworets retTyTwo0V, TyTworets retTyTwo1V, TyTworets retTyTwo2V,
        Targs arg0V, Targs arg1V, Targs arg2V);

    /// @brief Apply graph
    /// 5 inputs, 5+6 outputs
    /// Pure magneto-hydro models returning also physical space velocity
    /// @tparam TyOnerets
    /// @tparam TyTworets
    /// @tparam Targs
    /// @param retTyOne0V
    /// @param retTyOne1V
    /// @param retTyOne2V
    /// @param retTyOne3V
    /// @param retTyOne4V
    /// @param retTyTwo0V
    /// @param retTyTwo1V
    /// @param retTyTwo2V
    /// @param retTyTwo3V
    /// @param retTyTwo4V
    /// @param retTyTwo5V
    /// @param arg0V
    /// @param arg1V
    /// @param arg2V
    /// @param arg3V
    /// @param arg4V
    template <class TyOnerets, class TyTworets, class Targs>
    void apply(
        TyOnerets retTyOne0V, TyOnerets retTyOne1V, TyOnerets retTyOne2V,
        TyOnerets retTyOne3V, TyOnerets retTyOne4V,
        TyTworets retTyTwo0V, TyTworets retTyTwo1V, TyTworets retTyTwo2V,
        TyTworets retTyTwo3V, TyTworets retTyTwo4V, TyTworets retTyTwo5V,
        Targs arg0V, Targs arg1V, Targs arg2V, Targs arg3V, Targs arg4V);
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

    // #ifndef NDEBUG
    // mlir::tracing::InstallDebugHandler installDebugHandler(_ctx,
    //                                                _debugConfig);
    // #endif
}

namespace details
{

/// @brief Set options for  ViewWrapperPass in mlir
/// @tparam RANK
/// @param physDims physical dimensions
/// order v012, RThetaPhi
/// @param modsDims spectral sapce dimensions
/// order v012, NLM
/// @param layout
/// @param stage
/// @param dim mlir uses v210 with layer in first position (v021)
/// @param lay
template <std::uint32_t RANK>
void setOpt(const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::array<std::array<std::string, 2>, RANK> layout,
    const Stage stage, std::vector<std::int64_t>& dim, std::string& lay)
{
    dim.resize(RANK);
    switch(stage)
    {
        /// MLIR convention has layer first!
        case Stage::PPP :
            // v210
            dim[1] = physDims[2];
            dim[2] = physDims[1];
            dim[0] = physDims[0];
            lay = layout[0][0];
            break;;
        case Stage::MPP :
            // v210
            dim[1] = modsDims[2];
            dim[2] = physDims[1];
            dim[0] = physDims[0];
            lay = layout[0][1];
            break;
        case Stage::PPM :
            // v102
            dim[1] = physDims[1];
            dim[2] = physDims[0];
            dim[0] = modsDims[2];
            lay = layout[1][0];
            break;
        case Stage::MPM :
            // v102
            dim[1] = modsDims[1];
            dim[2] = physDims[0];
            dim[0] = modsDims[2];
            lay = layout[1][1];
            break;
        case Stage::PMM :
            // v021
            dim[1] = physDims[0];
            dim[2] = modsDims[2];
            dim[0] = modsDims[1];
            lay = layout[2][0];
            break;
        case Stage::MMM :
            // v021
            dim[1] = modsDims[0];
            dim[2] = modsDims[2];
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
    mlir::OpPassManager &nestedFuncPmPre = pmPre.nest<mlir::func::FuncOp>();
    nestedFuncPmPre.addPass(mlir::quiccir::createTransformContractionPass());

    std::vector<std::vector<std::int64_t>> dimArgs(1);
    std::vector<std::string> layArgs(1);
    details::setOpt<RANK>(physDims, modsDims, lay, inStage, dimArgs[0], layArgs[0]);
    _opt.wrap.dimArgs = dimArgs;
    _opt.wrap.layArgs = layArgs;
    std::vector<std::vector<std::int64_t>> dimRets(1);
    std::vector<std::string> layRets(1);
    details::setOpt<RANK>(physDims, modsDims, lay, outStage, dimRets[0], layRets[0]);
    // Set to defeault if not set
    if (_opt.wrap.dimRets.size() == 0)
    {
        _opt.wrap.dimRets = dimRets;
    }
    if (_opt.wrap.layRets.size() == 0)
    {
        _opt.wrap.layRets = layRets;
    }

    pmPre.addPass(mlir::quiccir::createViewWrapperPass(_opt.wrap));

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
    // mlir uses v210 with layer in first position (v021)
    std::array<std::int64_t, RANK> phys{physDims[0], physDims[2], physDims[1]};
    std::array<std::int64_t, RANK> mods{modsDims[0], modsDims[2], modsDims[1]};
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
void Jit<RANK>::setMeta(const std::vector<View::ViewBase<std::uint32_t>>& meta)
{
    for (std::size_t i = 0; i < _metaMap.size(); ++i)
    {
        _metaMap[i] = nullptr;
    }
    if (meta.size() > 0)
    {
        //  std::array<MemRefDescriptor<std::uint32_t, 1>*, 6> _metaMap;
        // ptr, ptr, offset, {size}, {stride}
        // idx Fr
        _metaStore[0] = {meta[0].data(), meta[0].data(), 0, {static_cast<intptr_t>(meta[0].size())}, {1}};
        _metaMap[0] = &_metaStore[0];
        // ptr Fr
        _metaStore[1] = {meta[1].data(), meta[1].data(), 0, {static_cast<intptr_t>(meta[1].size())}, {1}};
        _metaMap[1] = &_metaStore[1];
        // idx AL
        _metaStore[2] = {meta[2].data(), meta[2].data(), 0, {static_cast<intptr_t>(meta[2].size())}, {1}};
        _metaMap[2] = &_metaStore[2];
        // ptr AL
        _metaStore[3] = {meta[3].data(), meta[3].data(), 0, {static_cast<intptr_t>(meta[3].size())}, {1}};
        _metaMap[3] = &_metaStore[3];
        // idx JW
        _metaStore[4] = {meta[4].data(), meta[4].data(), 0, {static_cast<intptr_t>(meta[4].size())}, {1}};
        _metaMap[4] = &_metaStore[4];
        // ptr JW
        _metaStore[5] = {meta[5].data(), meta[5].data(), 0, {static_cast<intptr_t>(meta[5].size())}, {1}};
        _metaMap[5] = &_metaStore[5];
    }
}

template <std::uint32_t RANK>
void Jit<RANK>::setMap(const std::shared_ptr<Memory::memory_resource> mem)
{
    // Top level (module) pass manager
    mlir::PassManager pm(&_ctx);
    // Add implptr to mlir
    pm.addPass(mlir::quiccir::createSetImplptrPass());
    if (mlir::failed(pm.run(*_module))) {
        throw std::runtime_error("Failed to add implementation pointer.");
    }
    // _module->dump();

    // setup ops map and store
    _storeOp = std::move(QuICC::Graph::MapOps(*_module, _physParams, mem));
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
    // _module->dump();

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

template <std::uint32_t RANK>
Jit<RANK>::Jit(const std::string modStr,
    const std::shared_ptr<Memory::memory_resource> mem,
    const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::array<std::array<std::string, 2>, RANK> layOpt,
    const Stage outStage, const Stage inStage,
    const std::vector<View::ViewBase<std::uint32_t>>& meta,
    const PhysicalParameters<double> physParams,
    const PipelineOptions options) : _physParams(physParams), _opt(options)
{
    setDialects();
    // Load module
    _module = mlir::parseSourceString<mlir::ModuleOp>(modStr, &_ctx);
    // Passes
    setWrappers(mem, physDims, modsDims, layOpt, outStage, inStage, meta);
    lower();
    // Jit
    setEngineAndJit();
}

template <std::uint32_t RANK>
Jit<RANK>::Jit(const llvm::SourceMgr sourceMgr,
    const std::shared_ptr<Memory::memory_resource> mem,
    const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::array<std::array<std::string, 2>, RANK> layOpt,
    const Stage outStage, const Stage inStage,
    const std::vector<View::ViewBase<std::uint32_t>>& meta,
    const PhysicalParameters<double> physParams,
    const PipelineOptions options) : _physParams(physParams), _opt(options)
{
    setDialects();
    // Load module
    _module = mlir::parseSourceFile<mlir::ModuleOp>(sourceMgr, &_ctx);
    // Passes
    setWrappers(mem, physDims, modsDims, layOpt, outStage, inStage, meta);
    lower();
    // Jit;
    setEngineAndJit();
}

template <std::uint32_t RANK>
void Jit<RANK>::setWrappers(const std::shared_ptr<Memory::memory_resource> mem,
    const std::array<std::uint32_t, RANK> physDims,
    const std::array<std::uint32_t, RANK> modsDims,
    const std::array<std::array<std::string, 2>, RANK> layOpt,
    const Stage outStage, const Stage inStage,
    const std::vector<View::ViewBase<std::uint32_t>>& meta)
{
    insertWrapper(physDims, modsDims, layOpt, outStage, inStage);
    setDimensions(physDims, modsDims);
    setLayout(layOpt);
    setMap(mem);
    setMeta(meta);
}

template <std::uint32_t RANK>
template <class... Targs>
void Jit<RANK>::apply(Targs...)
{
    throw std::logic_error("Graph::apply not implemented");
}

template <std::uint32_t RANK>
template <class Trets, class Targs>
void Jit<RANK>::apply(Trets ret0V, Targs arg0V)
{
    // Map view to MLIR view descritor
    auto arg0 = getViewDescriptor(arg0V);
    auto ret0 = getViewDescriptor(ret0V);

    // Get operators map
    auto thisArr = _storeOp.getThisArr();

    // Apply graph
    auto fun = (void (*)(void*, void*, decltype(ret0)*, decltype(arg0)*))_funSym;
    fun(_metaMap.data(), thisArr.data(), &ret0, &arg0);
}

template <std::uint32_t RANK>
template <class Trets, class Targs>
void Jit<RANK>::apply(Trets ret0V, Targs arg0V, Targs arg1V, Targs arg2V)
{
    // Map view to MLIR view descritor
    auto arg0 = getViewDescriptor(arg0V);
    auto arg1 = getViewDescriptor(arg1V);
    auto arg2 = getViewDescriptor(arg2V);
    auto ret0 = getViewDescriptor(ret0V);

    // Get operators map
    auto thisArr = _storeOp.getThisArr();

    // Apply graph
    auto fun = (void (*)(void*, void*, decltype(ret0)*,
        decltype(arg0)*,
        decltype(arg1)*,
        decltype(arg2)*
        ))_funSym;
    fun(_metaMap.data(), thisArr.data(), &ret0, &arg0, &arg1, &arg2);
}

template <std::uint32_t RANK>
template <class Trets, class Targs>
void Jit<RANK>::apply(Trets ret0V, Trets ret1V, Trets ret2V, Targs arg0V, Targs arg1V, Targs arg2V)
{
    // Map view to MLIR view descritor
    auto arg0 = getViewDescriptor(arg0V);
    auto arg1 = getViewDescriptor(arg1V);
    auto arg2 = getViewDescriptor(arg2V);
    auto ret0 = getViewDescriptor(ret0V);
    auto ret1 = getViewDescriptor(ret1V);
    auto ret2 = getViewDescriptor(ret2V);

    // Get operators map
    auto thisArr = _storeOp.getThisArr();

    // Apply graph
    auto fun = (void (*)(void*, void*,
        decltype(ret0)*,
        decltype(ret1)*,
        decltype(ret2)*,
        decltype(arg0)*,
        decltype(arg1)*,
        decltype(arg2)*
        ))_funSym;
    fun(_metaMap.data(), thisArr.data(),
        &ret0, &ret1, &ret2,
        &arg0, &arg1, &arg2);

}

template <std::uint32_t RANK>
template <class TyOnerets, class TyTworets, class Targs>
void Jit<RANK>::apply(
    TyOnerets retTyOne0V, TyOnerets retTyOne1V, TyOnerets retTyOne2V,
    TyTworets retTyTwo0V, TyTworets retTyTwo1V, TyTworets retTyTwo2V,
    Targs arg0V, Targs arg1V, Targs arg2V)
{
    // Map view to MLIR view descritor
    auto arg0 = getViewDescriptor(arg0V);
    auto arg1 = getViewDescriptor(arg1V);
    auto arg2 = getViewDescriptor(arg2V);
    auto ret0 = getViewDescriptor(retTyOne0V);
    auto ret1 = getViewDescriptor(retTyOne1V);
    auto ret2 = getViewDescriptor(retTyOne2V);
    auto ret3 = getViewDescriptor(retTyTwo0V);
    auto ret4 = getViewDescriptor(retTyTwo1V);
    auto ret5 = getViewDescriptor(retTyTwo2V);

    // Get operators map
    auto thisArr = _storeOp.getThisArr();

    // Apply graph
    auto fun = (void (*)(void*, void*,
        decltype(ret0)*,
        decltype(ret1)*,
        decltype(ret2)*,
        decltype(ret3)*,
        decltype(ret4)*,
        decltype(ret5)*,
        decltype(arg0)*,
        decltype(arg1)*,
        decltype(arg2)*
        ))_funSym;
    fun(_metaMap.data(), thisArr.data(),
        &ret0, &ret1, &ret2,
        &ret3, &ret4, &ret5,
        &arg0, &arg1, &arg2);

}

template <std::uint32_t RANK>
template <class TyOnerets, class TyTworets, class Targs>
void Jit<RANK>::apply(
    TyOnerets retTyOne0V, TyOnerets retTyOne1V, TyOnerets retTyOne2V,
    TyOnerets retTyOne3V, TyOnerets retTyOne4V,
    TyTworets retTyTwo0V, TyTworets retTyTwo1V, TyTworets retTyTwo2V,
    TyTworets retTyTwo3V, TyTworets retTyTwo4V, TyTworets retTyTwo5V,
    Targs arg0V, Targs arg1V, Targs arg2V,
    Targs arg3V, Targs arg4V)
{
    // Map view to MLIR view descritor
    auto arg0 = getViewDescriptor(arg0V);
    auto arg1 = getViewDescriptor(arg1V);
    auto arg2 = getViewDescriptor(arg2V);
    auto arg3 = getViewDescriptor(arg3V);
    auto arg4 = getViewDescriptor(arg4V);
    auto ret0 = getViewDescriptor(retTyOne0V);
    auto ret1 = getViewDescriptor(retTyOne1V);
    auto ret2 = getViewDescriptor(retTyOne2V);
    auto ret3 = getViewDescriptor(retTyOne3V);
    auto ret4 = getViewDescriptor(retTyOne4V);
    auto ret5 = getViewDescriptor(retTyTwo0V);
    auto ret6 = getViewDescriptor(retTyTwo1V);
    auto ret7 = getViewDescriptor(retTyTwo2V);
    auto ret8 = getViewDescriptor(retTyTwo3V);
    auto ret9 = getViewDescriptor(retTyTwo4V);
    auto ret10 = getViewDescriptor(retTyTwo5V);

    // Get operators map
    auto thisArr = _storeOp.getThisArr();

    // Apply graph
    auto fun = (void (*)(void*, void*,
        decltype(ret0)*,
        decltype(ret1)*,
        decltype(ret2)*,
        decltype(ret3)*,
        decltype(ret4)*,
        decltype(ret5)*,
        decltype(ret6)*,
        decltype(ret7)*,
        decltype(ret8)*,
        decltype(ret9)*,
        decltype(ret10)*,
        decltype(arg0)*,
        decltype(arg1)*,
        decltype(arg2)*,
        decltype(arg3)*,
        decltype(arg4)*
        ))_funSym;
    fun(_metaMap.data(), thisArr.data(),
        &ret0, &ret1, &ret2,
        &ret3, &ret4,
        &ret5, &ret6, &ret7,
        &ret8, &ret9, &ret10,
        &arg0, &arg1, &arg2,
        &arg3, &arg4);
}


} // namespace Graph
} // namespace QuICC
