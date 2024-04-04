
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

// Passes
#include <Quiccir/IR/QuiccirDialect.h>
#include <Quiccir/Transforms/QuiccirPasses.h>
#include <Quiccir/Pipelines/Passes.h>
#include <Quiccir-c/Utils.h>

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

// QuICC
#include "Graph/MlirShims.hpp"
#include "Graph/BackendsMap.hpp"
#include "Graph/OpsMap.hpp"
#include "Memory/Cpu/NewDelete.hpp"
#include "Memory/Cuda/Malloc.hpp"
#include "Memory/Memory.hpp"


#include "Profiler/Interface.hpp"

int main(int argc, char **argv)
{
  // Register any command line options.
  // mlir::registerAsmPrinterCLOptions();
  // mlir::registerMLIRContextCLOptions();
  // mlir::registerPassManagerCLOptions();

  // llvm::cl::ParseCommandLineOptions(argc, argv, "quiccir jitter\n");

  QuICC::Profiler::Initialize();

  Catch::Session session; // There must be exactly one instance

  auto returnCode = session.run();

  QuICC::Profiler::Finalize();

  return returnCode;
}

TEST_CASE("One Dimensional Loop", "[OneDimLoop]")
{
  mlir::DialectRegistry registry;
  mlir::func::registerAllExtensions(registry);
  // Add the following to include *all* MLIR Core dialects, or selectively
  // include what you need like above. You only need to register dialects that
  // will be *parsed* by the tool, not the one generated
  registerAllDialects(registry);

  mlir::MLIRContext ctx(registry);
  // Load our Dialect in this MLIR Context.
  ctx.loadDialect<mlir::quiccir::QuiccirDialect>();

  std::string modStr =
    "!type_umod = !quiccir.view<5x6x3xf64, \"C_DCCSC3D_t\">\n"
    "!type_uval = !quiccir.view<5x11x3xf64, \"R_DCCSC3D_t\">\n"
    "!type_tumod = tensor<5x6x3xf64, \"C_DCCSC3D_t\">\n"
    "!type_tuval = tensor<5x11x3xf64, \"R_DCCSC3D_t\">\n"
    "func.func @entry(%thisArr: !llvm.ptr<array<2 x ptr>> {llvm.noalias}, %uout: !type_umod, %umod: !type_umod) {\n"
    "  %tumod = builtin.unrealized_conversion_cast %umod : !type_umod to !type_tumod\n"
    "  %tuval = quiccir.fr.prj %tumod : !type_tumod -> !type_tuval attributes{implptr = 0 :i64}\n"
    "  %ret = quiccir.fr.int %tuval : !type_tuval -> !type_tumod attributes{implptr = 1 :i64}\n"
    "  quiccir.materialize %ret in %uout : (!type_tumod, !type_umod)\n"
    "  return\n"
    "}\n";

  // Load
  mlir::OwningOpRef<mlir::ModuleOp> module;
  module = mlir::parseSourceString<mlir::ModuleOp>(modStr, &ctx);

  // Dump
  // module->dump();

  // setup ops map and store
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  QuICC::Graph::MapOps storeOp(*module, mem, /*physDims=*/{11, 3, 5},
    /*modsDims=*/{6, 3, 5});

  // Top level (module) pass manager
  mlir::PassManager pm(&ctx);
  mlir::quiccir::quiccLibCallPipelineBuilder(pm);

  // Lower
  if (mlir::failed(pm.run(*module))) {
    CHECK(false);
  }

  // Dump
  // module->dump();

  // JIT

  // Initialize LLVM targets.
  llvm::InitializeNativeTarget();
  llvm::InitializeNativeTargetAsmPrinter();

  // Register the translation from MLIR to LLVM IR, which must happen before we
  // can JIT-compile.
  mlir::registerBuiltinDialectTranslation(*module->getContext());
  mlir::registerLLVMDialectTranslation(*module->getContext());

  // An optimization pipeline to use within the execution engine.
  auto optPipeline = mlir::makeOptimizingTransformer(
      /*optLevel=*/0, /*sizeLevel=*/0,
      /*targetMachine=*/nullptr);

  // Create an MLIR execution engine. The execution engine eagerly JIT-compiles
  // the module.
  mlir::ExecutionEngineOptions engineOptions;
  engineOptions.transformer = optPipeline;
  // engineOptions.sharedLibPaths = executionEngineLibs;
  auto maybeEngine = mlir::ExecutionEngine::create(*module, engineOptions);
  assert(maybeEngine && "failed to construct an execution engine");
  auto &engine = maybeEngine.get();

  // Invoke the JIT-compiled function.
  std::string symbol = "entry";
  auto funSym = engine->lookup(symbol);
  if (auto E = funSym.takeError()) {
    assert(false && "JIT invocation failed");
  }

  // setup metadata
  constexpr std::size_t modsM = 6;
  constexpr std::size_t N = 3;
  constexpr std::size_t K = 5;
  std::array<std::uint32_t, 3> modsDimensions {modsM, N, K};

  std::array<std::vector<std::uint32_t>, 3> pointers {{{},{0, 2, 2, 2, 3, 4},{}}};
  std::array<std::vector<std::uint32_t>, 3> indices {{{},{0, 1, 0, 0},{}}};

  // host mem block
  std::size_t modsS = modsM*indices[1].size();
  QuICC::Memory::MemBlock<std::complex<double>> modsIn(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  using namespace QuICC::Graph;
  C_DCCSC3D_t modsInView({modsIn.data(), modsIn.size()}, modsDimensions, pointers, indices);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, modsDimensions, pointers, indices);

  // set input modes
  std::complex<double> val = {1.0, 0.0};
  for(std::size_t m = 0; m < modsInView.size(); ++m)
  {
    modsInView[m] = val;
  }

  // map input/output views
  using view3_cd_t = ViewDescriptor<std::complex<double>, std::uint32_t, 3>;
  view3_cd_t viewRef_in{{modsM, N, K},
    pointers[1].data(), (std::uint32_t)pointers[1].size(),
    indices[1].data(), (std::uint32_t)indices[1].size(),
    modsIn.data(), (std::uint32_t)modsIn.size()};
  view3_cd_t viewRef_out{{modsM, N, K},
    pointers[1].data(), (std::uint32_t)pointers[1].size(),
    indices[1].data(), (std::uint32_t)indices[1].size(),
    modsOut.data(), (std::uint32_t)modsOut.size()};

  // get operators map
  auto thisArr = storeOp.getThisArr();

  // Apply graph
  auto fun = (void (*)(void*, view3_cd_t*, view3_cd_t*))funSym.get();
  fun(thisArr.data(), &viewRef_out, &viewRef_in);

  // check
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(modsOutView[m] == val);
  }

}

TEST_CASE("Simple Tree", "[SimpleTree]")
{
  mlir::DialectRegistry registry;
  mlir::func::registerAllExtensions(registry);
  // Add the following to include *all* MLIR Core dialects, or selectively
  // include what you need like above. You only need to register dialects that
  // will be *parsed* by the tool, not the one generated
  registerAllDialects(registry);

  mlir::MLIRContext ctx(registry);
  // Load our Dialect in this MLIR Context.
  ctx.loadDialect<mlir::quiccir::QuiccirDialect>();

  // Load
  mlir::OwningOpRef<mlir::ModuleOp> module;
  std::string inputFilename = "./simple-tree.mlir";
  // Otherwise, the input is '.mlir'.
  llvm::ErrorOr<std::unique_ptr<llvm::MemoryBuffer>> fileOrErr =
      llvm::MemoryBuffer::getFileOrSTDIN(inputFilename);
  if (std::error_code EC = fileOrErr.getError()) {
    llvm::errs() << "Could not open input file: " << EC.message() << "\n";
    CHECK(false);
  }

  // Parse the input mlir.
  llvm::SourceMgr sourceMgr;
  sourceMgr.AddNewSourceBuffer(std::move(*fileOrErr), llvm::SMLoc());
  module = mlir::parseSourceFile<mlir::ModuleOp>(sourceMgr, &ctx);

  // Dump
  // module->dump();

  // Inline and set view attributes
  mlir::PassManager pmPre(&ctx);
  pmPre.addPass(mlir::createInlinerPass());
  std::array<std::array<std::string, 2>, 3> layOpt;
  layOpt[0] = {"R_DCCSC3D_t", "C_DCCSC3D_t"};
  layOpt[1] = {"C_DCCSC3D_t", "C_S1CLCSC3D_t"};
  layOpt[2] = {"C_DCCSC3D_t", "C_DCCSC3D_t"};
  mlir::OpPassManager &nestedFuncPmPre = pmPre.nest<mlir::func::FuncOp>();
  nestedFuncPmPre.addPass(mlir::quiccir::createSetViewLayoutPass(layOpt));

  if (mlir::failed(pmPre.run(*module))) {
    CHECK(false);
  }

  // Dump
  // module->dump();

  // Grid dimensions
  constexpr std::uint32_t rank = 3;
  std::array<std::uint32_t, rank> physDims{10, 6, 3};
  std::array<std::uint32_t, rank> modsDims{6, 3, 2};

  // setup ops map and store
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
  QuICC::Graph::MapOps storeOp(*module, mem, physDims, modsDims);

  // Lower to library call
  mlir::PassManager pm(&ctx);
  mlir::quiccir::quiccLibCallPipelineBuilder(pm);

  // Lower
  if (mlir::failed(pm.run(*module))) {
    CHECK(false);
  }

  // Dump
  // module->dump();

  // JIT

  // Initialize LLVM targets.
  llvm::InitializeNativeTarget();
  llvm::InitializeNativeTargetAsmPrinter();

  // Register the translation from MLIR to LLVM IR, which must happen before we
  // can JIT-compile.
  mlir::registerBuiltinDialectTranslation(*module->getContext());
  mlir::registerLLVMDialectTranslation(*module->getContext());

  // An optimization pipeline to use within the execution engine.
  auto optPipeline = mlir::makeOptimizingTransformer(
      /*optLevel=*/0, /*sizeLevel=*/0,
      /*targetMachine=*/nullptr);

  // Create an MLIR execution engine. The execution engine eagerly JIT-compiles
  // the module.
  mlir::ExecutionEngineOptions engineOptions;
  engineOptions.transformer = optPipeline;
  // engineOptions.sharedLibPaths = executionEngineLibs;
  auto maybeEngine = mlir::ExecutionEngine::create(*module, engineOptions);
  assert(maybeEngine && "failed to construct an execution engine");
  auto &engine = maybeEngine.get();

  // Invoke the JIT-compiled function.
  std::string symbol = "entry";
  auto funSym = engine->lookup(symbol);
  if (auto E = funSym.takeError()) {
    assert(false && "JIT invocation failed");
  }

  // setup metadata
  auto M = physDims[0];
  auto N = physDims[1];
  auto K = physDims[2];
  auto modsM = modsDims[0];
  auto modsN = modsDims[1];
  auto modsK = modsDims[2];
  std::array<std::uint32_t, 3> physDimensions {M, N, K};
  std::array<std::uint32_t, 3> modsDimensions {modsM, modsN, modsK};

  // Populate meta for fully populated tensor
  std::vector<std::uint32_t> ptr(K+1);
  std::vector<std::uint32_t> idx(K*N);
  ptr[0] = 0;
  for (std::size_t i = 1; i < ptr.size(); ++i) {
      ptr[i] = ptr[i-1]+N;
  }
  for (std::size_t i = 0; i < idx.size(); ++i) {
      idx[i] = i % N;
  }

  std::array<std::vector<std::uint32_t>, 3> pointers {{{}, ptr, {}}};
  std::array<std::vector<std::uint32_t>, 3> indices {{{}, idx, {}}};

  // host mem block
  std::size_t modsS = modsM*indices[1].size();
  QuICC::Memory::MemBlock<double> R(modsS, mem.get());
  QuICC::Memory::MemBlock<std::complex<double>> modsOut(modsS, mem.get());

  // host view
  using namespace QuICC::Graph;
  R_DCCSC3D_t RView({R.data(), R.size()}, physDims, pointers, indices);
  C_DCCSC3D_t modsOutView({modsOut.data(), modsOut.size()}, modsDims, pointers, indices);

  // set input modes
  double val = 1.0;
  for(std::size_t m = 0; m < RView.size(); ++m)
  {
    RView[m] = val;
  }

  // map input/output views
  using view3_rd_t = ViewDescriptor<double, std::uint32_t, 3>;
  using view3_cd_t = ViewDescriptor<std::complex<double>, std::uint32_t, 3>;

  view3_rd_t RviewDes{{M, N, K},
    pointers[1].data(), (std::uint32_t)pointers[1].size(),
    indices[1].data(), (std::uint32_t)indices[1].size(),
    R.data(), (std::uint32_t)R.size()};
  view3_rd_t PhiviewDes{{M, N, K},
    pointers[1].data(), (std::uint32_t)pointers[1].size(),
    indices[1].data(), (std::uint32_t)indices[1].size(),
    R.data(), (std::uint32_t)R.size()};
  view3_rd_t ThetaviewDes{{M, N, K},
    pointers[1].data(), (std::uint32_t)pointers[1].size(),
    indices[1].data(), (std::uint32_t)indices[1].size(),
    R.data(), (std::uint32_t)R.size()};
  view3_cd_t viewRef_out{{modsM, N, K},
    pointers[1].data(), (std::uint32_t)pointers[1].size(),
    indices[1].data(), (std::uint32_t)indices[1].size(),
    modsOut.data(), (std::uint32_t)modsOut.size()};

  // get operators map
  auto thisArr = storeOp.getThisArr();

  // Apply graph
  auto fun = (void (*)(void*, view3_cd_t*, view3_rd_t*, view3_rd_t*, view3_rd_t*))funSym.get();
  fun(thisArr.data(), &viewRef_out, &RviewDes, &PhiviewDes, &ThetaviewDes);

  // check
  // for(std::size_t m = 0; m < modsOutView.size(); ++m)
  // {
  //   CHECK(modsOutView[m] == val);
  // }

}