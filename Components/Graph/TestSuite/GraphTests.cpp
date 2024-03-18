
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
  auto mem = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
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

  // map views
  using view3_cd_t = ViewDescriptor<std::complex<double>, std::uint32_t, 3>;
  view3_cd_t viewRef_in{{modsM, N, K}, pointers[1].data(), pointers[1].size(), indices[1].data(),
    indices[1].size(), modsIn.data(), modsIn.size()};
  view3_cd_t viewRef_out{{modsM, N, K}, pointers[1].data(), pointers[1].size(), indices[1].data(),
    indices[1].size(), modsOut.data(), modsOut.size()};

  // get operators and map
  std::array<void*, 2> thisArr;
  using namespace QuICC::Transform::Fourier;
  using backend_t = QuICC::Graph::viewCpu_t;
  using Tin = C_DCCSC3D_t;
  using Tout = R_DCCSC3D_t;
  using backendFft_t = Fft_t<backend_t, Tout, Tin>;
  using backendDiff_t = MixedDiff_t<backend_t, Tin, 0, bwd_t,
    QuICC::Transform::Fourier::none_m>;
  using op_t = Mixed::Projector::DOp<Tout, Tin, backendFft_t,
  backendDiff_t>;

  op_t prjOp(mem);
  thisArr[0] = (void*)&prjOp;

  // Apply graph
  auto fun = (void (*)(void*, view3_cd_t*, view3_cd_t*))funSym.get();
  fun(thisArr.data(), &viewRef_out, &viewRef_in);

  // check
  for(std::size_t m = 0; m < modsOutView.size(); ++m)
  {
    CHECK(modsOutView[m] == val);
  }

}