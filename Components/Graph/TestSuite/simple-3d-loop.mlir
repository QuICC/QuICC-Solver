// mod -> phys
// phys -> mod
// Nr x Nphi x Ntheta
// !type_uval = !quiccir.view<6x10x10xf64, "DCCSC3D">
// L x N x M  ( ..., radial, ...)
// !type_umod = !quiccir.view<7x3x6xcomplex<f64>, "DCCSC3D">


func.func private @bwd(%Pol: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xf64>) {
  // Pol
  %Pol1 = quiccir.jw.prj %Pol : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{kind = "P"}
  %Pol1T = quiccir.transpose %Pol1 permutation = [1, 2, 0] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>>
  %Pol2 = quiccir.al.prj %Pol1T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{kind = "P"}
  %Pol2T = quiccir.transpose %Pol2 permutation = [1, 2, 0] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>>
  %Pol3 = quiccir.fr.prj %Pol2T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xf64> attributes{kind = "P"}
  return %Pol3 : tensor<?x?x?xf64>
}

func.func private @fwd(%R: tensor<?x?x?xf64>) -> (tensor<?x?x?xcomplex<f64>>) {
  // R
  %R1 = quiccir.fr.int %R : tensor<?x?x?xf64> -> tensor<?x?x?xcomplex<f64>> attributes{kind = "P"}
  %R1T = quiccir.transpose %R1 permutation = [2, 0, 1] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>>
  %R2 = quiccir.al.int %R1T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{kind = "P"}
  %R2T = quiccir.transpose %R2 permutation = [2, 0, 1] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>>
  %R3 = quiccir.jw.int %R2T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{kind = "P"}
  return %R3 : tensor<?x?x?xcomplex<f64>>
}

func.func @entry(%Polur: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xcomplex<f64>>) {
  %tmp = call @bwd(%Polur) : (tensor<?x?x?xcomplex<f64>>) -> tensor<?x?x?xf64>
  %PolNewur = call @fwd(%tmp) : (tensor<?x?x?xf64>) -> tensor<?x?x?xcomplex<f64>>
  return %PolNewur : tensor<?x?x?xcomplex<f64>>
}

// ./bin/quiccir-opt ../examples/simple-3d-loop.mlir --inline --quiccir-view-wrapper='dim-rets=7,3,6 dim-args=7,3,6 lay-args=DCCSC3D lay-rets=DCCSC3D' --inline --set-quiccir-dims='phys=6,10,10 mods=3,6,7' --set-quiccir-view-lay='lay-ppp2mpp=DCCSC3D,DCCSC3D lay-pmp2mmp=DCCSC3D,S1CLCSC3D lay-pmm2mmm=DCCSC3D,DCCSC3D' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize
