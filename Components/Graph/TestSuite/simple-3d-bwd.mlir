// mod -> phys
// Nr x Nphi x Ntheta
// !type_uval = !quiccir.view<3x10x6xf64, "DCCSC3D">
// L x N x M  ( ..., radial, ...)
// !type_umod = !quiccir.view<6x2x6xxomplex<f64>, "DCCSC3D">

func.func @entry(%Pol: tensor<?x?x?xcomplex<f64>>) -> (tensor<?x?x?xf64>) {
  // Pol
  %Pol1 = quiccir.jw.prj %Pol : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 0 :i64}
  %Pol1T = quiccir.transpose %Pol1 permutation = [1, 2, 0] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 1 :i64}
  %Pol2 = quiccir.al.prj %Pol1T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 2 :i64}
  %Pol2T = quiccir.transpose %Pol2 permutation = [1, 2, 0] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 3 :i64}
  %Pol3 = quiccir.fr.prj %Pol2T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xf64> attributes{implptr = 4 :i64}
  return %Pol3 : tensor<?x?x?xf64>
}

// ./bin/quiccir-opt ../simple-3d-bwd.mlir --inline --quiccir-view-wrapper='dim-rets=7,3,6 dim-args=7,3,6 lay-args=DCCSC3D lay-rets=DCCSC3D' --inline --set-quiccir-dims='phys=6,10,10 mods=3,6,7' --set-quiccir-view-lay='lay-ppp2mpp=R_DCCSC3D,DCCSC3D lay-pmp2mmp=DCCSC3D,S1CLCSC3D lay-pmm2mmm=DCCSC3D,DCCSC3D' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize