// phys -> mod
// Nr x Nphi x Ntheta
// !type_uval = !quiccir.view<3x10x6xf64, "DCCSC3D">
// L x N x M  ( ..., radial, ...)
// !type_umod = !quiccir.view<6x2x6xcomplex<f64>, "DCCSC3D">

func.func @entry(%R: tensor<?x?x?xf64>) -> (tensor<?x?x?xcomplex<f64>>) {
  // R
  %R1 = quiccir.fr.int %R : tensor<?x?x?xf64> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 10 :i64, kind = "P"}
  %R1T = quiccir.transpose %R1 permutation = [2, 0, 1] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 11 :i64, kind = "P"}
  %R2 = quiccir.al.int %R1T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 12 :i64, kind = "P"}
  %R2T = quiccir.transpose %R2 permutation = [2, 0, 1] : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 13 :i64, kind = "P"}
  %R3 = quiccir.jw.int %R2T : tensor<?x?x?xcomplex<f64>> -> tensor<?x?x?xcomplex<f64>> attributes{implptr = 14 :i64, kind = "P"}
  return %R3 : tensor<?x?x?xcomplex<f64>>
}

// ./bin/quiccir-opt ../simple-3d-fwd.mlir --inline --quiccir-view-wrapper='dim-rets=7,3,6 dim-args=6,10,10 lay-args=DCCSC3D lay-rets=DCCSC3D' --inline --set-quiccir-dims='phys=6,10,10 mods=3,6,7' --set-quiccir-view-lay='lay-ppp2mpp=DCCSC3D,DCCSC3D lay-pmp2mmp=DCCSC3D,S1CLCSC3D lay-pmm2mmm=DCCSC3D,DCCSC3D' --convert-quiccir-to-call --quiccir-view-deallocation --lower-quiccir-alloc --canonicalize --finalize-quiccir-view --convert-func-to-llvm --canonicalize
