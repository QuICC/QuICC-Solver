# Some notes on how to modify the QuICC code

## Adding a backward transform path

*The following will be replaced soon by a new version, in some models it already is. This section was written when developing the `AnelasticShellRTC` model.*

### Example: curl 

#### The backward transform

The curl is used in the computation of the nonlinear term in `SphericalSelfAdvection.cpp` or in the `MomentumKernel.cpp` of the `BoussinesqShellRTC` model.

The curl is defined via a transform path, operating directly on the toroidal-poloidal potentials, and is written in the file `ShellTransformSteps.cpp`, via the function `backwardCurl`.

The differential operators acting on each leg of the path such as:

    
    transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));

    transform.back().addEdge(Backward::Slaplr::id());

    transform.back().addEdge(Backward::OversinDphi::id());
    
    transform.back().addEdge(Backward::P::id(), FieldComponents::Physical::THETA, Arithmetics::Sub::id());
    
The operators `Slaplr` and `OversinDphi` are hashes, declared in the file `Components/Framework/cmake.d/setup/RegisterTransformBackward.cmake`. The hashes are 'place-holders', defined in files that can be found in  `Components/Framework/src/Transform/`. Namely:
- `OversinDphi` is defined in  `Components/Framework/src/Transform/ALegendreTransform.cpp`
- `Slaplr` (for the shell) is defined in `Components/Framework/src/Transform/DefaultShellChebyshevMap.cpp`

In these files one finds the actual projectors and integrators these operators refer to. For the purposes of the calculation of the curl, we have to look at Projectors.

For example:
- `OversinDphi` is defined through the `DivS1Dp` Projector, The relevant paths are:
    -   `Components/Transform/include/QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1Dp.hpp`
    - `Components/Transform/include/QuICC/Transform/Poly/ALegendre/Projector/DivS1Dp.hpp` (relevant?)
    -  `Components/Transform/src/Poly/ALegendre/Integrator/Base/DivS1Dp.cpp`

which lead to:
- `Components/Transform/include/QuICC/Transform/Poly/ALegendre/Projector/Base/DivS1.hpp`
- `Components/Transform/src/Poly/ALegendre/Projector/Base/DivS1.cpp`

and eventually to:
- `Components/Polynomial/include/QuICC/Polynomial/ALegendre/sin_1Plm.hpp` 
- `Components/Polynomial/include/QuICC/Polynomial/ALegendre/ALegendreBase.hpp`
- `Components/Polynomial/src/ALegendre/ALegendreBase.cpp`
where the actual recurrence relations seem to be implemented. 

There is a test for these operators, in `Components/Framework/TestSuite/Tests/Framework/TransformConfigurators/TransformTree/TorPolTreeTest.cpp`