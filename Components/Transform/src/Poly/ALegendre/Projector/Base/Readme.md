# Projector operators





## Tests

The tests are performed as follows:

1. Create a `build` folder in the `QuICC` root directory. It is best to call it exactly just `build`. Not `build_test` or anything like that. Unless one were to modify the Mathematica script to place the references in a specific build folder (see below)
2. Add the new operator to `Components/Transform/TestSuite/Tests/Transform/ALegendre/Projector/CMakeLists.txt`.

3. Generate reference with the mathematica notebook `Components/Transform/TestSuite/Mathematica/ALegendreReference.nb`
4. Build the test
5. Run the test

See below for details on each step:


### Mathematica reference

The easiest thing to do is to modify an existing reference. Some caveats/things to consider:

- Be sure that the reference is placed in the `build` directory. Some Mathematica ref files will not have a general path generation but it's hardcoded into the notebook and needs to be changed. 

    This is not currently (January 2025) a problem for the notebook being discussed now

- It is expensive to generate tests for all operators implemented in the notebook. Unless that is what you want, you might want to comment out operators that are not to be tested (comment out these lines in the `generateProjectors` Module definition).

- For the projection operator, only run the `generateProjectors` function. Not the whole notebook. You don't need, for example, `generateIntegrators` for this example.




### Build and run the test

In the `build` directory (or wherever you put the tests):

-       cmake ..  -DQUICC_TESTSUITE_TRANSFORM=ON

    Furthermore...supposedly configuring as 

            cmake ..  -DQUICC_TESTSUITE_TRANSFORM=ON -DQUICC_TESTSUITE_Llm1D1=ON

    is faster, but I haven't seen a difference

-       make -j 6 TransformALegendreTests

-       ctest -R Llm1D1 --output-on-failure




#### Some troubleshooting
