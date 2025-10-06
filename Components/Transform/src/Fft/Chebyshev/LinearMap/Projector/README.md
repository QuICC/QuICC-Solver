# Projector operators

These operators perform a projection from one spectral field to another spectral field via a given operator involving some spatial derivative.

Consider a field $v$ whose spectral coefficients are $\hat{\bf v} = (\hat{v}_0, \hat{v}_1, \dots, \hat{v}_{N-1})^T$ and say we want to act on this field via the operator ```D2Y1```. This operator is meant to perform the following operation:
$$
\texttt{D2Y1} \rightarrow \partial^2_r(r v),
$$
in spectral space. Expanding $v$ in a Chebyshev basis:
$$
\partial^2_r(r v) = \sum_n t_n \hat{v}_n\partial^2_r\left(r T_n\right),
$$
where $t_n$ is a factor that takes care of the different normalisation of $T_0$ with respect of $T_n$ with $n>0$.

The ```D2Y1``` itself converts the $\hat{\bf v}$ into a set of $\hat{\bf w}$ so that:

$$
\sum_n t_n \hat{v}_n\partial^2_r\left(r T_n\right) = \sum_n t_n \hat{w}_n T_n.
$$

As a shorthand notation, we could write:
$$
w = \texttt{D2Y1} \cdot v.
$$


# Adding a new operator
The procedure to follow to define and test a new operator is:

1. Implement new operator (for example copying an existing operator and modifying it)
2. Add the operator to `Components/Transform/src/Fft/Chebyshev/LinearMap/Projector/CMakeList.txt`
3. Test the operator

See below for some details.

# Example: D2Y1 operator

The ```D2Y1``` operator was created modifying the (slightly simpler) ```D1Y1```operator. See comments in `D2Y1.cpp` and `D2Y1.hpp`. The subtelty is that the $w$ coefficients are not obtained directly from the differential operators, but are obtained as follows:
$$
w = \texttt{D2Y1} \cdot v \quad \Rightarrow \quad \texttt{I2} \cdot w = \texttt{Y1} \cdot v 
$$
and then solving for $w$. This is possible because $\texttt{I2}\cdot \texttt{D2}$ is an identity. Care must be taken in setting the correct sizes of the matrices, but it makes mathematical sense: we'd rather solve a linear system than calculate derivatives via recurrence relations such as:
$
(1-x^2) T'_n = -n x T_n + n T_{n-1},;
$
because then we'd have to deal with a division by $1-x^2$ (which, numerically, you'd like to avoid).



## Test the new D2Y1 operator

The tests are performed as follows:

1. Create a `build` folder in the `QuICC` root directory. It is best to call it exactly just `build`. Not `build_test` or anything like that. Unless one were to modify the Mathematica script to place the references in a specific build folder (see below)
2. Add the new operator to `Components/Transform/TestSuite/Tests/Transform/Chebyshev/LinearMap/Projector/CMakeLists.txt`.

    Remeber to add the ULP and MPULP for the new test. This is essentially the maximum allowed numerical difference between the reference and the test values.

    **If the test fails because of ULP issues, this value can be set to a higher value to make the test pass**

3. Generate reference with the mathematica notebook `Components/Transform/TestSuite/Mathematica/ChebyshevReference.nb`
4. Build the test
5. Run the test

See below for details on each step:



### Mathematica reference

The easiest thing to do is to modify an existing reference. Some caveats/things to consider:

- Be sure that the reference is placed in the `build` directory. Some Mathematica ref files will not have a general path generation but it's hardcoded into the notebook and needs to be changed. 

    This is not currently (January 2025) a problem for the notebook being discussed now

- It is expensive to generate tests for all operators implemented in the notebook. Unless that is what you want, you might want to comment out operators that are not to be tested (comment out these lines in the `generateProjectors` Module definition).

- For the `D2Y1` operator, only run the `generateProjectors` function. Not the whole notebook. You don't need, for example, `generateIntegrators` for this example.




### Build and run the test

In the `build` directory (or wherever you put the tests):

-       cmake ..  -DQUICC_TESTSUITE_TRANSFORM=ON

    Furthermore...supposedly configuring as 

            cmake ..  -DQUICC_TESTSUITE_TRANSFORM=ON -DQUICC_TESTSUITE_D2Y1=ON

    is faster, but I haven't seen a difference

-       make -j 6 TransformChebyshevTests

-       ctest -R D2Y1 --output-on-failure




#### Some troubleshooting

- If the test fails because of an error that looks like this:

        /QuICC/QuICC/Components/Transform/TestSuite/include/QuICC/TestSuite/Transform/TesterBase.hpp:428: FAILED:
        CHECK( std::get<0>(err) )
        with expansion:
        false
        with messages:
        type: projector
        id: 2
        n: 0
        position: 12 / 39
        refData: 2.3597691173755408e+01
        checked normal value
        outData: 2.3597691173753361e+01
        max ulp: 250
        measured ulp: 390.547

        ===============================================================================
        test cases:   1 |   0 passed | 1 failed
        assertions: 365 | 358 passed | 7 failed



        0% tests passed, 1 tests failed out of 1

        Total Test time (real) =   0.10 sec

        The following tests FAILED:
        5 - TransformChebyshevTests_D2Y1_projector_ulp250 (Failed)
        Errors while running CTest


    You can set the ULP in `Components/Transform/TestSuite/Tests/Transform/Chebyshev/LinearMap/Projector/CMakeLists.txt` to a higher value. This obviously makes sense if the test generally passes, except for a few values (that do not fail by much). 
    
    The rule is then to set the ULP to 10% higher than needed and round up to closest 5.







# Mathematica reference: calculation and output

For these operators the tests:

- take as input a given spectra $\hat{v}$ : a `minSpectrum` (only $T_1$ coefficient is nonzero), a `maxSpectrum` (only $T_n$ coefficients for the max $n$ is nonzero) and a `unitSpectrum` (each coefficient is $1-2\textrm{i}$). Each spectra is identified by the id of the output files.

- gives in output the value of the projected $w$ in physical space, evaluated at the quadrature points.

## Calculation

The analytical forms of the operators involves the functions, e.g., `dyChebyshev` and similar. These are the mathematical forms of the operators, expressed with $x$ as the spatial variable. To convert, use:
$$
r = ax+b.
$$
From which:
$$
\partial_x = a \ \partial_r;\quad dr = a \ dx
$$

From which, for example, the form of `dyChebyshev` is derived:
$$
\partial_r (r T_n) = T_n +(ax+b) a^{-1} T_n'
$$
as in the Mathematica notebook.

## Output

The output files are, for each operator:
- `*_meta*`, which contains, in order, outprec, nN, Ng, (-a+b), (a+b), {ls}

- `*_in*`, which contains the input spectral coefficients, with their real and imaginary part separated

        inData = Join[Re[spec],Im[spec]]

- `*_ref*`, the values of the output, in physical space, also ordered by their real and imaginary part. One column for each of the ls. Along each column we have the real parts on the grid points, and then the imaginary parts on the grid points.

