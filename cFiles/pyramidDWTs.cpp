/**
 * implimentation of the Haar wavelet transform and inverse using the MODWT method outlined in Percival
 * and Walden (2000) "Wavelet Methods for Time Series Analysis." (PW2000) See page 177 for the pseudocode.
 * 
 * compile with:
 * g++ -std=c++20 -openmp -shared -fPIC pyramidDWTs.cpp -o pyramidDWTs.so
 */

#include <vector>
#include <algorithm>
#include <cmath>

#include <iostream>
#include <stdio.h>


/**
 * @brief Structure to hold the MODWT wavelet filters and scaling filters for the Haar wavelet transform.
 * 
 * @param haarWaveletFilter The Haar wavelet filter coefficients
 * @param haarScalingFilter The Haar scaling filter coefficients
 */
struct MODWTFilters
{
    std::vector<double> waveletFilter;
    std::vector<double> scalingFilter;
};


/**
 * @brief Values for the Haar wavelet and scaling filters used in the MODWT transform.
 * See page 163 of Percival and Walden (2000) "Wavelet Methods for Time Series Analysis." for the 
 * definition of the MODWT filters, page 60 for the definition of the Haar wavelet filter, and equation
 * 75a for the definition of the scaling filter.
 */
const MODWTFilters haarFilters
{
    {0.5, -0.5}, // Haar wavelet filter
    {0.5, 0.5}  // Haar scaling filter
};


/**
 * @brief Values for the Daubechies 2 (db2) wavelet and scaling filters used in the MODWT transform.
 * See page 163 of Percival and Walden (2000) "Wavelet Methods for Time Series Analysis." for the
 * definition of the db2 MODWT filters, equation 59a for the definition of the db2 wavelet filter, and
 * equation 75a for the definition of the scaling filter.
 */
const MODWTFilters haarMODWTDb2Filters
{
    {(1-sqrt(3))/8, (-3+sqrt(3))/8, (3+sqrt(3))/8, (-1-sqrt(3))/8},
    {(1+sqrt(3))/8, (3+sqrt(3))/8, (3-sqrt(3))/8, (1-sqrt(3))/8}
};


/**
 * @brief Perform the MODWT transform on the input signal using the specified filters.
 * 
 * @param filters The MODWT filters to use for the transform
 * @param input The input signal to transform
 * @param output The output vector to hold the wavelet coefficients for each level of decomposition
 * @param intermediate The output vector to hold the final scaling coefficients after the transform
 */
void internalMODWTTransform
(
    const MODWTFilters& filters,
    const std::vector<double>& input,
    std::vector<std::vector<double>>& output,
    std::vector<double>& intermediate,
    int jMax,
    int N
)
{
    // get the length of the filters
    int L = filters.waveletFilter.size();
    // intitialise the current and next intermediate vectors
    std::vector<double> VCurrent = input; // set the initial scaling coefficients to the input signal
    std::vector<double> VNext(N, 0.0); // initialise the next intermediate vector
    // MODWT algorithm for J levels of decomposition
    // j loop for each level of decomposition (for j=1,...,J)
    for (int j = 1; j <= jMax; ++j)
    {
        //#pragma omp parallel for schedule(static)
        // t loop for each time point (for t=0,...,N-1)
        for (int t = 0; t < N; ++t)
        {
            // set k to t
            int k = t;
            // set W_{j,t} to h_0 * V_{j-1,k}
            output[j][t] = filters.waveletFilter[0] * VCurrent[k];
            // set V_{j,t} to g_0 * V_{j-1,k}
            VNext[t] = filters.scalingFilter[0] * VCurrent[k];
            // inner loop (for n=1,...,L-1)
            for (int n = 1; n < L; ++n)
            {
                // decrement k by 2^{j-1}
                k -= pow(2, j - 1);
                // if k<0, set k to k+N
                if (k < 0)
                    k += N;
                // increment W_{j,t} by h_n * V_{j-1,k}
                output[j][t] += filters.waveletFilter[n] * VCurrent[k];
                // increment V_{j,t} by g_n * V_{j-1,k}
                VNext[t] += filters.scalingFilter[n] * VCurrent[k];
            }
        }
        VCurrent = VNext;
    }
    intermediate = VCurrent; // set the final scaling coefficients to the intermediate vector
}


/**
 * @brief Perform the MODWT transform on the input signal using the specified filters.
 * 
 * @param signal The input signal to transform
 * @param signalLength The length of the input signal
 * @param transformType The type of transform to perform (0 for Haar, 1 for Db2)
 * @param waveletOutput The output array to hold the wavelet coefficients for each level of decomposition
 * @param scalingOutput The output array to hold the final scaling coefficients after the transform
 * @return 0 on success,
 * @return 1 on null pointer error,
 * @return 2 on invalid size error,
 * @return 3 on invalid transform type error
 */
extern "C" int MODWTTransform
(
    const double* signal,
    int signalLength,
    int transformType,
    double* waveletOutput,
    double* scalingOutput
)
{
    // basic input validation
    if (signal == nullptr || waveletOutput == nullptr || scalingOutput == nullptr){return 1;} // null pointer error
    if (signalLength <= 0){return 2;} // invalid size error
    // convert the input arrays to vectors for easier manipulation
    std::vector<double> signalVector(signal, signal + signalLength);
    std::vector<double> intermediateVector(signalLength, 0.0);
    // get the number of decomposition levels
    int jMax = floor(log2(signalLength)) - 1;
    // create a 2D vector to hold the wavelet coefficients for each level of decomposition
    std::vector<std::vector<double>> waveletOutputVector(
        jMax + 1,
        std::vector<double>(signalLength, 0.0)
    );
    // perform the MODWT transform using the correct filters based on the transformType
    if (transformType == 0)
    {
        internalMODWTTransform
        (
            haarFilters,
            signalVector,
            waveletOutputVector,
            intermediateVector,
            jMax,
            signalLength
        );
    }
    else if (transformType == 1)
    {
        internalMODWTTransform
        (
            haarMODWTDb2Filters,
            signalVector,
            waveletOutputVector,
            intermediateVector,
            jMax,
            signalLength
        );
    }
    else
    {
        return 3; // invalid transform type error
    }

    for (int j = 1; j <= jMax; ++j)
    {
        std::copy(
            waveletOutputVector[j].begin(),
            waveletOutputVector[j].end(),
            waveletOutput + ((j - 1) * signalLength)
        );
    }
    std::copy(
        intermediateVector.begin(),
        intermediateVector.end(),
        scalingOutput
    );

    return 0;
}


/**
 * @brief Perform the inverse MODWT transform on the wavelet coefficients to reconstruct the original signal.
 * 
 * @param filters The MODWT filters to use for the inverse transform
 * @param input The wavelet coefficients for each level of decomposition
 * @param output The output vector to hold the reconstructed signal
 * @param intermediate The input vector holding the final scaling coefficients after the transform
 */
void internalMODWTInverseTransform
(
    const MODWTFilters& filters,
    const std::vector<std::vector<double>>& input,
    std::vector<double>& output,
    std::vector<double>& intermediate,
    int levelCount,
    int N
)
{
    // get the length of the filters
    int L = filters.waveletFilter.size();
    // intitialise the current and previous intermediate vectors
    std::vector<double> VCurrent = intermediate; // set the initial scaling coefficients to the intermediate vector
    std::vector<double> VPrevious(N, 0.0);

    for (int j = levelCount - 1; j >= 0; --j)
    {
        // outer loop (for t=0,...,N-1)
        //#pragma omp parallel for schedule(static)
        for (int t = 0; t < N; ++t)
        {
            // set k to t
            int k = t;
            // set V_{j-1,t} to h_0 * W_{j,k} + g_0 * V_{j,k}
            VPrevious[t] = filters.waveletFilter[0] * input[j][k] + filters.scalingFilter[0] * VCurrent[k];
            // inner loop (for n=1,...,L-1)
            for (int n = 1; n < L; ++n)
            {
                // increment k by 2^j
                k += (1 << j);
                // if k>=N, set k to k-N
                if (k >= N)
                    k -= N;
                // increment V_{j-1,t} by h_n * W_{j,k} + g_n * V_{j,k}
                VPrevious[t] += filters.waveletFilter[n] * input[j][k] + filters.scalingFilter[n] * VCurrent[k];
            }
        }
        VCurrent = VPrevious; // set the current intermediate vector to the previous intermediate vector
    }
    output = VPrevious; // set the reconstructed signal to the previous intermediate vector
}


/**
 * @brief Perform the inverse MODWT transform on the wavelet coefficients to reconstruct the original signal.
 * 
 * @param waveletInput The wavelet coefficients for each level of decomposition
 * @param waveletInputLength The length of the wavelet input array
 * @param scalingInput The final scaling coefficients after the transform
 * @param transformType The type of transform to perform (0 for Haar, 1 for Db2)
 * @param output The output array to hold the reconstructed signal
 * @return 0 on success,
 * @return 1 on null pointer error,
 * @return 2 on invalid size error,
 * @return 3 on invalid transform type error
 */
extern "C" int MODWTInverseTransform
(
    const double* waveletInput,
    int waveletInputLength,
    int waveletLevelCount,
    const double* scalingInput,
    int transformType,
    double* output
)
{
    // basic input validation
    if (waveletInput == nullptr || scalingInput == nullptr || output == nullptr){return 1;} // null pointer error
    if (waveletInputLength <= 0 || waveletLevelCount <= 0){return 2;} // invalid size error
    if (waveletInputLength % waveletLevelCount != 0){return 2;} // invalid size/shape error

    // convert the input arrays to vectors for easier manipulation
    std::vector<double> waveletInputVector(waveletInput, waveletInput + waveletInputLength);
    int signalLength = waveletInputLength / waveletLevelCount;
    std::vector<double> scalingInputVector(scalingInput, scalingInput + signalLength);
    std::vector<double> outputVector(signalLength, 0.0);

    // create a 2D vector to hold the wavelet coefficients for each level of decomposition
    std::vector<std::vector<double>> waveletInput2DVector(
        waveletLevelCount,
        std::vector<double>(signalLength, 0.0)
    );
    for (int j = 0; j < waveletLevelCount; ++j)
    {
        std::copy(
            waveletInputVector.begin() + (j * signalLength),
            waveletInputVector.begin() + ((j + 1) * signalLength),
            waveletInput2DVector[j].begin()
        );
    }
    // perform the inverse MODWT transform using the correct filters based on the transformType
    if (transformType == 0)
    {
        internalMODWTInverseTransform
        (
            haarFilters,
            waveletInput2DVector,
            outputVector,
            scalingInputVector,
            waveletLevelCount,
            signalLength
        );
    }
    else if (transformType == 1)
    {
        internalMODWTInverseTransform
        (
            haarMODWTDb2Filters,
            waveletInput2DVector,
            outputVector,
            scalingInputVector,
            waveletLevelCount,
            signalLength
        );
    }
    else
    {
        return 3; // invalid transform type error
    }

    std::copy(
        outputVector.begin(),
        outputVector.end(),
        output
    );

    return 0;
}


#ifdef MODWT_STANDALONE_DEMO
int main()
{
    std::vector<double> signal = {1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0};
    int J = 2;

    std::vector<double> waveletOutput(signal.size() * J, 0.0);
    std::vector<double> scalingOutput(signal.size(), 0.0);

    int status = MODWTTransform
    (
        signal.data(),
        static_cast<int>(signal.size()),
        0, // Haar transform
        waveletOutput.data(),
        scalingOutput.data()
    );

    if (status != 0)
    {
        std::cerr << "Error performing MODWT transform: " << status << std::endl;
        return status;
    }

    for (int j = 0; j < J; ++j)
    {
        std::cout << "W" << (j + 1) << ": ";
        for (int t = 0; t < static_cast<int>(signal.size()); ++t)
        {
            std::cout << waveletOutput[j * signal.size() + t] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "V" << J << ": ";
    for (double val : scalingOutput)
    {
        std::cout << val << " ";
    }
    std::cout << "\n";

    return 0;
}
#endif

