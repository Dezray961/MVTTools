#include <vector> // dynamic arrays
#include <math.h> // mathematical functions
#include <numeric> // for std::accumulate
#include <algorithm> // for std::max
#include "HaarCoefficient.hpp" // include the header file for the HaarCoefficient class

#include <stdio.h> // for printf debugging


HaarCoefficient::HaarCoefficient
(
    const std::vector<double>& rate,
    const std::vector<double>& time,
    const std::vector<double>& rateErr,
    const std::vector<double>& timeInBins,
    int lengthOfData,
    std::vector<int> scaleSet,
    int shiftOffset
)
    : rate(rate),
      time(time),
      rateErr(rateErr),
      timeInBins(timeInBins)
{
    // initialize the data members with the provided arguments
    this->lengthOfData = lengthOfData;
    this->scaleSet = scaleSet;
    this->shiftOffset = shiftOffset;
    // reserve space for the results vector
    results.reserve(lengthOfData * std::max<std::size_t>(1, scaleSet.size()));
    // calculate the Haar coefficients and store the results
    findHaarCoefficients();
};


double HaarCoefficient::sliceMean
(
    const std::vector<double>& slice,
    int sliceSize
)
{
    double sum = std::accumulate(slice.begin(), slice.end(), 0.0);
    return sum / sliceSize;
}


double HaarCoefficient::findCoefficient
(
    int startIndex,
    int halfBlockSize
)
{
    double leftMean = 0.0;
    double rightMean = 0.0;
    for (int i = 0; i < halfBlockSize; i++)
    {
        int leftIndex = (startIndex + i + shiftOffset) % lengthOfData;
        int rightIndex = (startIndex + halfBlockSize + i + shiftOffset) % lengthOfData;
        leftMean += rate[leftIndex];
        rightMean += rate[rightIndex];
    }
    leftMean /= halfBlockSize;
    rightMean /= halfBlockSize;
    return rightMean - leftMean;
}


double HaarCoefficient::findCoefficientVariance
(
    int startIndex,
    int halfBlockSize
)
{
    double leftVariance = 0.0;
    double rightVariance = 0.0;
    for (int i = 0; i < halfBlockSize; i++)
    {
        int leftIndex = (startIndex + i + shiftOffset) % lengthOfData;
        int rightIndex = (startIndex + halfBlockSize + i + shiftOffset) % lengthOfData;
        leftVariance += rateErr[leftIndex];
        rightVariance += rateErr[rightIndex];
    }
    leftVariance /= (halfBlockSize * halfBlockSize);
    rightVariance /= (halfBlockSize * halfBlockSize);
    return leftVariance + rightVariance;
}


double HaarCoefficient::tauIJ
(
    int blockSize,
    int startIndex
)
{
    /// calculate the total time in a given block of data by summing the timeInBins vector
    return std::accumulate(timeInBins.begin() + startIndex, timeInBins.begin() + startIndex + blockSize - 1, 0.0);
}


double HaarCoefficient::convertCoefficientToPower
(
    double coefficientValue,
    double coefficientVariance
)
{
    return (coefficientValue * coefficientValue) - coefficientVariance;
}


void HaarCoefficient::findHaarCoefficients()
{
    std::size_t estimatedResults = 0;
    for (int blockSize : scaleSet)
    {
        if (blockSize > 1 && blockSize <= lengthOfData)
        {
            estimatedResults += static_cast<std::size_t>(lengthOfData - blockSize + 1);
        }
    }
    results.reserve(estimatedResults);

    for (int blockSize : scaleSet)
    {
        if (blockSize < 2 || blockSize > lengthOfData || (blockSize % 2) != 0)
        {
            continue;
        }
        int halfBlockSize = blockSize / 2;
        for (int index = 0; index + blockSize <= lengthOfData; index++)
        {
            /// calculate the Haar coefficient, its variance, and the power for the current block
            double coefficientValue = findCoefficient(index, halfBlockSize);
            double variance = findCoefficientVariance(index, halfBlockSize);
            double power = convertCoefficientToPower(coefficientValue, variance);
            double logTauIJ = log10(tauIJ(blockSize, index));
            /// store the results in the results vector
            results.emplace_back(std::initializer_list<double>{logTauIJ, power});
        }
    }
}    
