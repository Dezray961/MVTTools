#include <vector>
#include <math.h>
#include <numeric>
#include <algorithm>
#include "HaarCoefficient.hpp"


HaarCoefficient::HaarCoefficient
(
    std::vector<double> rate,
    std::vector<double> time,
    std::vector<double> rateErr,
    std::vector<double> timeInBins,
    int lengthOfData,
    std::vector<int> scaleSet
)
{
    this->rate = rate;
    this->time = time;
    this->rateErr = rateErr;
    this->timeInBins = timeInBins;
    this->lengthOfData = lengthOfData;
    this->scaleSet = scaleSet;
    results.reserve(lengthOfData * std::max<std::size_t>(1, scaleSet.size()));
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
    auto leftBegin = rate.begin() + startIndex;
    auto rightBegin = leftBegin + halfBlockSize;
    auto rightEnd = rightBegin + halfBlockSize;
    double leftMean = std::accumulate(leftBegin, rightBegin, 0.0) / halfBlockSize;
    double rightMean = std::accumulate(rightBegin, rightEnd, 0.0) / halfBlockSize;
    return rightMean - leftMean;
}


double HaarCoefficient::findCoefficientVariance
(
    int startIndex,
    int halfBlockSize
)
{
    /**
     *  because the input block is the squared propagated error of the log of the count rate values, we can calculate 
     * the variance of the Haar coefficient using the sliceMean and just dividing by the halfBlockSize. This is the same
     * as 1/n^2 * sum(sigma^2) 
     */
    auto leftBegin = rateErr.begin() + startIndex;
    auto rightBegin = leftBegin + halfBlockSize;
    auto rightEnd = rightBegin + halfBlockSize;
    double leftVariance = std::accumulate(leftBegin, rightBegin, 0.0) / (halfBlockSize * halfBlockSize);
    double rightVariance = std::accumulate(rightBegin, rightEnd, 0.0) / (halfBlockSize * halfBlockSize);
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
            results.push_back({logTauIJ, power});
        }
    }
}    
