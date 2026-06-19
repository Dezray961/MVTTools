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
    std::vector<int> scaleSet,
    int kMax
)
{
    this->rate = rate;
    this->time = time;
    this->rateErr = rateErr;
    this->timeInBins = timeInBins;
    this->lengthOfData = lengthOfData;
    this->scaleSet = scaleSet;
    this->kMax = kMax;
};


double HaarCoefficient::sliceMean
(
    std::vector<double> slice,
    int sliceSize
)
{
    double sum = std::accumulate(slice.begin(), slice.end(), 0.0);
    return sum / sliceSize;
}


double HaarCoefficient::coefficient
(
    std::vector<double> block,
    int halfBlockSize
)
{
    std::vector<double> leftSlice(block.begin(), block.begin() + halfBlockSize);
    std::vector<double> rightSlice(block.begin() + halfBlockSize, block.end());
    double leftMean = sliceMean(leftSlice, halfBlockSize);
    double rightMean = sliceMean(rightSlice, halfBlockSize);
    return rightMean - leftMean;
}


double HaarCoefficient::coefficientVariance
(
    std::vector<double> block,
    int halfBlockSize
)
{
    std::vector<double> leftSlice(block.begin(), block.begin() + halfBlockSize);
    std::vector<double> rightSlice(block.begin() + halfBlockSize, block.end());
    /**
     *  because the input block is the squared propagated error of the log of the count rate values, we can calculate 
     * the variance of the Haar coefficient using the sliceMean and just dividing by the halfBlockSize. This is the same
     * as 1/n^2 * sum(sigma^2) 
     */
    double leftVariance = sliceMean(leftSlice, halfBlockSize) / halfBlockSize;
    double rightVariance = sliceMean(rightSlice, halfBlockSize) / halfBlockSize;
    return leftVariance + rightVariance;
}


double HaarCoefficient::tauIJ
(
    int blockSize,
    int startIndex
)
{
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


