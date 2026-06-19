/// include staements
#include "permuteAnalysis.hpp"
#include <math.h>
#include <numeric>
#include <algorithm>
#include "HaarCoefficient.hpp"

/// public member functions
/// constructor
PermuteAnalysis::PermuteAnalysis
(
    double *rateArray,
    double *timeArray,
    double *rateErrArray,
    int lenghtOfData,
    int numberOfUniformKValues = 0
)
{
    this->lengthOfData = lenghtOfData;
    for (int i = 0; i < lenghtOfData; i++)
    {
        this->rate.push_back(rateArray[i]);
        this->time.push_back(timeArray[i]);
        this->rateErr.push_back(rateErrArray[i]);
    }
    findKMax();
    if (numberOfUniformKValues > 0)
    {
        generateUniformKSet(numberOfUniformKValues);
    }
    else
    {
        findKSet();
    }
    this->kSetSize = kSet.size();
    getTimeInBins();
}


/// private member functions
void PermuteAnalysis::findKMax()
{
    int exponent = floor(log2(lengthOfData));
    int powerOfTwo = pow(2, exponent);
    this->kMax = powerOfTwo;
};


void PermuteAnalysis::findKSet()
{
    for (int number = 1; number <= kMax; number *= 2)
    {
        this->kSet.push_back(number);
    }
};


void PermuteAnalysis::generateUniformKSet
(
    int numberOfUniformKValues
)
{
    /// find the logarithm of the maximum and minimum Δt values
    double logDeltaTMax = log10(time[lengthOfData - 1] - time[0]);
    double logDeltaTMin = log10(time[1] - time[0]);
    /// find the step size for the logarithm of Δt values
    double logDeltaTStep = (logDeltaTMax - logDeltaTMin) / (numberOfUniformKValues - 1);
    /// generate the uniform k set based on the logarithm of Δt values
    for (int i = 0; i < numberOfUniformKValues; i++)
    {
        double logDeltaT = logDeltaTMin + i * logDeltaTStep;
        double deltaT = pow(10, logDeltaT);
        int k = 0;
        for (int j = 1; j < lengthOfData; j++)
        {
            if (time[j] - time[0] >= deltaT)
            {
                k = j;
                break;
            }
        }
        if (k > 0 && k <= kMax)
        {
            this->kSet.push_back(k);
        }
    }
};


void PermuteAnalysis::getTimeInBins()
{
    for (int i = 0; i < lengthOfData - 1; i++)
    {
        timeInBins.push_back(time[i] - time[i + 1]);
    }
}



std::vector<std::vector<double>> PermuteAnalysis::shiftData
(
    int index
)
{
    std::vector<double> permutedRate;
    std::vector<double> permutedSigma;
    for (int i = 0; i < lengthOfData; i++)
    {
        int permutedIndex = (i + index) % lengthOfData;
        permutedRate.push_back(rate[permutedIndex]);
        permutedSigma.push_back(rateErr[permutedIndex]);
    }
    return {permutedRate, permutedSigma};
}


std::vector<double> PermuteAnalysis::analysePermutation
(
    int index
)
{
    /// generate a permutation of the dataset based on the given index
    std::vector<std::vector<double>> permutedData = shiftData(index);
    /// run the MVTAnalysis on the permuted dataset
    HaarCoefficient haarCoefficient(permutedData[0], time, permutedData[1], timeInBins, lengthOfData, kSet, kMax);
    return haarCoefficient.results;
}


std::vector<std::vector<double>> PermuteAnalysis::runMVTAnalysisOnPermutations()
{
    for (int shiftIndex = 0; shiftIndex < lengthOfData; shiftIndex++)
    {
        std::vector<double> permutationResults = analysePermutation(shiftIndex);
        
    }
}