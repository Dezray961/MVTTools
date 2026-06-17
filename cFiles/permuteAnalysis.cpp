/// include staements
#include "permuteAnalysis.hpp"
#include <math.h>
#include <numeric>
#include <algorithm>
#include "MVTClass.hpp"

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
    findDeltaTSet();
    this->kSetSize = kSet.size();
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


void PermuteAnalysis::findDeltaTSet()
{
    // itterate through the k set and calculate Δt for each k
    for (int k : kSet)
    {
        std::vector<double> differenceList; // list to store the differences between time values
        // calculate the differences between time values for the given k
        for (int i = k; i < lengthOfData; i++)
        {
            double difference = time[i] - time[i - k];
            differenceList.push_back(difference);
        }
        // calculate the mean of the differences
        double meanDeltaT = std::reduce(differenceList.begin(), differenceList.end(), 0.0) / (lengthOfData - k);
        this->deltaT.push_back(meanDeltaT);
        // calculate the standard deviation of the differences
        // create a new vector to store the differences between each difference and the mean
        std::vector<double> deltaTDifferenceList(differenceList.size());
        std::transform // applies a function to each element in the input range and stores the result in the output range
        (
            differenceList.begin(), // input range
            differenceList.end(), // input range
            deltaTDifferenceList.begin(), // output range
            [meanDeltaT](double difference) { return difference - meanDeltaT; } // lambda function
        );
        // calculate the sum of the squared differences
        double sumOfSquaredDifferences = std::inner_product(
            deltaTDifferenceList.begin(),
            deltaTDifferenceList.end(),
            deltaTDifferenceList.begin(),
            0.0
        );
        // calculate the standard deviation
        double standardDeviation = sqrt(sumOfSquaredDifferences / deltaTDifferenceList.size());
        // store the standard deviation in the deltaTError vector
        this->deltaTError.push_back(standardDeviation);
    }
};


std::vector<double> PermuteAnalysis::generatePermutation
(
    int index
)
{

}


std::vector<std::vector<double>> PermuteAnalysis::runMVTAnalysisOnPermutations()
{

}