#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <numeric>
#include <string>
#include "MVTClass.hpp"
#include <algorithm>

namespace
{
    constexpr int kParallelThreshold = 64;
}

// constructor
MVTAnalysis::MVTAnalysis
(
    double *rateArray,
    double *timeArray,
    double *rateErrArray,
    int lenghtOfData,
    bool generateEvenKSetFlag
)
{
    this->lenghtOfData = lenghtOfData;
    for (int i = 0; i < lenghtOfData; i++)
    {
        this->rate.push_back(rateArray[i]);
        this->time.push_back(timeArray[i]);
        this->rateErr.push_back(rateErrArray[i]);
    }
    findKMax();
    if (generateEvenKSetFlag)
    {
        generateEvenKSet();
    }
    else
    {
        findKSet();
    }
    findDeltaTSet();
    this->kSetSize = kSet.size();
    this->VTSet = findVTForAllK();
    this->VTSetArray = getVTSetAsArray();

    localAverageOverKBinsForDataSet(4);
};
        

// private member functions
// function to find kMax
void MVTAnalysis::findKMax()
{
    int exponent = floor(log2(lenghtOfData));
    int powerOfTwo = pow(2, exponent);
    this->kMax = powerOfTwo;
};


// function to find the k set
void MVTAnalysis::findKSet()
{
    this->kSet = std::vector<int>();
    int number = 1;
    while (number <= kMax)
    {
        this->kSet.push_back(number);
        number *= 2;
    }
};


// function to generate a kSet if the user wants all even numbers between 1 and kMax
void MVTAnalysis::generateEvenKSet()
{
    this->kSet = std::vector<int>();
    for (int number = 2; number <= kMax; number += 2)
    {
        this->kSet.push_back(number);
    }
};


// function to find Δt for a given k
void MVTAnalysis::findDeltaTSet()
{
    // itterate through the k set and calculate Δt for each k
    for (int k : kSet)
    {
        std::vector<double> differenceList; // list to store the differences between time values
        // calculate the differences between time values for the given k
        for (int i = k; i < lenghtOfData; i++)
        {
            double difference = time[i] - time[i - k];
            differenceList.push_back(difference);
        }
        // calculate the mean of the differences
        double meanDeltaT = std::reduce(differenceList.begin(), differenceList.end(), 0.0) / (lenghtOfData - k);
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
}


// function to find the local average over k samples
double MVTAnalysis::localAverageOverKBins
(
    int index,
    int k
)
{
    std::vector<double> summationList;
    for (int i = 0; i < k; i++)
    {
        int interiorIndex = index - i;
        double interiorRate = rate[interiorIndex];
        summationList.push_back(rate[interiorIndex]);
    }
    return std::reduce(summationList.begin(), summationList.end(), 0.0) / k;
};


// function to find the local average over k samples for the entire data set
std::vector<double> MVTAnalysis::localAverageOverKBinsForDataSet
(
    int k
)
{
    const int dataPointCount = lenghtOfData - k;
    std::vector<double> localAverages(dataPointCount);
    if (dataPointCount <= 0)
    {
        return localAverages;
    }

    const bool useParallel = k >= kParallelThreshold;

#ifdef _OPENMP
#pragma omp parallel for if(useParallel) schedule(static)
#endif
    for (int i = 0; i < dataPointCount; i++)
    {
        localAverages[i] = localAverageOverKBins(i + k, k);
    }
    return localAverages;
};

// function to find the square of the difference between local averages
std::vector<double> MVTAnalysis::squaredDifference
(
    int k
)
{
    std::vector<double> localAverages = localAverageOverKBinsForDataSet(k);
    const int dataPointCount = lenghtOfData - k;
    std::vector<double> squaredDifferences(dataPointCount);
    if (dataPointCount <= 0)
    {
        return squaredDifferences;
    }

    const bool useParallel = k >= kParallelThreshold;

#ifdef _OPENMP
#pragma omp parallel for if(useParallel) schedule(static)
#endif
    for (int i = 0; i < dataPointCount; i++)
    {
        double difference = rate[i + k] - localAverages[i];
        squaredDifferences[i] = difference * difference;
    }
    return squaredDifferences;
};


// function to find the VT from the square of the difference between local averages
double MVTAnalysis::findVT
(
    int k
)
{
    std::vector<double> squaredDifferences = squaredDifference(k);
    double sumOfSquaredDifferences = std::reduce(squaredDifferences.begin(), squaredDifferences.end(), 0.0);
    return sqrt(sumOfSquaredDifferences / (lenghtOfData - k));
};


// function to find the VT for all k in the k set
std::vector<double> MVTAnalysis::findVTForAllK()
{
    std::vector<double> VTSet;
    for (int k : kSet)
    {
        VTSet.push_back(findVT(k));
    }
    return VTSet;
};


// function to conbvert the VT set to an array of doubles
double* MVTAnalysis::getVTSetAsArray()
{
    int size = VTSet.size();
    double* VTSetArray = new double[size];
    for (int i = 0; i < size; i++)
    {
        VTSetArray[i] = VTSet[i];
    }
    return VTSetArray;
};


