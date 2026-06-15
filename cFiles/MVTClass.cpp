#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <numeric>
#include <string>
#include "MVTClass.hpp"

// constructor
MVTAnalysis::MVTAnalysis
(
    double *rateArray,
    double *timeArray,
    double *rateErrArray,
    int lenghtOfData
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
    findKSet();
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
    return std::accumulate(summationList.begin(), summationList.end(), 0.0) / k;
};


// function to find the local average over k samples for the entire data set
std::vector<double> MVTAnalysis::localAverageOverKBinsForDataSet
(
    int k
)
{
    std::vector<double> localAverages;
    for (int i = k; i < lenghtOfData; i++)
    {
        double localAverage = localAverageOverKBins(i, k);
        localAverages.push_back(localAverage);
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
    std::vector<double> squaredDifferences;
    for (int i = k; i < lenghtOfData; i++)
    {
        double difference = rate[i] - localAverages[i - k];
        squaredDifferences.push_back(pow(difference, 2));
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
    double sumOfSquaredDifferences = std::accumulate(squaredDifferences.begin(), squaredDifferences.end(), 0.0);
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


