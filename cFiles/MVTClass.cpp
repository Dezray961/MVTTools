#include <vector>
#include <math.h>
#include <numeric>
#include "MVTClass.hpp"


namespace
{
    constexpr int kParallelThreshold = 64;
}

// constructor
MVTAnalysis::MVTAnalysis
(
    std::vector<double> rateArray,
    std::vector<double> timeArray,
    std::vector<double> rateErrArray,
    int lenghtOfData,
    std::vector<int> kSet,
    int kMax
)
{
    this->lenghtOfData = lenghtOfData;
    this->rate = rateArray;
    this->time = timeArray;
    this->rateErr = rateErrArray;
    this->kSet = kSet;
    findVTForAllK();
};
        

// private member functions
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
void MVTAnalysis::findVTForAllK()
{
    for (int k : kSet)
    {
        this->VTSet.push_back(findVT(k));
    }
};

