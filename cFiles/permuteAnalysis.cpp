/// include staements
#include "permuteAnalysis.hpp"
#include <math.h>
#include <numeric>
#include <algorithm>
#include <limits>
#include "HaarCoefficient.hpp"

/// public member functions
/// constructor
PermuteAnalysis::PermuteAnalysis
(
    double *rateArray,
    double *timeArray,
    double *rateErrArray,
    int lenghtOfData,
    int numberOfTimeBins
)
{
    this->lengthOfData = lenghtOfData;
    this->numberOfTimeBins = numberOfTimeBins;
    this->rate.reserve(lenghtOfData);
    this->time.reserve(lenghtOfData);
    this->rateErr.reserve(lenghtOfData);
    for (int i = 0; i < lenghtOfData; i++)
    {
        this->rate.push_back(rateArray[i]);
        this->time.push_back(timeArray[i]);
        this->rateErr.push_back(rateErrArray[i]);
    }
    findKMax();
    findKSet();
    this->kSetSize = kSet.size();
    getTimeInBins();
    generateOutputBins();
    this->powerSetSums.assign(numberOfTimeBins, 0.0);
    this->powerSetSumSquares.assign(numberOfTimeBins, 0.0);
    this->powerSetCounts.assign(numberOfTimeBins, 0);
    this->totalPermutations = lengthOfData;    
}


void PermuteAnalysis::runAnalysis()
{
    permutationsCompleted.store(0, std::memory_order_relaxed);
    analysisComplete.store(false, std::memory_order_relaxed);
    runMVTAnalysisOnPermutations();
    analysisComplete.store(true, std::memory_order_relaxed);
}


int PermuteAnalysis::getPermutationsCompleted() const
{
    return permutationsCompleted.load(std::memory_order_relaxed);
}


int PermuteAnalysis::getTotalPermutations() const
{
    return totalPermutations;
}


bool PermuteAnalysis::isAnalysisComplete() const
{
    return analysisComplete.load(std::memory_order_relaxed);
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
    this->kSet.clear();
    for (int number = 1; number <= kMax; number *= 2)
    {
        this->kSet.push_back(number);
    }
};


void PermuteAnalysis::generateOutputBins()
{
    /// the Haar output is indexed by tau values derived from the bin widths,
    /// so the output bins must live on the same scale.
    double minTau = std::numeric_limits<double>::infinity();
    double maxTau = 0.0;
    for (double deltaT : timeInBins)
    {
        if (deltaT > 0.0)
        {
            minTau = std::min(minTau, deltaT);
            maxTau += deltaT;
        }
    }
    if (!std::isfinite(minTau) || minTau <= 0.0 || maxTau <= minTau)
    {
        minTau = 1.0;
        maxTau = 2.0;
    }
    double logMinTime = log10(minTau);
    double logMaxTime = log10(maxTau);
    /// generate a vector of evenly spaced values between the log of the minimum and maximum tau values
    double binWidth = (logMaxTime - logMinTime) / numberOfTimeBins;
    this->logBinEdges.clear();
    this->logBinCenters.clear();
    this->logBinEdges.reserve(numberOfTimeBins + 1);
    this->logBinCenters.reserve(numberOfTimeBins);
    /// generate the logBinEdges and logBinCenters vectors this is the 
    for (int i = 0; i <= numberOfTimeBins; i++)
    {
        logBinEdges.push_back(logMinTime + i * binWidth);
    }
    /// generate the logBinCenters vector
    for (int i = 0; i < numberOfTimeBins; i++)
    {
        logBinCenters.push_back((logBinEdges[i] + logBinEdges[i + 1]) / 2);
    }
}



void PermuteAnalysis::getTimeInBins()
{
    timeInBins.clear();
    timeInBins.reserve(lengthOfData - 1);
    for (int i = 0; i < lengthOfData - 1; i++)
    {
        timeInBins.push_back(std::abs(time[i + 1] - time[i]));
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


std::vector<std::vector<double>> PermuteAnalysis::analysePermutation
(
    int index
)
{
    /// generate a permutation of the dataset based on the given index
    std::vector<std::vector<double>> permutedData = shiftData(index);
    /// run the MVTAnalysis on the permuted dataset
    HaarCoefficient haarCoefficient(permutedData[0], time, permutedData[1], timeInBins, lengthOfData, kSet);
    return haarCoefficient.results;
}


void PermuteAnalysis::binResults
(
    const std::vector<std::vector<double>>& results,
    std::vector<double>& powerSums,
    std::vector<double>& powerSumSquares,
    std::vector<int>& powerCounts
) const
{
    /// find the correct time bin for each result and add the result to the corresponding bin
    for (std::size_t i = 0; i < results.size(); i++)
    {
        /// find the current logTauIJ and power values
        double logTauIJ = results[i][0];
        double power = results[i][1];
        /// find the correct bin for the current logTauIJ value
        for (std::size_t j = 0; j + 1 < logBinEdges.size(); j++)
        {
            if (logTauIJ >= logBinEdges[j] && logTauIJ < logBinEdges[j + 1])
            {
                /// check if this is a burst or pre-burst analysis
                if (!burst)
                {
                    powerSums[j] += power;
                    powerSumSquares[j] += power * power;
                    powerCounts[j] += 1;
                }
                else
                {
                 /// only add the power value to the bin if it is greater than 3 times the pre-burst standard deviation 
                    if (power > 3 * preBurstPowerSetStdDev[j])
                    {
                        powerSums[j] += power;
                        powerSumSquares[j] += power * power;
                        powerCounts[j] += 1;
                    }
                }
                break;
            }
        }
    }
}


void PermuteAnalysis::findAveragePowerInBins()
{
    powerSetAverages.clear();
    powerSetAverages.reserve(powerSetCounts.size());
    for (int i = 0; i < powerSetCounts.size(); i++)
    {
        double average = 0.0;
        if (powerSetCounts[i] > 0)
        {
            average = powerSetSums[i] / powerSetCounts[i];
        }
        powerSetAverages.push_back(average);
    }
}


void PermuteAnalysis::findStdDevOfPowerInBins()
{
    powerSetStdDevs.clear();
    powerSetStdDevs.reserve(powerSetCounts.size());
    for (int i = 0; i < powerSetCounts.size(); i++)
    {
        double mean = powerSetAverages[i];
        double stdDev = 0.0;
        if (powerSetCounts[i] > 0)
        {
            double variance = (powerSetSumSquares[i] / powerSetCounts[i]) - (mean * mean);
            if (variance < 0.0)
            {
                variance = 0.0;
            }
            stdDev = sqrt(variance);
        }
        powerSetStdDevs.push_back(stdDev);
    }
}



void PermuteAnalysis::runMVTAnalysisOnPermutations()
{
    std::fill(powerSetSums.begin(), powerSetSums.end(), 0.0);
    std::fill(powerSetSumSquares.begin(), powerSetSumSquares.end(), 0.0);
    std::fill(powerSetCounts.begin(), powerSetCounts.end(), 0);

    std::vector<double> globalPowerSums(numberOfTimeBins, 0.0);
    std::vector<double> globalPowerSumSquares(numberOfTimeBins, 0.0);
    std::vector<int> globalPowerCounts(numberOfTimeBins, 0);

    #pragma omp parallel
    {
        std::vector<double> localPowerSums(numberOfTimeBins, 0.0);
        std::vector<double> localPowerSumSquares(numberOfTimeBins, 0.0);
        std::vector<int> localPowerCounts(numberOfTimeBins, 0);

        #pragma omp for schedule(dynamic)
        for (int shiftIndex = 0; shiftIndex < lengthOfData; shiftIndex++)
        {
            std::vector<std::vector<double>> permutationResults = analysePermutation(shiftIndex);
            binResults(permutationResults, localPowerSums, localPowerSumSquares, localPowerCounts);
            permutationsCompleted.fetch_add(1, std::memory_order_relaxed);
        }

        #pragma omp critical
        {
            for (int i = 0; i < numberOfTimeBins; i++)
            {
                globalPowerSums[i] += localPowerSums[i];
                globalPowerSumSquares[i] += localPowerSumSquares[i];
                globalPowerCounts[i] += localPowerCounts[i];
            }
        }
    }

    powerSetSums = std::move(globalPowerSums);
    powerSetSumSquares = std::move(globalPowerSumSquares);
    powerSetCounts = std::move(globalPowerCounts);
    findAveragePowerInBins();
    findStdDevOfPowerInBins();
}