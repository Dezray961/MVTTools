#ifndef HARRCOEFFICIENT_HPP
#define HARRCOEFFICIENT_HPP
#include <vector>

// MVTAnalysis class declaration
class HaarCoefficient
{
    public:
        // data members
        std::vector<double> rate;
        std::vector<double> time;
        std::vector<double> rateErr;
        std::vector<double> timeInBins;
        int lengthOfData;
        int kMax;
        std::vector<int> scaleSet;
        double deltaT;
        double coefficientValue;
        double coefficientVariance;
        int HaarLevel;
        std::vector<double> results;


        /**
         * Calculates the Haar coefficient for the given data and parameters.
         * @param rate A vector of the log of the count rate values.
         * @param time A vector of time values.
         * @param rateErr A vector of propagated rate error values.
         * @param timeInBins A vector of time values in bins.
         * @param lengthOfData The length of the data vectors.
         * @param scaleSet A vector of scale values to use for the Haar coefficient calculation.
         * @param kMax The maximum value of k to use for the Haar coefficient calculation.
         * @return A vector containing the results of the Haar coefficient calculation.
         */
        HaarCoefficient
        (
            std::vector<double> rate,
            std::vector<double> time,
            std::vector<double> rateErr,
            std::vector<double> timeInBins,
            int lengthOfData,
            std::vector<int> kSet,
            int kMax
        );


    private:
        /**
         * Calculates the mean of a given slice of data.
         * @param slice A vector of data values to calculate the mean of.
         * @return The mean of the given slice of data.
         */
        double sliceMean
        (
            std::vector<double> slice,
            int sliceSize
        );


        /**
         * Calculates the Haar coefficient for a given block of data.
         * @param block A vector of data values to calculate the Haar coefficient for.
         * @return The Haar coefficient for the given block of data.
         */
        double coefficient
        (
            std::vector<double> block,
            int halfBlockSize
        );


        /**
         * Calculates the variance of the Haar coefficient for a given block of data.
         * @param block A vector of data values to calculate the Haar coefficient variance for. This must be the squared propagated error of the log of the count rate values.
         * @return The variance of the Haar coefficient for the given block of data.
         */
        double coefficientVariance
        (
            std::vector<double> block,
            int halfBlockSize
        );


        /**
         * Calculates the time difference between the first and last time values in a given block of data.
         * @param timeInBinsBlock A vector of time values to calculate the time difference for.
         * @return The total time in a given block of data.
         */
        double tauIJ
        (
            int blockSize,
            int startIndex
        );


        /**
         * Converts a Haar coefficient value to a power value.
         * @param coefficientValue The Haar coefficient value to convert to a power value.
         * @param coefficientVariance The variance of the Haar coefficient value.
         * @return The power value corresponding to the given Haar coefficient value.
         */
        double convertCoefficientToPower
        (
            double coefficientValue,
            double coefficientVariance
        );
};

#endif // HARRCOEFFICIENT_HPP