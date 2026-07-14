#ifndef HARRCOEFFICIENT_HPP
#define HARRCOEFFICIENT_HPP
#include <vector>

/**
 * The HaarCoefficient class calculates the Haar coefficients for a given dataset. It takes in vectors of rate, time, rate error, and time in bins, as well as the length of the data and a set of scales to use for the calculation. The class provides methods to calculate the Haar coefficient, its variance, and the power for each block of data, and stores the results in a vector.
 */
class HaarCoefficient
{
    public:
        // data members
        const std::vector<double>& rate;
        const std::vector<double>& time;
        const std::vector<double>& rateErr;
        const std::vector<double>& timeInBins;
        int lengthOfData;
        std::vector<int> scaleSet;
        int shiftOffset;
        double deltaT;
        double coefficientValue;
        double coefficientVariance;
        std::vector<std::vector<double>> results;


        /**
         * Calculates the Haar coefficient for the given data and parameters.
         * @param rate A vector of the log of the count rate values.
         * @param time A vector of time values.
         * @param rateErr A vector of propagated rate error values.
         * @param timeInBins A vector of time values in bins.
         * @param lengthOfData The length of the data vectors.
         * @param scaleSet A vector of scale values to use for the Haar coefficient calculation.
         * @return A vector containing the results of the Haar coefficient calculation.
         */
        HaarCoefficient
        (
            const std::vector<double>& rate,
            const std::vector<double>& time,
            const std::vector<double>& rateErr,
            const std::vector<double>& timeInBins,
            int lengthOfData,
            std::vector<int> scaleSet,
            int shiftOffset
        );


    private:
        /**
         * Calculates the mean of a given slice of data.
         * @param slice A vector of data values to calculate the mean of.
         * @return The mean of the given slice of data.
         */
        double sliceMean
        (
            const std::vector<double>& slice,
            int sliceSize
        );


        /**
         * Calculates the Haar coefficient for a given block of data.
         * @param block A vector of data values to calculate the Haar coefficient for.
         * @return The Haar coefficient for the given block of data.
         */
        double findCoefficient
        (
            int startIndex,
            int halfBlockSize
        );


        /**
         * Calculates the variance of the Haar coefficient for a given block of data.
         * @param block A vector of data values to calculate the Haar coefficient variance for. This must be the squared propagated error of the log of the count rate values.
         * @return The variance of the Haar coefficient for the given block of data.
         */
        double findCoefficientVariance
        (
            int startIndex,
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


        /**
         * Calculates the Haar coefficients for the given data and parameters.
         * @return A vector containing the results of the Haar coefficient calculation.
         */
        void findHaarCoefficients();
};

#endif // HARRCOEFFICIENT_HPP