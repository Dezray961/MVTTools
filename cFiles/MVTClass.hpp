#ifndef FINDMVT_HPP
#define FINDMVT_HPP
#include <vector>

// MVTAnalysis class declaration
class MVTAnalysis
{
    public:
        // data members
        std::vector<double> rate;
        std::vector<double> time;
        std::vector<double> rateErr;
        int lenghtOfData;
        std::vector<double> VTSet;
        std::vector<int> kSet;

        // constructor
        MVTAnalysis
        (
            std::vector<double> rateArray,
            std::vector<double> timeArray,
            std::vector<double> rateErrArray,
            int lenghtOfData,
            std::vector<int> kSet,
            int kMax
        );
        /**
         * Constructs a MVTAnalysis object with the given data arrays and parameters.
         * @param rateArray Pointer to an array of rate values.
         * @param timeArray Pointer to an array of time values.
         * @param rateErrArray Pointer to an array of rate error values.
         * @param lenghtOfData The length of the data arrays.
         * @param kSet A vector of k values to use for the analysis.
         * @param kMax The maximum value of k to use for the analysis.
         * @return A MVTAnalysis object with the given data and parameters.
         */


    private:
        // private member functions
        double localAverageOverKBins(int index, int k);
        /**
         * Finds the local average over k samples for a given index.
         * @param index The index of the data point to calculate the local average for.
         * @param k The number of samples to use for the local average.
         * @return The local average over k samples for the given index.
         */


        std::vector<double> localAverageOverKBinsForDataSet(int k);
        /**
         * Finds the local average over k samples for the entire data set.
         * @param k The number of samples to use for the local average.
         * @return A vector of local averages over k samples for the entire data set.
         */


        std::vector<double> squaredDifference(int k);
        /**
         * Finds the square of the difference between local averages for a given k.
         * @param k The number of samples to use for the local average.
         * @return A vector of squared differences between local averages for the given k.
         */

         
        double findVT(int k);
        /**
         * Finds the VT from the square of the difference between local averages for a given k.
         * @param k The number of samples to use for the local average.
         * @return The VT value for the given k.
         */


        void findVTForAllK();
        /**
         * Finds the VT for all k in the k set.
         * @return A vector of VT values for all k in the k set.
         */
};

#endif // FINDMVT_HPP