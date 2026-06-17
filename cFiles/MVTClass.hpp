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
        int kMax;
        std::vector<int> kSet;
        std::vector<double> VTSet;
        double* VTSetArray;
        int kSetSize;
        std::vector<double> deltaT;
        std::vector<double> deltaTError;

        // constructor
        MVTAnalysis
        (
            double *rateArray,
            double *timeArray,
            double *rateErrArray,
            int lenghtOfData,
            bool generateEvenKSetFlag
        );
        // function to find the VT for all k in the k set
        std::vector<double> findVTForAllK
        (
            std::vector<std::vector<double>> dataSet
        );
        // function to conbvert the VT set to an array of doubles
        double* getVTSetAsArray();
    private:
        // private member functions
        // function to find kMax
        void findKMax();
        // function to find the k set
        void findKSet();
        // function to generate a kSet if the user wants all even numbers between 1 and kMax
        void generateEvenKSet();
        // function to find Δt for a given k
        void findDeltaTSet
        (
            std::vector<std::vector<double>> dataSet
        );
        // function to find the local average over k samples
        double localAverageOverKBins
        (
            std::vector<std::vector<double>> dataSet,
            int index,
            int k
        );
        // function to find the local average over k samples for the entire data set
        std::vector<double> localAverageOverKBinsForDataSet
        (
            std::vector<std::vector<double>> dataSet,
            int k
        );
        // function to find the square of the difference between local averages
        std::vector<double> squaredDifference
        (
            std::vector<std::vector<double>> dataSet,
            int k
        );
        // function to find the VT from the square of the difference between local averages
        double findVT
        (
            std::vector<std::vector<double>> dataSet,
            int k
        );
        // function to permute the data set
        std::vector<std::vector<double>> permuteDataSet
        (
            int index
        );
        // function to loop trhough N permutations of the data set and find the VT for each permutation
        void findVTForPermutations();
};

#endif // FINDMVT_HPP