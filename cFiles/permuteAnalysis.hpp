#ifndef PERMUTEANALYSIS_HPP
#define PERMUTEANALYSIS_HPP

#include <vector>

// PermuteAnalysis class declaration
class PermuteAnalysis
{
    public:
        // data members
        std::vector<double> rate;
        std::vector<double> time;
        std::vector<double> rateErr;
        int lengthOfData;
        int kMax;
        std::vector<int> kSet;
        std::vector<double> VTSet;
        double* VTSetArray;
        int kSetSize;
        std::vector<double> deltaT;
        std::vector<double> deltaTError;


        PermuteAnalysis
        (
            double *rateArray,
            double *timeArray,
            double *rateErrArray,
            int lenghtOfData,
            int numberOfUniformKValues = 0
        );
        /**
         * Constructs a PermuteAnalysis object with the given data arrays and parameters.
         * @param rateArray Pointer to an array of rate values.
         * @param timeArray Pointer to an array of time values.
         * @param rateErrArray Pointer to an array of rate error values.
         * @param lenghtOfData The length of the data arrays.
         * @param numberOfUniformKValues The number of uniformly spaced k values to generate. If 0, a default set of k
         * values will be used.
         */

    private:
        // private member functions
        /**
         * There are a number of functions that are currently defined in the MVTAnalysis class that should be moved here.
         * They only need to be calculated once for the entire dataset so the results should be handed into the
         * MVTAnalysis class to avoid recalculating them for each permutation. These functions include:
         * findKMax()
         * findKSet()
         * generateEvenKSet - this should be changed to generateUniformKSet and should be able to generate a number of
         *                    uniformly spaced k values (log scale) based on a user defined number of k values to generate.
         * findDeltaTSet() - I need to look at the analysis to see if this is going to be the same for an average over 
         *                   all the permutations
         * */


        void findKMax();
        /**
         * Finds the maximum value of k based on the length of the data. Sets the kMax member variable to the largest
         * power of 2 that is less than or equal to the length of the data.
         */


        void findKSet();
        /**
         * Finds the set of k values based on the maximum value of k. Sets the kSet member variable to a vector 
         * containing the powers of 2 from 1 to kMax. 
         */


        void generateUniformKSet
        (
            int numberOfUniformKValues
        );
        /**
         * Generates a uniform set of k values based on a user-defined number of k values to generate. The k values are
         * generated on a log scale between 1 and kMax. Sets the kSet member variable to a vector containing the
         * uniformly spaced k values. This function replaces the generateEvenKSet function in the MVTAnalysis class, if
         * a non-zero value is passed to the generateUniformKSetFlag parameter in the constructor.
         */


        void findDeltaTSet();
        /**
         * Finds the set of delta t values based on the k set. This function calculates the mean and standard deviation
         * of the differences between time values for each k in the k set. Sets the deltaT and deltaTError member
         * variables to vectors containing the mean and standard deviation of the delta t values for each k in the k set.
         */


        std::vector<double> generatePermutation
        (
            int index
        );
        /**
         * Generates a permutation of the dataset based on the given index and runs the MVTAnalysis on it.
         * @param index The index of the permutation to generate.
         * @return A vector containing the results of the MVTAnalysis.
         */

         

        std::vector<std::vector<double>> runMVTAnalysisOnPermutations();
        /**
         * Manages multithreading and runs the MVTAnalysis on each permutation of the dataset.
         * @return A vector of vectors containding the average and standard deviation of the VT values for each value of
         * k in the kSet.
         */
};

#endif // PERMUTEANALYSIS_HPP