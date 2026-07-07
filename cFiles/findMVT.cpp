/*
g++ -fPIC -shared -std=c++17 -fopenmp findMVT.cpp MVTClass.hpp MVTClass.cpp -o findMVT.so
*/


# include "permuteAnalysis.hpp"

// extern "C" block to expose the PermuteAnalysis class to Python via ctypes
extern "C" {
    /**
     * Allocates a new PermuteAnalysis object with the given data arrays and parameters.
     * @param rate Pointer to an array of rate values.
     * @param time Pointer to an array of time values.
     * @param rateErr Pointer to an array of rate error values.
     * @param lenghtOfData The length of the data arrays.
     * @param numberOfTimeBins The number of time bins to use for the analysis.
     * @return A pointer to the newly allocated PermuteAnalysis object.
     */
    PermuteAnalysis* allocatePermuteAnalysis(
        double *rate,
        double *time,
        double *rateErr,
        int lenghtOfData,
        int numberOfTimeBins,
        int permutationStep
    )
    {
        return new PermuteAnalysis
        (
            rate,
            time,
            rateErr,
            lenghtOfData,
            numberOfTimeBins,
            permutationStep
        );
    }

    /**
     * Runs the permutation analysis on the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to run.
     */
    void runPermuteAnalysis(PermuteAnalysis* analysis)
    {
        analysis->runAnalysis();
    }

    /**
     * Gets the number of permutations completed for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @return The number of permutations completed.
     */
    int getPermutationsCompleted(PermuteAnalysis* analysis)
    {
        return analysis->getPermutationsCompleted();
    }

    /**
     * Gets the total number of permutations for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @return The total number of permutations.
     */
    int getTotalPermutations(PermuteAnalysis* analysis)
    {
        return analysis->getTotalPermutations();
    }

    /**
     * Checks if the analysis is complete for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @return 1 if the analysis is complete, 0 otherwise.
     */
    int isAnalysisComplete(PermuteAnalysis* analysis)
    {
        return analysis->isAnalysisComplete() ? 1 : 0;
    }

    /**
     * Frees the memory allocated for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to free.
     */
    void freePermuteAnalysis(PermuteAnalysis* analysis)
    {
        delete analysis;
    }

    /**
     * Gets a pointer to the log bin edges for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @return A pointer to the log bin edges array.
     */
    double* getLogBinEdges(PermuteAnalysis* analysis)
    {
        return analysis->logBinEdges.data();
    }

    /**
     * Gets a pointer to the log bin centers for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @return A pointer to the log bin centers array.
     */
    double* getLogBinCenters(PermuteAnalysis* analysis)
    {
        return analysis->logBinCenters.data();
    }

    /**
     * Gets a pointer to the power set averages for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @param index The index of the power set average to retrieve.
     * @return A pointer to the power set average at the given index.
     */
    double* getPowerSetAverages(PermuteAnalysis* analysis, int index)
    {
        return &analysis->powerSetAverages[index];
    }

    /**
     * Gets a pointer to the power set standard deviations for the given PermuteAnalysis object.
     * @param analysis Pointer to the PermuteAnalysis object to query.
     * @param index The index of the power set standard deviation to retrieve.
     * @return A pointer to the power set standard deviation at the given index.
     */
    double* getPowerSetStdDevs(PermuteAnalysis* analysis, int index)
    {
        return &analysis->powerSetStdDevs[index];
    }
}

int main()
{
    return 0;
}