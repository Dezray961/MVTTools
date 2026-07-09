from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT
from analysisTools.haarDenoise import haarDenoise
from numpy import ceil, log, ndarray, arange, zeros, newaxis, sqrt, zeros_like as zerosLike, clip
from numpy import tile, where, intersect1d, sum as sumElements, asarray
from scipy.stats import chi2
from scipy.ndimage import uniform_filter1d as uniformFilter1D
import pandas as pd
from loadConfig import getInitialBinSize
from warnings import filterwarnings
import matplotlib.pyplot as plt
from tqdm import tqdm

# suppress warnings for divide by zero and invalid value encountered in the logarithm and square root calculations
filterwarnings('ignore', r'divide by zero encountered')
filterwarnings('ignore', r'invalid value encountered')

class HaarMTVFinder:
    """
    A class to find the Haar MTV (Maximum Time Variation) in a light curve.
    Attributes:
        data (LightCurveData): The light curve data for which to find the Haar MTV.


    """
    def __init__(
            self,
            data: LightCurveData,
            chi2CriticalTablePointer: ndarray,
            plot1: bool = False,
    ) -> None:
        """docstring here"""

        # initialise class attributes from constructor arguments
        self.data: LightCurveData = data # The light curve data for which to find the Haar MTV.
        self.__plot1: bool = plot1 # Whether to plot the Haar MTV finding process.
        self.__chi2CriticalTable: ndarray = chi2CriticalTablePointer # The chi-squared critical values for the Haar MTV finding process. 



        # TODO logging...
        self.__loglevel: str = config.generalSettings.logging.logLevel
        if self.__loglevel != "NONE": # Whether to create a log file for the Haar MTV finding process.
            self.__logFile: str = f"{config.generalSettings.directories.logPath}/{self.data.name}/haarMTVFinderLog.txt"
        self.__silent: bool = config.generalSettings.logging.silentMode # Whether to print verbose output to the console.




        # Set the geometric bins for the Haar MTV finding process.
        self.__setGeometricBins()

        # initialise arrays
        self.__initialiseArrays()

        # initialise scalars
        self.__initialiseValues()

        # pad the signal timeline to prepare for the Haar wavelet transform
        self.__padSignalTimeline()

        # process the Haar wavelet transform
        self.__processHaarWaveletTransform()



    def __setGeometricBins(
            self
            ) -> None:
        """Sets the geometric bins for the Haar MTV finding process."""
        print("Setting geometric bins for Haar MTV finding process...")
        # define the minimum and maximum delta time values
        minimumDeltaTime: float = 1.0e-4
        maximumDeltaTime: float = 1.0e3

        # define the number of bins and the bin factor
        self.__binningFactor: float = 2.0
        self.__numberOfBins: int = int(
            self.__binningFactor * ceil(log(maximumDeltaTime / minimumDeltaTime) / log(2.0)))

        # log the minimum and maximum delta time values
        logMinimumDeltaTime: float = log(minimumDeltaTime) / log(2.0)
        logMaximumDeltaTime: float = log(maximumDeltaTime) / log(2.0)

        # geometrically spaced bins
        self.__timeBinStart: ndarray = 2.0 ** (
            logMinimumDeltaTime + (logMaximumDeltaTime - logMinimumDeltaTime)
            * arange(self.__numberOfBins) / (self.__numberOfBins - 1)
        )
        self.__timeBinEnd: ndarray = self.__shift(self.__timeBinStart, 1)

        # subtract 1 from the number of bins because we are using deltaTime rather than bin edges
        self.__numberOfBins -= 1

        # truncate the timeBinStart/End arrays to the number of bins
        self.__timeBinStart = self.__timeBinStart[:self.__numberOfBins]
        self.__timeBinEnd = self.__timeBinEnd[:self.__numberOfBins]
        print("Geometric bins set.")


    def __initialiseArrays(
            self
    ) -> None:
        """Initialises the arrays for the Haar MTV finding process. Seperate method to keep the constructor clean. Using float32 to save memory and improve performance."""
        self.__binChi2Sum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # sum of chi-squared values for each bin.
        self.__binAdjustedWeightSum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # sum of adjusted weights for each bin.
        self.__binRawWeightSum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # sum of raw weights for each bin.
        self.__noiseBaseline: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # noise baseline for each bin.
        self.__powerSpectrum: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # power spectrum for each bin.
        self.__powerSpectrumError: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # error in the power spectrum for each bin
        self.__binTermCounts: ndarray = zeros(
            (self.__numberOfBins),
            dtype='float32') # number of terms in each bin.
        
        # extract the data from the LightCurveData object and store it in arrays
        self.__time: ndarray = self.data.burstData['time'].to_numpy(dtype='float32') # time values for the light curve data.
        self.__timeBinDuration: ndarray = self.data.burstData['timeInBin'].to_numpy(dtype='float32') # duration of each time bin for the light curve data.
        self.__rate: ndarray = self.data.burstData['rate'].to_numpy(dtype='float32') # rate values for the light curve data.
        self.__deltaRate: ndarray = self.data.burstData['error'].to_numpy(dtype='float32') # error in the rate values for the light curve data.
        self.__logRate: ndarray = log(self.__rate) # logarithm of the rate values for the light curve data.
        self.__deltaLogRate: ndarray = self.__deltaRate / self.__rate # error in the logarithm of the rate values for the light curve data.
        self.__totcounts: ndarray = self.data.burstData['totcounts'].to_numpy(dtype='int') # total counts for the light curve data.


    def __initialiseValues(
            self
    ) -> None:
        """Initialises the variables for the Haar MTV finding process. Seperate method to keep the constructor clean."""
        self.__minimumDeltaTimeArray: float = self.__timeBinStart.max() # initalised to the maximum value of the timeBinStart array, will be updated during the Haar MTV finding process.
        self.__maximumDeltaTimeArray: float = 0.0 # initalised to 0, will be updated during the Haar MTV finding process.
        self.__degreesOfFreedom: int = len(self.__chi2CriticalTable) # initalised to the length of the chi2CriticalTable, will be updated during the Haar MTV finding process.
        self.__numberOfRepetitions: int = 1 # techincally this is only 1 for the Haar transform, if DB2, or some other WT, is implimented then it will need to be changed.
        self.__initialBinSize: float = float(getInitialBinSize(config, self.data.source)) # get the initial bin size from the configuration file based on the source of the data.
        self.__SNRThreshold: float = float(config.mvtAnalysisConfig.signalToNoiseThreshold) # get the signal to noise threshold from the configuration file. 
        self.__significanceSigmaThreshold: float = float(config.mvtAnalysisConfig.significanceSigmaThreshold) # get the significance sigma threshold from the configuration file.


    @staticmethod
    def __shift(
        list: list,
        n: int
    ) -> list:
        """ Shift the elements of a list by n positions to the left.

        Args:
            list list[any]: The list to be shifted.
            n int: The number of positions to shift.

        Returns:
            list[any]: The shifted list.
        """
        return list[n:] + list[:n]


    def __MODWTWrapper(
            self
    ) -> tuple[ndarray, ndarray, ndarray, ndarray, ndarray]:
        """
        Wrapper for the MODWT class to perform the Haar wavelet transform on the light curve data. Does the following:
            1. Performs the MODWT on the light curve data.
            2. Calculates the geometric variables for the Haar wavelet coefficients.
        
        Returns:
            tuple: A tuple containing the Haar wavelet coefficients, scaling coefficients, time intervals, raw coefficient errors, and adjusted coefficients.

        """
        print("Performing Haar wavelet transform...")
        # perform the Haar wavelet transform on the light curve data
        waveletCoefficients, scalingCoefficients = MODWT(
            self.__logRate
        )
        print("Haar wavelet transform complete.")
        print("Calculating geometric variables for Haar wavelet coefficients...")
        numberOfLevels, totalDataPoints = waveletCoefficients.shape

        # match and reconstruct variable spaces
        haarCoefficientMatrix: ndarray = waveletCoefficients * sqrt(2.0) # "un-normalize" the wavelet coefficients

        # allocate arrays for the time intervals, raw coefficients, and adjusted coefficients
        timeIntervalsMatrix: ndarray = zerosLike(haarCoefficientMatrix)
        rawCoefficientErrorMatrix: ndarray = zerosLike(haarCoefficientMatrix)
        adjustedCoefficientsMatrix: ndarray = zerosLike(haarCoefficientMatrix)

        # loop over the MODWT levels to dynamically calculate the geometric variables
        for levelIndex in range(numberOfLevels):
            # convert level to its actual geometric binwidth level = 2^(levelIndex + 1)
            scaleWidthInBins: int = 2 ** (levelIndex + 1)
            halfStride: float = int(scaleWidthInBins / 2.0)

            # track physicsal time duration for this specific level of the MODWT
            physicalDuration: float = scaleWidthInBins * self.__initialBinSize
            timeIntervalsMatrix[levelIndex, :] = physicalDuration

            # propagate the instrumental error
            # compute rolling squared mean over the active filter window size
            squaredErrorMeans: ndarray = uniformFilter1D(
                self.__deltaLogRate ** 2,
                size = halfStride
            )
            rawCoefficientErrorMatrix[levelIndex, :] = sqrt(squaredErrorMeans * 2.0) # propagate the error through the Haar wavelet transform

            # propagate adjusted source error
            # calculate empirical variance of the signal to find intrinsic source variability
            localMeanOfSignal: ndarray = uniformFilter1D(
                self.__logRate,
                size = scaleWidthInBins
            )
            localMeanOfSquared: ndarray = uniformFilter1D(
                self.__logRate ** 2,
                size = scaleWidthInBins
            )
            localSignalVariance: ndarray = localMeanOfSquared - localMeanOfSignal ** 2
            localSignalVariance = clip(localSignalVariance, 0.0, None) # ensure non-negative variance

            # add localised background noise penalty to the instrumental foundation
            adjustedVariance: ndarray = (
                (rawCoefficientErrorMatrix[levelIndex, :] ** 2)
                + (localSignalVariance / scaleWidthInBins))
            adjustedCoefficientsMatrix[levelIndex, :] = sqrt(adjustedVariance) # propagate the error through the Haar wavelet transform

        print("Geometric variables for Haar wavelet coefficients calculated.")

        return haarCoefficientMatrix, scalingCoefficients, timeIntervalsMatrix, rawCoefficientErrorMatrix, adjustedCoefficientsMatrix


    def __padSignalTimeline(
            self
    ) -> None:
        """First step of processing the light curve data. Does the following:
            1. Denoises the signal.
            2. Pads the signal to the next power of 2 in order to negate edge effects.
            3. Performs the Haar wavelet transform.
            4. Stores the Haar wavelet coefficients, scaling coefficients, time intervals, raw coefficient errors, and adjusted coefficients as class attributes.
        """
        # create a dataframe to hand into the haarDenoise function
        data: pd.DataFrame = pd.DataFrame({
            'rate': self.__logRate,
            'totcounts': self.__totcounts
        })

        # denoise the signal timeline using the MODWT
        self.__logRate = haarDenoise(data) # NOTE Golkhu's code also hands in the error to this step

        print("Padding the signal timeline to the next power of 2 for the Haar wavelet transform...")

        # calculate the span of a single repetition of the signal timeline
        maximumTimeDifference: float = self.__time.max() - self.__time.min()
        totalCycles: int = self.__numberOfRepetitions + 1
        
        # tile the measurement arrays
        self.__timeBinDuration = tile(self.__timeBinDuration, totalCycles)
        self.__logRate = tile(self.__logRate, totalCycles)
        self.__deltaLogRate = tile(self.__deltaLogRate, totalCycles)

        # construct an advancing time array for the tiled signal timeline
        timeOffset: float = maximumTimeDifference * arange(totalCycles)
        self.__time = (self.__time[:, newaxis] + timeOffset).ravel(order = 'F')

        print("Signal timeline padded.")
        print("Performing Haar wavelet transform on the padded signal timeline...")

        # perform the Haar wavelet transform on the tiled signal timeline
        MODWTResults: tuple = self.__MODWTWrapper()
        self.__haarCoefficientMatrix: ndarray = MODWTResults[0]
        self.__scalingCoefficients: ndarray = MODWTResults[1]
        self.__timeIntervalsMatrix: ndarray = MODWTResults[2]
        self.__rawCoefficientErrorMatrix: ndarray = MODWTResults[3]
        self.__adjustedCoefficientsMatrix: ndarray = MODWTResults[4]

        print("Haar wavelet transform on the padded signal timeline complete.")


    def __processHaarWaveletTransform(
            self
    ) -> None:
        """
        Second step of processing the light curve data. Does the following:
            1. Calculates the power spectrum and its error.
            2. Filters the Haar wavelet data into specific time-scale bins and aggregates their statistics
            3. Calculates the weighted power spectrum and the expected noise baseline.
        """
        # calculate the power spectrum and its error
        waveletPower: ndarray = self.__haarCoefficientMatrix ** 2
        adjustedCoefficientVariance: ndarray = (self.__adjustedCoefficientsMatrix ** 2
                                                * (self.__numberOfRepetitions + 1)
                                                / self.__binningFactor)
        rawCoefficientVariance: ndarray = self.__rawCoefficientErrorMatrix ** 2

        """
        Filter the Haar wavelet data into specific time-scale bins and aggregate their statistics.
        For each time-scale bin, sum the Chi-squared test statistic and store the weighting factors
        Track the minimum and maximum time-scale values encountered during the binning process.
        """
        for binIndex in tqdm(range(self.__numberOfBins), desc="Filtering bins", unit="bin"):
            # find the indices of the Haar wavelet coefficients that fall within the current time-scale bin
            validLowerBound: ndarray = where(
                self.__timeIntervalsMatrix >= self.__timeBinStart[binIndex])
            validUpperBound: ndarray = where(
                self.__timeIntervalsMatrix < self.__timeBinEnd[binIndex])
            matchingIndices: ndarray = intersect1d(validLowerBound, validUpperBound)
            numberOfMatchingIndices: int = len(matchingIndices)

            if numberOfMatchingIndices > 1:
                # unroll the matrices to 1D arrays to prevent memory issues when summing the Chi-squared test statistic and weighting factors
                waveletPower = waveletPower.ravel()[matchingIndices]
                rawCoefficientVariance = rawCoefficientVariance.ravel()[matchingIndices]
                adjustedCoefficientVariance = adjustedCoefficientVariance.ravel()[matchingIndices]

                # aggregate the Chi-squared test statistic for the current time-scale bin
                self.__binChi2Sum[binIndex] = sumElements(
                    waveletPower / rawCoefficientVariance)
                self.__binAdjustedWeightSum[binIndex] = sumElements(
                    1.0 / adjustedCoefficientVariance)
                self.__binRawWeightSum[binIndex] = sumElements(
                    1.0 / rawCoefficientVariance)
                self.__binTermCounts[binIndex] = numberOfMatchingIndices

                # update the minimum and maximum time-scale values encountered during the binning process
                if self.__timeBinStart[binIndex] < self.__minimumDeltaTimeArray:
                    self.__minimumDeltaTimeArray = self.__timeBinStart[binIndex]
                if self.__timeBinEnd[binIndex] > self.__maximumDeltaTimeArray:
                    self.__maximumDeltaTimeArray = self.__timeBinEnd[binIndex]

        # find bins that actually have data points
        validWeightIndices: ndarray = where(self.__binRawWeightSum > 0.0)
        numberOfValidWeightBins: int = len(validWeightIndices[0])
        if numberOfValidWeightBins > 0:
            # calculate the weighted power spectrum and the expected noise baseline
            self.__powerSpectrum[validWeightIndices] = (
                self.__binChi2Sum[validWeightIndices]
                / self.__binRawWeightSum[validWeightIndices])
            self.__noiseBaseline[validWeightIndices] = (
                self.__binTermCounts[validWeightIndices]
                / self.__binRawWeightSum[validWeightIndices])

            # calculate the statistical uncertainty in the power spectrum
            weightedTermsRatio: ndarray = (
                self.__binAdjustedWeightSum[validWeightIndices]
                * self.__binRawWeightSum[validWeightIndices]
                / self.__binTermCounts[validWeightIndices] ** 2)
            self.__powerSpectrumError[validWeightIndices] = (
                sqrt(2.0) / sqrt(weightedTermsRatio))
            
        # identify bins where the power is statistically insignificant (less than 2 sigma)
        insignificantBins: ndarray = where(
            (abs(self.__powerSpectrum) < 2.0 * self.__powerSpectrumError)
            & (self.__powerSpectrumError > 0))
        self.__numberOfInsignificantBins: int = len(insignificantBins[0])


    def __adaptiveRebinBySNR(
            self
    ) -> None:
        
        # prevent division by zero by clipping the adjusted weights to a minimum of 1.0
        safeAdjustedWeights: ndarray = clip(self.__binAdjustedWeightSum, 1.0, None)
        error: ndarray = sqrt(
            self.__binRawWeightSum * 2.0 * self.__binTermCounts / safeAdjustedWeights)
        
        # adaptively merge bins to satisfy the target SNR threshold
        rateRebinned: tuple = rateRebin(
            self.__timeBinStart,
            self.__timeBinEnd,
            self.__timeBinDuration,
            self.__binChi2Sum,
            error,
            self.__SNRThreshold
        )
        rebinMatrix: ndarray = asarray(rateRebinned)

        # extract the new, merged bin configurations
        self.__rebinnedTimeBinStart: ndarray = rebinMatrix[0]
        self.__rebinnedTimeBinEnd: ndarray = rebinMatrix[1]
        rebinnedTimeBinDuration: ndarray = rebinMatrix[2]
        self.__rebinnedSignalSum: ndarray = rebinMatrix[3]
        rebinnedErrorSum: ndarray = rebinMatrix[4]
        rebinnedMappingIndex: ndarray = rebinMatrix[5]

        # calculate the noralised power spectral density and its associated error
        safeDurations: ndarray = clip(rebinnedTimeBinDuration, 1.0, None)
        self.__powerSpectralDensity: ndarray = self.__rebinnedSignalSum / safeDurations
        self.__powerSpectralDensityError: ndarray = rebinnedErrorSum / safeDurations

        # colapse the statistical information to match the new bin structures
        rebinnedRawWeightSum: ndarray = doRebin(
            self.__binRawWeightSum,
            rebinnedMappingIndex,
        )
        rebinnedTermCounts: ndarray = doRebin(
            self.__binTermCounts,
            rebinnedMappingIndex,
        )
        rebinnedAdjustedWeightSum: ndarray = doRebin(
            self.__binAdjustedWeightSum,
            rebinnedMappingIndex,
        )

        # final normalisation to calculate tracking means per term inside the new bins
        safeTermCounts: ndarray = clip(rebinnedTermCounts, 1.0, None)
        safeRawWeights: ndarray = clip(rebinnedRawWeightSum, 1.0, None)

        self.__meanNoiseBaseline: ndarray = rebinnedTermCounts / safeRawWeights
        self.__meanVariabilitySignal: ndarray = self.__rebinnedSignalSum / safeTermCounts
        self.__meanAdjustedWeights: ndarray = rebinnedAdjustedWeightSum / safeTermCounts

        # find indices where the power spectrum exceeds the significance threshold
        self.__significantSignalIndices: ndarray = where(
            (self.__powerSpectralDensity > self.__significanceSigmaThreshold * self.__powerSpectralDensityError)
            & (self.__powerSpectralDensityError > 0)
        )

        # count how many bins have significant signals
        self.__numberOfSignificantBins: int = len(self.__significantSignalIndices[0])


# TODO implement the rateRebin and doRebin functions!




















    def __evaluteSignificanceAndFit(
            self
    ) -> None:
        ...


def getChi2LookupTable(
        data: LightCurveData,
    ) -> ndarray:
    """
    Generates a chi-squared lookup table for a given LightCurveData object. The lookup table is used to determine the critical chi-squared value for a given number of degrees of freedom.

    Args:
        data (LightCurveData): The light curve data for which to generate the lookup table.

    Returns:
        ndarray: The chi-squared lookup table.
    """
    expectedDegreesOfFreedom: int = len(data.burstData) + 1
    return chi2.ppf(0.95, arange(0, expectedDegreesOfFreedom))


if __name__ == "__main__":
    from loadConfig import importConfiguration
    config = importConfiguration()
    print(config.generalSettings.directories.dataPath)
    data: LightCurveData = LightCurveData("GRB080319B", "Swift")
    chi2LookupTable: ndarray = getChi2LookupTable(data)
    analysis: HaarMTVFinder = HaarMTVFinder(data, chi2LookupTable, plot1 = True)