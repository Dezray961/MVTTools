from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT
from analysisTools.haarDenoise import haarDenoise
from numpy import ceil, log, log10, ndarray, arange, zeros, newaxis, sqrt, zeros_like as zerosLike, clip
from numpy import tile, where, intersect1d, sum as sumElements, asarray, bincount, size, array, append
from numpy import exp, isnan, isinf
from scipy.stats import chi2
from scipy.ndimage import uniform_filter1d as uniformFilter1D
import pandas as pd
from loadConfig import getInitialBinSize
from warnings import filterwarnings
import matplotlib.pyplot as plt
from tqdm import tqdm
from os import path, makedirs

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

        # adaptively rebin the power spectrum by signal-to-noise ratio
        self.__adaptiveRebinBySNR()

        # evaluate the significance of the power spectrum and fit the results
        self.__evaluteSignificanceAndFit()


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
            self.__rate
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
                self.__deltaRate ** 2,
                size = halfStride
            )
            rawCoefficientErrorMatrix[levelIndex, :] = sqrt(squaredErrorMeans * 2.0) # propagate the error through the Haar wavelet transform
            rawCoefficientErrorMatrix[levelIndex, :] = clip(rawCoefficientErrorMatrix[levelIndex, :], 1.0e-10, None) # ensure non-zero error

            # propagate adjusted source error
            # calculate empirical variance of the signal to find intrinsic source variability
            localMeanOfSignal: ndarray = uniformFilter1D(
                self.__rate,
                size = scaleWidthInBins
            )
            localMeanOfSquared: ndarray = uniformFilter1D(
                self.__rate ** 2,
                size = scaleWidthInBins
            )
            localSignalVariance: ndarray = localMeanOfSquared - localMeanOfSignal ** 2
            localSignalVariance = clip(localSignalVariance, 0.0, None) # ensure non-negative variance

            # add localised background noise penalty to the instrumental foundation
            adjustedVariance: ndarray = (
                (rawCoefficientErrorMatrix[levelIndex, :] ** 2)
                + (localSignalVariance / scaleWidthInBins))
            adjustedCoefficientsMatrix[levelIndex, :] = sqrt(adjustedVariance) # propagate the error through the Haar wavelet transform
            adjustedCoefficientsMatrix[levelIndex, :] = clip(adjustedCoefficientsMatrix[levelIndex, :], 1.0e-10, None) # ensure non-zero error

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
            'rate': self.__rate,
            'totcounts': self.__totcounts
        })

        # denoise the signal timeline using the MODWT
        self.__rate = haarDenoise(data) # NOTE Golkhu's code also hands in the error to this step

        print("Padding the signal timeline to the next power of 2 for the Haar wavelet transform...")

        # calculate the span of a single repetition of the signal timeline
        maximumTimeDifference: float = self.__time.max() - self.__time.min()
        totalCycles: int = self.__numberOfRepetitions + 1
        
        # tile the measurement arrays
        self.__timeBinDuration = tile(self.__timeBinDuration, totalCycles)
        self.__rate = tile(self.__rate, totalCycles)
        self.__deltaRate = tile(self.__deltaRate, totalCycles)

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

    """
    this isn't right... This bins the data into the geometric bins by fixed SNR. There should be thousands
    new bins not 47...
    """

    @staticmethod
    def __rateRebin(
            timeBinStart: ndarray,
            timeBinEnd: ndarray,
            timeBinDuration: ndarray,
            netVariabilitySignal: ndarray,
            error: ndarray,
            targetSNR: float,
            minSignalClamp: float = 0.0,
            minErrorClamp: float = 2.0,
            fallbackValue: float = 0.0
    ) -> tuple[ndarray, ndarray, ndarray, ndarray, ndarray, ndarray]:
        """
        Rebins the input data based on the target signal-to-noise ratio (SNR).
        Args:
            timeBinStart (ndarray): Start times of the original bins.
            timeBinEnd (ndarray): End times of the original bins.
            timeBinDuration (ndarray): Durations of the original bins.
            netVariabilitySignal (ndarray): Net variability signal for each bin.
            error (ndarray): Error associated with each bin.
            targetSNR (float): Target signal-to-noise ratio for rebinning.
            minSignalClamp (float, optional): Minimum signal value to consider. Defaults to 0.0.
            minErrorClamp (float, optional): Minimum error value to consider. Defaults to 2.0.
            fallbackValue (float, optional): Fallback value for bins that do not meet criteria. Defaults to 0.0.

        Returns:
            tuple: Rebinned time bin start, end, duration, net variability signal, error, and mapping index.
        """
        # ensure that the input arrays are 1D numpy arrays
        startArray: ndarray = asarray(timeBinStart, dtype='float32').ravel()
        endArray: ndarray = asarray(timeBinEnd, dtype='float32').ravel()
        durationArray: ndarray = asarray(timeBinDuration, dtype='float32').ravel()
        signalArray: ndarray = asarray(netVariabilitySignal, dtype='float32').ravel()
        errorArray: ndarray = asarray(error, dtype='float32').ravel()

        # initialise lists to hold the rebinned data
        newStarts: list = []
        newEnds: list = []
        newDurations: list = []
        newSignals: list = []
        newErrors: list = []
        
        # mapping index to track which original bins contribute to each new bin
        numberOfOriginalBins: int = len(startArray)
        mappingIndex: ndarray = zeros(numberOfOriginalBins, dtype='int32') 
        
        currentNewBinIndex: int = 0
        i: int = 0

        # main sequential squeezing loop to merge bins until the target SNR is achieved
        while i < numberOfOriginalBins:
            # initialise a new consolidated tracking block
            blockStart: float = startArray[i]
            blockEnd: float = endArray[i]
            blockDuration: float = durationArray[i]
            blockSignal: float = signalArray[i]
            blockErrorSquared: float = errorArray[i] ** 2 # error propagates quadratically in variance space

            mappingIndex[i] = currentNewBinIndex

            j: int = i + 1
            # keep merging bins until the target SNR is achieved or we run out of bins
            while j < numberOfOriginalBins:
                # calculate the current block's SNR
                currentBlockError: float = sqrt(blockErrorSquared)

                # protect against nan values
                if isnan(blockSignal) or isnan(currentBlockError):
                    activeSNR: float = fallbackValue

                # guard against division by zero
                if currentBlockError <= minSignalClamp:
                    activeSNR: float = fallbackValue
                else:
                    activeSNR: float = blockSignal / currentBlockError
                
                # check if the current block's SNR meets the target
                if activeSNR >= targetSNR:
                    break # target SNR achieved, exit the inner loop

                # merge the next bin into the current block
                blockEnd = endArray[j]
                blockDuration += durationArray[j]
                blockSignal += signalArray[j]
                blockErrorSquared += errorArray[j] ** 2

                mappingIndex[j] = currentNewBinIndex
                j += 1
            
            # append the successfully completed block to the new rebinned lists
            newStarts.append(blockStart)
            newEnds.append(blockEnd)
            newDurations.append(blockDuration)
            newSignals.append(blockSignal)
            newErrors.append(sqrt(blockErrorSquared))

            # advance the current bin index to the next unprocessed original bin
            currentNewBinIndex += 1
            i = j

        # return as a tuple of numpy arrays for the rebinned data
        return (
            asarray(newStarts, dtype='float32'),
            asarray(newEnds, dtype='float32'),
            asarray(newDurations, dtype='float32'),
            asarray(newSignals, dtype='float32'),
            asarray(newErrors, dtype='float32'),
            mappingIndex
        )


    @staticmethod
    def __doRebin(
            originalArray: ndarray,
            mappingIndex: ndarray
    ) -> ndarray:
        """
        Rebins the original array based on the provided mapping index.
        Args:
            originalArray (ndarray): The original array to be rebinned.
            mappingIndex (ndarray): The mapping index indicating which new bin each original bin belongs to.

        Returns:
            ndarray: The rebinned array.
        """
        # ensure inputs are numpy arrays
        originalArray = asarray(originalArray, dtype='float32').ravel()
        mappingIndex = asarray(mappingIndex, dtype='int32').ravel()

        # validate the lengths of the input arrays match
        if originalArray.shape != mappingIndex.shape:
            raise ValueError("Original array and mapping index must have the same length.")
        
        # numpy bincount to sum the original array values according to the mapping index
        # weights=originalArray ensures that the values are summed correctly for each new bin
        rebinnedArray: ndarray = bincount(mappingIndex, weights=originalArray)

        return rebinnedArray


    def __adaptiveRebinBySNR(
            self
    ) -> None:
        """
        Third step of processing the light curve data. Does the following:
            1. Calculates the error for each bin based on the raw weights, adjusted weights, and term counts.
            2. Adaptively merges bins to satisfy the target SNR threshold.
            3. Extracts the new, merged bin configurations.
            4. Calculates the normalised power spectral density and its associated error.
            5. Collapses the statistical information to match the new bin structures.
            6. Calculates the tracking means per term inside the new bins.
            7. Finds indices where the power spectrum exceeds the significance threshold.
            8. Counts how many bins have significant signals.
        """        
        # prevent division by zero by clipping the adjusted weights to a minimum of 1.0
        safeAdjustedWeights: ndarray = clip(self.__binAdjustedWeightSum, 1.0, None)
        error: ndarray = sqrt(
            self.__binRawWeightSum * 2.0 * self.__binTermCounts / safeAdjustedWeights)
        print("Rebinning the power spectrum by signal-to-noise ratio...")

        # adaptively merge bins to satisfy the target SNR threshold
        rateRebinned: tuple = self.__rateRebin(
            self.__timeBinStart,
            self.__timeBinEnd,
            self.__timeBinDuration,
            self.__binChi2Sum,
            error,
            self.__SNRThreshold
        )
        print("Rebinning complete.")

        # extract the new, merged bin configurations
        self.__rebinnedTimeBinStart: ndarray = rateRebinned[0]
        self.__rebinnedTimeBinEnd: ndarray = rateRebinned[1]
        rebinnedTimeBinDuration: ndarray = rateRebinned[2]
        self.__rebinnedSignalSum: ndarray = rateRebinned[3]
        rebinnedErrorSum: ndarray = rateRebinned[4]
        rebinnedMappingIndex: ndarray = rateRebinned[5]
        print("Rebinned time bin start")
        print(self.__rebinnedTimeBinStart)
        print("Rebinned time bin end")
        print(self.__rebinnedTimeBinEnd)
        print("Rebinned time bin duration")
        print(rebinnedTimeBinDuration)
        print("Rebinned signal sum")
        print(self.__rebinnedSignalSum)
        print("Rebinned error sum")
        print(rebinnedErrorSum)
        print("Rebinned mapping index")
        print(rebinnedMappingIndex)

        # convert to logarithmic space
        self.__rebinnedSignalSum = log(self.__rebinnedSignalSum)

        # calculate the noralised power spectral density and its associated error
        safeDurations: ndarray = clip(rebinnedTimeBinDuration, 1.0, None)
        self.__powerSpectralDensity: ndarray = self.__rebinnedSignalSum / safeDurations
        self.__powerSpectralDensityError: ndarray = rebinnedErrorSum / safeDurations

        print("Collapsing statistical information to match the new bin structures...")
        # colapse the statistical information to match the new bin structures
        rebinnedRawWeightSum: ndarray = self.__doRebin(
            self.__binRawWeightSum,
            rebinnedMappingIndex,
        )
        rebinnedTermCounts: ndarray = self.__doRebin(
            self.__binTermCounts,
            rebinnedMappingIndex,
        )
        rebinnedAdjustedWeightSum: ndarray = self.__doRebin(
            self.__binAdjustedWeightSum,
            rebinnedMappingIndex,
        )
        print("Collapsing complete.")

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


    def __doPlot1(
            self
    ) -> None:
        """
        Plots the Haar scalogram. Only called if the plot1 flag is set to True in the config file.
        """
        print("Plotting Haar scalogram...")
        # initialize log-log figure axis wrapper
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xscale("log")
        ax.set_yscale("log")

        # plot the master data track
        ax.errorbar(
            self.__binCentreTimes,
            self.__allBinPSDSignals, 
            xerr=self.__binHalfWidths,
            yerr=self.__adjustedPSDErrors,
            fmt='.k',
            ecolor='gray',
            elinewidth=1,
            capsize=2
        )

        # explicitly mark the confirmed significant flux variations
        if (self.__numberOfSignificantBins > 0):
            ax.plot(
                self.__binCentreTimes[self.__significantSignalIndices], 
                self.__allBinPSDSignals[self.__significantSignalIndices], 
                'bv',
                markersize=8
            )

        # generate reference background power-law slopes
        xAxisLimits: ndarray = array([1.e-9, 1.e9])
        logMinDt: int = int(log10(min(self.__binCentreTimes)) * 2.0 - 4.0)
        logMaxDt: int = int(log10(max(self.__binCentreTimes)) * 2.0)

        for i in range(logMinDt, logMaxDt):
            slopeYValues: ndarray = self.__minimumYPlotLimit * xAxisLimits * exp(-i * log(10.0) / 2.0)
            ax.plot(xAxisLimits, slopeYValues, 'c:', alpha=0.5)

        # add the 1:1 reference line
        ax.plot(xAxisLimits, xAxisLimits, 'r-.', alpha=0.7, label='1:1 Scale Trend')

        # apply limits, titles, labels, and formatting
        xLowerLimit: float = min(self.__binCentreTimes) / 2.0
        xUpperLimit: float = max(self.__binCentreTimes) * 2.0
        ax.set_xlim((xLowerLimit, xUpperLimit))

        maximumYValue: float = max(append(self.__allBinPSDSignals[self.__significantSignalIndices], 1.0)) * 2.0
        ax.set_ylim((self.__minimumYPlotLimit, maximumYValue))

        ax.set_title('Haar Wavelet Flux Variability Profile', fontsize=14, fontweight='bold')
        ax.set_xlabel(r'$\Delta$T [seconds]', fontsize=12)
        ax.set_ylabel(r'Flux Variation $\sigma_{X,\Delta t}$ [%]', fontsize=12)
        ax.legend(loc='upper right')
        ax.grid(True, which="both", ls="--", alpha=0.3)

        # dynamic File Saving Architecture
        outputFileName: str = f"{self.data.name}FluxVariance.png"
        destinationDirectory: str = config.generalSettings.directories.processedDataPath + f"{self.data.name}/plots/"
        if not path.exists(destinationDirectory):
            makedirs(destinationDirectory)
        plt.savefig(f"{destinationDirectory}{outputFileName}", format='png', dpi=300)
        plt.show()
        plt.close(fig) # Memory efficient alternative to clf()


    def __evaluteSignificanceAndFit(
            self
    ) -> None:
        print(self.__numberOfSignificantBins)
        if self.__numberOfSignificantBins > 0:
            # isolate all bins that have valid error calculations
            validErrorIndices: ndarray = where(self.__powerSpectralDensityError > 0)
            self.__numberOfSignificantBins = len(validErrorIndices[0])

            # extract central plotting coordinates and geometric widths
            self.__binCentreTimes: ndarray = 0.5 * (
                self.__rebinnedTimeBinStart[validErrorIndices]
                + self.__rebinnedTimeBinEnd[validErrorIndices])
            self.__binHalfWidths: ndarray = 0.5 * (
                self.__rebinnedTimeBinEnd[validErrorIndices]
                - self.__rebinnedTimeBinStart[validErrorIndices])
            
            # pull master signal sets for the valid bins
            self.__allBinPSDSignals: ndarray = self.__powerSpectralDensity[validErrorIndices]
            basePSDErrors: ndarray = self.__powerSpectralDensityError[validErrorIndices]
            self.__allBinNoiseFloor: ndarray = self.__meanNoiseBaseline[validErrorIndices] # unused
            self.__allBinMeanWeights: ndarray = self.__meanAdjustedWeights[validErrorIndices]

            # apply safety-clipped statistical scaling to the error bars
            errorCorrectionFactor: ndarray = (1.0 + 2.0 * self.__rebinnedSignalSum[validErrorIndices])
            safeCorrectionFactor: ndarray = where(errorCorrectionFactor < 0.0, 0.0, errorCorrectionFactor)
            self.__adjustedPSDErrors: ndarray = basePSDErrors * sqrt(safeCorrectionFactor)

            # catogorise indices into insignificant noise vs significant signals
            insignificantBinIndices: ndarray = where(
                self.__allBinPSDSignals < self.__significanceSigmaThreshold * basePSDErrors)
            self.__numberOfInsignificantBins: int = size(insignificantBinIndices)
            self.__significantBinIndices: ndarray = where(
                self.__allBinPSDSignals >= self.__significanceSigmaThreshold * basePSDErrors)
            self.__numberOfSignificantBins: int = size(self.__significantBinIndices)

            if self.__numberOfInsignificantBins > 0:
                # calculate a threshold-shifted power limit for the insignificant points
                shiftedPowerLimits: ndarray = (
                    self.__allBinPSDSignals[insignificantBinIndices]
                    + self.__significanceSigmaThreshold * basePSDErrors[insignificantBinIndices])
                
                # clip negative values to zero and convert from power back to amplitude space
                safeShiftedPowerLimits: ndarray = where(shiftedPowerLimits < 0.0, 0.0, shiftedPowerLimits)
                self.__allBinPSDSignals[insignificantBinIndices] = sqrt(safeShiftedPowerLimits)
                self.__adjustedPSDErrors[insignificantBinIndices] = 0.0

            if self.__numberOfSignificantBins > 0:
                # isolate significant signal profiles
                self.__significantPSD: ndarray = self.__allBinPSDSignals[self.__significantBinIndices]
                self.__significantAdjustedPSD: ndarray = self.__adjustedPSDErrors[self.__significantBinIndices]
                self.__significantTimeScales: ndarray = self.__binCentreTimes[self.__significantBinIndices]
                self.__significantMeanWeights: ndarray = self.__allBinMeanWeights[self.__significantBinIndices]

                # transform unit space from power density to linear amplitude space
                self.__allBinPSDSignals[self.__significantBinIndices] = sqrt(self.__significantPSD)

                # propagate errors through the square root transform
                self.__adjustedPSDErrors[self.__significantBinIndices] = (
                    0.5 * self.__significantAdjustedPSD / sqrt(self.__significantPSD)) # BUG this is not correct because allBinPSDSignals has been modified on the previous line.
                
                # establish a visual lower boundary for plot scales
                self.__minimumYPlotLimit: float = min(self.__allBinPSDSignals[self.__significantBinIndices]) * 0.5

                if self.__plot1:
                    self.__doPlot1()
                








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
    data: LightCurveData = LightCurveData("GRB080319B", "Swift")
    chi2LookupTable: ndarray = getChi2LookupTable(data)
    analysis: HaarMTVFinder = HaarMTVFinder(data, chi2LookupTable, plot1 = True)


# BUG This isn't working, it doesn't find any significant signals. The logarithmic transform of the rate
# isn't working because the rate is sometimes negative. I have tried:
#   1. clamping the rate to a minimum, which is not ideal, it didn't work.
#   2. not using the logarithmic transform, which is not ideal, it didn't work.
#   3. adding a check to remove the nan/inf values before the rebinning, it didn't work.
# I might need to look at removing the background subtraction?
