
from numpy import *
import haar_nondec              # non-decimated Haar wavelet transform
import rate_rebin               # rebin to constant S/N
import do_rebin                 # rebin the statistical metadata weights to match the new bin structures
import haar_denoise             # Calling the haar_denoising modules
import matplotlib.pyplot as plt
import mu0MinimizeChi2Fmin      # chi^2 minimising function to find the best mu0 value for the power spectrum


def haar2_power_mod2a_Zach_denoising(
        inputFileName: str,
        plot1: bool = False
        ): 

    def shift(
            list: list[any],
            n: int
            ) -> list[any]:
        """ Shift the elements of a list by n positions to the left.

        Args:
            list list[any]: The list to be shifted.
            n int: The number of positions to shift.

        Returns:
            list[any]: The shifted list.
        """
        return list[n:] + list[:n]

    # open the output files for writing
    text_file1 = open("out2put_tMIN_10_23.txt", "a")
    text_file2 = open("out2put_chi2_10_23.txt", "a")
    text_file3 = open("out2put_name_10_23.txt", "a")
    text_file4 = open("out2put_tminMIN_10_23.txt", "a")
    text_file5 = open("out2put_tmin_array_10_23.txt", "a")
    text_file6 = open("out2put_cvSize_10_23.txt", "a")
    
    # define the minimum and maximum deltaT values
    minimumDeltaTime: float = 1.0e-4
    maximumDeltaTime: float = 1.0e3
    
    # define the number of bins
    binningFactor: float = 2.0
    numberOfBins: int = binningFactor * ceil(log(maximumDeltaTime /minimumDeltaTime ) / log(2.0))
    
    # define the deltaT and deltaT1 arrays
    logMinimumDeltaTime: float = log(minimumDeltaTime)/log(2.0)
    logMaximumDeltaTime: float = log(maximumDeltaTime)/log(2.0)
    
    # geometric spacing of deltaT bins
    timeBinStart: ndarray[float] = 2** (logMinimumDeltaTime +
                        (logMaximumDeltaTime - logMinimumDeltaTime)
                        *arange(numberOfBins) / (numberOfBins-1))

    numberOfBins = numberOfBins - 1 # subtract 1 because we are using deltaTime rather than bin edges


    timeBinEnd: ndarray[float] = shift(
        timeBinStart,1
        ) # shift the deltaTime array by 1 to create the deltaTime1 array
    
    # truncate the deltaTime and deltaTime1 arrays to the number of bins
    timeBinStart = timeBinStart[0:numberOfBins]
    timeBinEnd = timeBinEnd[0:numberOfBins]

    # initialize arrays
    binChi2Sum: ndarray[float] = zeros((numberOfBins),dtype = 'float32')
    binAdjustedWeightSum: ndarray[float] = zeros((numberOfBins),dtype = 'float32')
    binRawWeightSum: ndarray[float] = zeros((numberOfBins),dtype = 'float32')
    noiseBaseline: ndarray[float] = zeros((numberOfBins),dtype = 'float32')
    powerSpectrum: ndarray[float] = zeros((numberOfBins),dtype = 'float32') 
    powerSpectrumError: ndarray[float] = zeros((numberOfBins),dtype = 'float32')
    binTermCounts: ndarray[float] = zeros((numberOfBins),dtype = 'float32')

    # define the minimum and maximum deltaTimeArray values
    minimumDeltaTimeArray: float = timeBinStart.max()
    maximumDeltaTimeArray: float = 0.0

    # Load the lookup table for chi^2 critical values
    chi2ThresholdLookup: ndarray[float] = loadtxt('/Users/Vahid/python_codes/chi2_criticVal.txt',dtype = 'float32')
    degreesOfFreedom: ndarray[float] = chi2ThresholdLookup[:,0]
    chi2CriticalValues: ndarray[float] = chi2ThresholdLookup[:,1]
    
    # Load the data from the input file
    DATA: ndarray[float] = loadtxt(inputFileName, dtype = float32)

    # Check if the data is not empty
    if (DATA.size > 0):
        time: ndarray[float] = DATA[:,0]
        timeBinDurations: ndarray[float] = DATA[:,1]
        rate: ndarray[float] = DATA[:,2]
        deltaRate: ndarray[float] = DATA[:,3]
        logRate: ndarray[float] = log(rate)
        deltaLogRate: ndarray[float] = deltaRate/rate # propagate the error in log space
        
        # Denoising the data
        logRate: ndarray[float] = haar_denoise.haar_denoise(logRate,deltaLogRate)
        
        
        numberOfRepetitions: int = 1 
        """^^^^^^
        Technically, for a non-Haar wavelet transform, this would need to be L-1 where L is the length of the
        wavelet filter. For the Haar wavelet, L = 2, so L-1 = 1. Therefore, we can set numberOfRepetitions = 1.
        """

        maximumTimeDifference: float = time.max()-time.min()

        # make copies of the original time, deltaTime, logRate, and deltaLogRate arrays
        deltaTimeCopy: ndarray[float] = timeBinDurations
        timeCopy: ndarray[float] = time
        logRateCopy: ndarray[float] = logRate
        deltaLogRateCopy: ndarray[float] = deltaLogRate

        # 
        for k in range(0,numberOfRepetitions):
            time = concatenate(
                (time, timeCopy + (k+1) * maximumTimeDifference))
            timeBinDurations = concatenate(
                (timeBinDurations, deltaTimeCopy))
            logRate = concatenate(
                (logRate, logRateCopy))
            deltaLogRate = concatenate(
                (deltaLogRate, deltaLogRateCopy))
            

        
        # perform the Haar non-decimated wavelet transform on the logRate and deltaLogRate arrays
        waveletResults: tuple = haar_nondec.haar_nondec(
            time,
            logRate,
            deltaLogRate,
            16.0
            )
        """
        The imputs for this function do not match the inputs for the haar_nondec function in haar_nondec.py.
        The haar_nondec function in haar_nondec.py takes in a file path, number of bins, and two lists.
        The variable names below are best guesses based on the context of the code.
        """
        waveletArray: ndarray[float] = asarray(waveletResults)

        # extract physical quantities from the waveletArray
        timeIntervals: ndarray[float] = waveletArray[0,:]
        timeAverage: ndarray[float] = waveletArray[1,:]
        haarCoefficients: ndarray[float] = waveletArray[2,:]
        rawCoefficientError: ndarray[float] = waveletArray[3,:]
        adjustedCoefficientError: ndarray[float] = waveletArray[4,:]

        # calculate the power spectrum and its error
        waveletPower: ndarray[float] = haarCoefficients**2
        adjustedCoefficientVariance: ndarray[float] = adjustedCoefficientError**2*(numberOfRepetitions+1.0)/binningFactor
        rawCoefficientVariance: ndarray[float] = rawCoefficientError**2
        
        """
        Filter the Haar wavelet data into specific time-scale bins and aggregate their statistics.
        For each time-scale bin, sum the Chi-squared test statistic and store the weighting factors
        Track the minimum and maximum time-scale values encountered during the binning process.
        """
        for binIndex in range(0, int(numberOfBins)):
            validLowerBound: ndarray[int] = where(
                timeIntervals >=  timeBinStart[binIndex])
            validUpperBound: ndarray[int] = where(
                (timeIntervals < timeBinEnd[binIndex])
                * (adjustedCoefficientVariance > 0))
            matchingIndices: ndarray[int] = intersect1d(
                validLowerBound,
                validUpperBound)
            numberOfMatchingPoints: int = size(matchingIndices)

            if (numberOfMatchingPoints > 1):
                # aggregate the Chi-squared statistics for the current time-scale bin
                binChi2Sum[binIndex] = sum(
                    waveletPower[matchingIndices] / rawCoefficientVariance[matchingIndices])
                binAdjustedWeightSum[binIndex] = sum(
                    1.0 / adjustedCoefficientVariance[matchingIndices])
                binRawWeightSum[binIndex] = sum(
                    1.0 / rawCoefficientVariance[matchingIndices])
            
                binTermCounts[binIndex] = numberOfMatchingPoints

                if (timeBinStart[binIndex] < minimumDeltaTimeArray):
                    minimumDeltaTimeArray = timeBinStart[binIndex]
                if (timeBinEnd[binIndex] > maximumDeltaTimeArray):
                    maximumDeltaTimeArray = timeBinEnd[binIndex]
        
        
        
        # find bins that actually have data points
        validWeightIndices: ndarray[int] = where(binRawWeightSum > 0)
        numberOfValidWeightBins: int = size(validWeightIndices)
        if (numberOfValidWeightBins > 0):
            # calculate the weighted power spectrum and the expected noise baseline
            powerSpectrum[validWeightIndices] = (binChi2Sum[validWeightIndices]
                                        / binRawWeightSum[validWeightIndices])
            noiseBaseline[validWeightIndices] = (binTermCounts[validWeightIndices]
                                        / binRawWeightSum[validWeightIndices])
            
            # calculate the statistical uncertainty of the power spectrum
            weightedTermsRatio: ndarray[float] = (binAdjustedWeightSum[validWeightIndices]
                                                * binRawWeightSum[validWeightIndices]
                                                / binTermCounts[validWeightIndices])
            powerSpectrumError[validWeightIndices] = (sqrt(2.0) / sqrt(weightedTermsRatio))

        # identify bins where the power is statistically insignificant (less than 2 sigma)
        insigninficantBins: ndarray[int] = where(
            (abs(powerSpectrum) < 2.0 * powerSpectrumError) 
            & (powerSpectrumError > 0))
        numberOfInsignificantBins: int = size(insigninficantBins)




        # TODO: everything below this point needs to be implimented in findDeltat.py






        # setup variables for the rebinning process
        targetSignalToNoiseRatio: float = 3.0
        netVariabilitySignal: ndarray[float] = binChi2Sum  #-binTermCounts*1.0 

        # Prevent division by zero by clipping adjusted weights to a minimum of 1.0
        safeAdjustedWeights = clip(binAdjustedWeightSum, 1.0, None)
        error: ndarray[float] = sqrt(
            binRawWeightSum * 2.0 * binTermCounts / safeAdjustedWeights)
        timeBinDurations = binRawWeightSum
        
        # adaptively merge bins to satisfy the target Signal-to-Noise ratio
        rateRebin: tuple = rate_rebin.rate_rebin(
            timeBinStart,
            timeBinEnd,
            timeBinDurations,
            netVariabilitySignal,
            error,
            targetSignalToNoiseRatio,
            0.0,
            2.0,
            0.0
            )
        rebinMatrix: ndarray[float] = asarray(rateRebin)

        # Extract the new, merged bin configurations
        rebinnedTimeBinStart: ndarray[float] = rebinMatrix[0]
        rebinnedTimeBinEnd: ndarray[float] = rebinMatrix[1]
        rebinnedTimeBinDurations: ndarray[float] = rebinMatrix[2]
        rebinnedSignalSum: ndarray[float] = rebinMatrix[3]
        rebinnedErrorSum: ndarray[float] = rebinMatrix[4]
        rebinnedMappingIndex: ndarray[int] = rebinMatrix[5]

        # calculate the normalised power specral density and its associated error
        safeDurations = clip(rebinnedTimeBinDurations, 1.0, None)
        powerSpectralDensity: ndarray[float] = rebinnedSignalSum / safeDurations
        powerSpectralDensityError: ndarray[float] = rebinnedErrorSum / safeDurations
        
        # collapse the statistical metadata weights to match the new bin structures
        rebinnedRawWeightSum: ndarray[float] = do_rebin.do_rebin(
            binRawWeightSum,
            rebinnedMappingIndex
            )
        rebinnedTermCounts: ndarray[float] = do_rebin.do_rebin(
            binTermCounts,
            rebinnedMappingIndex
            )
        rebinnedAdjustedWeightSum: ndarray[float] = do_rebin.do_rebin(
            binAdjustedWeightSum,
            rebinnedMappingIndex
            )
        
        # final normalisation to calculate tracking means per term inside the new bins
        safeTermCounts: ndarray[float] = clip(rebinnedTermCounts, 1.0, None)
        safeRawWeights: ndarray[float] = clip(rebinnedRawWeightSum, 1.0, None)

        meanNoiseBaseline: ndarray[float] = rebinnedTermCounts / safeRawWeights
        meanVariabilitySignal: ndarray[float] = rebinnedSignalSum / safeTermCounts
        meanAdjustedWeights: ndarray[float] = rebinnedAdjustedWeightSum / safeTermCounts

        # rename for explicit statsitical context
        significanceSigmaThreshold: float = targetSignalToNoiseRatio

        # find indices where the power spectrum exceeds the 3-sigma confidence threshold
        significantSingnalIndices: ndarray[int] = where(
            (powerSpectralDensity > significanceSigmaThreshold * powerSpectralDensityError)
            & (powerSpectralDensityError > 0))
        
        # count how many bins have significant signals
        numberOfSignificantBins: int = size(significantSingnalIndices)











        if (numberOfSignificantBins>1):
            # isolate all bins that have valid error calculations
            validErrorIndices: ndarray[int] = where(powerSpectralDensityError > 0)
            numberOfSignificantBins: int = size(validErrorIndices)

            # extract central plotting coordinates and geometric widths
            binCenterTimes: ndarray[float] = 0.5 * (
                rebinnedTimeBinStart[validErrorIndices]
                + rebinnedTimeBinEnd[validErrorIndices])
            binHalfWidths: ndarray[float] = 0.5 * (
                rebinnedTimeBinEnd[validErrorIndices]
                - rebinnedTimeBinStart[validErrorIndices]) 
            
            # pull master signal sets for the valid bins
            allBinPSDSignals: ndarray[float] = powerSpectralDensity[validErrorIndices]
            basePSDErrors: ndarray[float] = powerSpectralDensityError[validErrorIndices]
            allBinNoiseFloor: ndarray[float] = noiseBaseline[validErrorIndices]
            allBinMeanWeights: ndarray[float] = meanAdjustedWeights[validErrorIndices]

            # apply safety-clipped statistical scaling to the error bars
            errorCorrectionFactor: ndarray[float] = (1.0 + 2 * rebinnedSignalSum[validErrorIndices]) 
            safeCorrectionFactor  = where(errorCorrectionFactor < 0, 0, errorCorrectionFactor)
            adjustedPSDErrors     = basePSDErrors * sqrt(safeCorrectionFactor)

            # catogorise indices into insignificant noise vs significant signals
            insigninficantBinIndices: ndarray[int] = where(allBinPSDSignals < significanceSigmaThreshold*basePSDErrors)
            numberOfInsignificantBins: int = size(insigninficantBinIndices)
            significantBinIndices: ndarray[int] = where(allBinPSDSignals >=  significanceSigmaThreshold*basePSDErrors)
            numberOfSignificantBins: int = size(significantBinIndices)






            if (numberOfInsignificantBins > 0): # check for the arguments >0
                # calculate a threshold-shifted power limit for the insignificant points
                shiftedPowerLimits: ndarray[float] = (
                    allBinPSDSignals[insigninficantBinIndices] 
                    + significanceSigmaThreshold
                    * basePSDErrors[insigninficantBinIndices])
                
                # clip negative values to zero and convert from power  back to amplitude space
                safePowerLimits = where(shiftedPowerLimits < 0, 0, shiftedPowerLimits)
                allBinPSDSignals[insigninficantBinIndices] = sqrt(safePowerLimits)
                adjustedPSDErrors[insigninficantBinIndices] = 0.0
                
            if (numberOfSignificantBins > 0):
                # isolate significant signal profiles
                significantPSD: ndarray[float] = allBinPSDSignals[significantBinIndices]
                significantAdjustedPSD: ndarray[float] = adjustedPSDErrors[significantBinIndices]
                significantTimeScales: ndarray[float] = binCenterTimes[significantBinIndices]
                significantMeanWeights: ndarray[float] = allBinMeanWeights[significantBinIndices]

                # cache the raw power before conversion if needed for standard error propagation
                # rawSignificantPower = allBinPSDSignals[significantBinIndices].copy()

                # transform unit space from power density to linear amplitude space
                allBinPSDSignals[significantBinIndices] = sqrt(allBinPSDSignals[significantBinIndices])

                # propagate error through the square root transforma
                adjustedPSDErrors[significantBinIndices] = (
                    0.5 * adjustedPSDErrors[significantBinIndices]
                    / sqrt( allBinPSDSignals[significantBinIndices])) # BUG! this is not correct because the allBinPSDSignals has been modified in the previous line. 
                
                # etablish a visual lower boundary for plot scales
                minimumYPlotLimit: float = min(allBinPSDSignals[significantBinIndices]) /2.0


                if plot1:
                    # initialize log-log figure axis wrapper
                    fig, ax = plt.subplots(figsize=(8, 6))
                    ax.set_xscale("log")
                    ax.set_yscale("log")

                    # plot the master data track
                    ax.errorbar(
                        binCenterTimes,
                        allBinPSDSignals, 
                        xerr=binHalfWidths,
                        yerr=adjustedPSDErrors,
                        fmt='.k',
                        ecolor='gray',
                        elinewidth=1,
                        capsize=2
                    )

                    # explicitly mark the confirmed significant flux variations
                    if (numberOfSignificantBins > 0):
                        ax.plot(
                            binCenterTimes[significantBinIndices], 
                            allBinPSDSignals[significantBinIndices], 
                            'bv',
                            markersize=8
                        )

                    # generate reference background power-law slopes
                    xAxisLimits = array([1.e-9, 1.e9])
                    logMinDt = int(log10(min(binCenterTimes)) * 2.0 - 4.0)
                    logMaxDt = int(log10(max(binCenterTimes)) * 2.0)

                    for i in range(logMinDt, logMaxDt):
                        slopeYValues = minimumYPlotLimit * xAxisLimits * exp(-i * log(10.0) / 2.0)
                        ax.plot(xAxisLimits, slopeYValues, 'c:', alpha=0.5)

                    # add the 1:1 reference line
                    ax.plot(xAxisLimits, xAxisLimits, 'r-.', alpha=0.7, label='1:1 Scale Trend')

                    # apply limits, titles, labels, and formatting
                    xLowerLimit = min(binCenterTimes) / 2.0
                    xUpperLimit = max(binCenterTimes) * 2.0
                    ax.set_xlim((xLowerLimit, xUpperLimit))

                    maximumYValue = max(append(allBinPSDSignals[significantBinIndices], 1.0)) * 2.0
                    ax.set_ylim((minimumYPlotLimit, maximumYValue))

                    ax.set_title('Haar Wavelet Flux Variability Profile', fontsize=14, fontweight='bold')
                    ax.set_xlabel(r'$\Delta$T [seconds]', fontsize=12)
                    ax.set_ylabel(r'Flux Variation $\sigma_{X,\Delta t}$ [%]', fontsize=12)
                    ax.legend(loc='upper right')
                    ax.grid(True, which="both", ls="--", alpha=0.3)

                    # dynamic File Saving Architecture
                    outputFileName = f"{inputFileName}_fluxVariance.png"
                    destinationDirectory = '/home/zach/project_wavelet/Flux/'
                    plt.savefig(f"{destinationDirectory}{outputFileName}", format='png', dpi=300)
                    plt.close(fig) # Memory efficient alternative to clf()

                
                # chi^2 minimization to find the best mu0 value for the power spectrum:
                optimalMu0, minimisedChi2 = mu0MinimizeChi2Fmin(
                    significantPSD,
                    significantAdjustedPSD,
                    significantTimeScales,
                    significantMeanWeights
                    )

                # save the results to the output files for further analysis and record-keeping
                # log the filename to the master dataset index registry
                text_file3.write(f"{inputFileName}\n")

                # append the minimized Chi-Square goodness-of-fit statistic
                text_file2.write(f"{minimisedChi2},\n")

                # record the minimum significant variability timescale (tau_min)
                # Index 0 represents the shortest timescale that survived the 3-sigma filter
                minimumSignificantTimeScale = significantTimeScales[0]
                text_file4.write(f"{inputFileName}:{minimumSignificantTimeScale}\n")

                # save the full list of active significant timescales for this signal
                text_file5.write(f"{inputFileName}:{list(significantTimeScales)},\n")

                # locate indices where the calculated Chi-Square falls into the noise floor
                maxComparisonLength: int = size(minimisedChi2)
                noiseFloorIndices = where(minimisedChi2 <=  chi2CriticalValues[:maxComparisonLength])
                if (size(noiseFloorIndices) > 0):
                    # log the size of the noise profile array to text_file6
                    text_file6.write(f"{inputFileName}:{size(noiseFloorIndices)}\n")

                    # generate sequential reference array to check for index continuity gaps
                    sequentialReference: ndarray[int] = arange(size(noiseFloorIndices))
                    continuityCheck = noiseFloorIndices - sequentialReference

                    # because continuityCheck is a flat 1D array, we unpack [0]
                    contiuityGaps = where(continuityCheck[0] !=  0)

                    # truncate using the flat array criteria
                    if (size(contiuityGaps[0])  ==  0):
                        # noise dominates from the stat; filter using the full noiseFloorIndices array
                        minimisedChi2 = array([minimisedChi2[index] for index in noiseFloorIndices])
                        maxNoiseIndex = int(max(noiseFloorIndices))
                        significantTimeScales = significantTimeScales[:maxNoiseIndex + 2]
                    else:
                        # signal was good initially; truncate at the first gap index
                        firstNoiseCutoffIndex = int(min(contiuityGaps))
                        minimisedChi2 = minimisedChi2[:firstNoiseCutoffIndex]
                        significantTimeScales = significantTimeScales[:firstNoiseCutoffIndex + 1]
                    
                    # ensure unifor numpy array formatting across truncated outputs
                    minimisedChi2 = asarray(minimisedChi2, dtype = 'float')
                    significantTimeScales = asarray(significantTimeScales, dtype = 'float')

                    # derive updated geometric coordinate spacings on the cropped data
                    shiftedTimeScales = significantTimeScales[1:]
                    doubleShiftedMidpoints = 0.5 * (shiftedTimeScales[1:] + shiftedTimeScales[:-1])
                    
                    # calculate absolute sigma step change on the cropped dataset
                    chi2SigmaDistanceChange = sqrt(abs(diff(minimisedChi2)))


                    if (chi2SigmaDistanceChange.size > 0):
                        sigmaThresholdLimit = 2.

                        # unpack the tuple to find wher the threshold condition is met
                        thresholdCrossings = where(chi2SigmaDistanceChange >=  sigmaThresholdLimit)[0]

                        if (thresholdCrossings.size  ==  0):
                            # scenario 1: no major breaks found; set to the last significant midpoint
                            minimumTrueTimeScale = doubleShiftedMidpoints[-1]
                            degreesOfFreedom = minimisedChi2.size - 1
                            reducedChi2 = minimisedChi2[-1] / degreesOfFreedom

                            logPayload: dict = {
                                "inputFileName": inputFileName,
                                "minimumTrueTimeScale": minimumTrueTimeScale,
                                "reducedChi2": reducedChi2,
                                "localReducedChi2": None,
                                "minimisedChi2": minimisedChi2[-1],
                                "degreesOfFreedom": degreesOfFreedom,
                                "chi2SigmaDistanceChange": chi2SigmaDistanceChange[-1],
                                "firstCrossingIndex": None,
                                "exitStatus": -1
                            }
                            text_file1.write(f"{logPayload}\n")
                        else:
                            # extract the target index wher the threshold was triggered
                            firstCrossingIndex = thresholdCrossings[0]

                            if (firstCrossingIndex > 0):
                                # scenario 2: perform linear interpolation to find the exact crossing point
                                yStart = chi2SigmaDistanceChange[firstCrossingIndex-1]
                                yEnd = chi2SigmaDistanceChange[firstCrossingIndex]

                                timeStart = doubleShiftedMidpoints[firstCrossingIndex-1]
                                timeEnd = doubleShiftedMidpoints[firstCrossingIndex]

                                # linear intersection calculation
                                minimumTrueTimeScale = (
                                    (sigmaThresholdLimit - yEnd)
                                    / (yEnd - yStart)
                                    * (timeEnd - timeStart)
                                    + timeEnd)
                                
                                localReducedChi2 = minimisedChi2[firstCrossingIndex] / (firstCrossingIndex)
                                
                                logPayload: dict = {
                                    "inputFileName": inputFileName,
                                    "minimumTrueTimeScale": minimumTrueTimeScale,
                                    "reducedChi2": None,
                                    "localReducedChi2": localReducedChi2,
                                    "minimisedChi2": minimisedChi2[firstCrossingIndex],
                                    "degreesOfFreedom": None,
                                    "chi2SigmaDistanceChange": chi2SigmaDistanceChange[firstCrossingIndex],
                                    "firstCrossingIndex": firstCrossingIndex,
                                    "exitStatus": 0
                                }
                                text_file1.write(f"{logPayload}\n")
                            else:
                                # scenario 3: threshold hit immediately
                                minimumTrueTimeScale = doubleShiftedMidpoints[firstCrossingIndex]
                                
                                logPayload: dict = {
                                    "inputFileName": inputFileName,
                                    "minimumTrueTimeScale": minimumTrueTimeScale,
                                    "reducedChi2": None,
                                    "localReducedChi2": None,
                                    "minimisedChi2": minimisedChi2[0],
                                    "degreesOfFreedom": None,
                                    "chi2SigmaDistanceChange": chi2SigmaDistanceChange[0],
                                    "firstCrossingIndex": None,
                                    "exitStatus": 1
                                }
                                text_file1.write(f"{logPayload}\n")
                                
                        print(f't_min = {minimumTrueTimeScale}')
                    
    text_file1.close()    
    text_file2.close()
    text_file3.close()
    text_file4.close()
    text_file5.close()
    text_file6.close()

    return
