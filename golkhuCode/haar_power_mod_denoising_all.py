
from numpy import *
import haar_nondec
import rate_rebin
import do_rebin
import haar_denoise            # Calling the haar_denoising modules
import matplotlib.pyplot as plt
import mu0_minimize_CHI2_fmin


def haar2_power_mod2a_Zach_denoising(inputFileName): 

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
                        *arange(numberOfBins) / (numberOfBins-1.))

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
            g00: ndarray[int] = where(powerSpectralDensityError > 0)
            numberOfSignificantBins: int = size(g00)
            all_t: ndarray[float] = 0.5*(rebinnedTimeBinStart[g00]+rebinnedTimeBinEnd[g00])
            all_dt: ndarray[float] = 0.5*(rebinnedTimeBinEnd[g00]-rebinnedTimeBinStart[g00]) 
            all_sig: ndarray[float] = powerSpectralDensity[g00]
            all_err0: ndarray[float] = powerSpectralDensityError[g00]
            all_err_chk: ndarray[float] = (1.+2*rebinnedSignalSum[g00]) 
            all_err: ndarray[float] = all_err0*sqrt(((all_err_chk<0).choose(all_err_chk,0))) #check for the argument being>0
            all_sig0: ndarray[float] = noiseBaseline[g00]
            sum2_2_2: ndarray[float] = meanAdjustedWeights[g00]
            validWeightIndices: ndarray[int] = where(all_sig < significanceSigmaThreshold*all_err0)
            numberOfInsignificantBins: int = size(validWeightIndices)
            g2: ndarray[int] = where(all_sig >=  significanceSigmaThreshold*all_err0)
            ng2: int = size(g2)
            if (numberOfInsignificantBins > 0): # check for the arguments >0
                all_sig_chk: ndarray[float] = all_sig[validWeightIndices]+significanceSigmaThreshold*all_err0[validWeightIndices]       
                all_sig[validWeightIndices] = sqrt(((all_sig_chk<0).choose(all_sig_chk,0)))
                all_err[validWeightIndices] = 0.
                
            if (ng2 > 0):
                # fitting mu_0 using chi2 minimization!
                pspec_p: ndarray[float] = all_sig[g2]
                dpspec_p: ndarray[float] = all_err[g2]
                tau_time: ndarray[float] = all_t[g2]
                sum22: ndarray[float] = sum2_2_2[g2]
                #
                #        
                all_sig[g2] = sqrt(all_sig[g2])
                #all_err[g2] = 0.5*all_err[g2]/sqrt(all_sig0[g2])
                all_err[g2] = 0.5*all_err[g2]/sqrt( all_sig[g2] )
                #
                miny: float = min(all_sig[g2])/2.
                ##subplot(111, xscale = "log", yscale = "log")
                #axis = [min_dta/2.,max_dta*2., miny,max(append(all_sig([g2],1.)))*2]
                #errorbar(all_t, all_sig, all_dt, all_err,'r.-')
                #
                """
                ax = plt.subplot(111)
                ax.set_xscale("log", nonposx = 'clip')
                ax.set_yscale("log", nonposy = 'clip')
                plt.errorbar(all_t, all_sig, xerr = all_dt, yerr = all_err, fmt = '.k')
                ax.set_xlim((min_dta/2.,max_dta*2.))
                maxy: float = max(append(all_sig[g2],1.))*2
                ax.set_ylim((miny, maxy))
                ax.set_title('Title')
                ax.set_xlabel(r'$\mathrm{\Delta T}$  [s]', fontsize = 12)
                ax.set_ylabel(r'Flux Variation  $\mathrm{\sigma_{X,\Delta t}}$  [%]', fontsize = 12)
                #ax.set_text(0.05, 0.9, 'Text goes here',
                #        fontsize = 14, transform = pl.gca().transAxes,
                #        ha = 'left', va = 'bottom')
                #plt.show()
                #       
                xx: ndarray[float] = array([1.e-9,1.e9])
                for i in xrange(int(log10(min_dt)*2.-4.), int(log10(max_dt)*2)):
                    plt.plot(xx, miny*xx*exp(-i*log(10.)/2.),'c:', markersize = 6)
                #
                if (ng > 0):
                    plt.plot(all_t[g1], all_sig[g1], 'bv')
                #
                plt.plot(xx,xx,'r-.')
                base1 = file_input+"_fluxVariance.png"
                base_address = '/home/zach/project_wavelet/Flux/'
                plt.savefig(base_address+base1,format = 'png')
                plt.clf()
                #plt.savefig('/Users/Vahid/Desktop/plots/file_input.png',format = 'png')
                #plt.savefig('/Users/Vahid/untitled/PLOTS/testplot.pdf',format = 'pdf')
                #Image.open('/Users/Vahid/untitled/PLOTS/testplot.png').save('/Users/Vahid/untitled/PLOTS/testplot.jpg','JPEG')
                #plt.show()
                """
                
                
                #   
                # fitting mu0:
                #import mu0_minimize_CHI2_fmin
                MU0, CHI2 = mu0_minimize_CHI2_fmin.mu0_minimize_CHI2_fmin(pspec_p,dpspec_p,tau_time,sum22)
                #
                #prnt_chi2 = file_input+'  : '+str(CHI2)
                #text_file2.write("%s,\n"%prnt_chi2)
                prnt_chi3: str = inputFileName
                text_file3.write("%s\n"%prnt_chi3)
                text_file2.write("%s,\n"%str(CHI2))
                prnt_tauMIN: str = inputFileName+':'+str(tau_time[0])
                text_file4.write("%s\n"%prnt_tauMIN)
                prnt_tauArray: str = inputFileName+':'+str(tau_time)
                text_file5.write("%s,\n"%prnt_tauArray)
                
                ####plt.figure(2)
                ###tau2_time = tau_time[1:]
                ###tau3_time = 0.5*(tau2_time[1:]+tau2_time[:-1])
                #prb = exp(-0.5*diff(CHI2))
                #plt.plot(tau3_time,prb,'bD-')
                ###chi2_diffTest = sqrt(diff(CHI2))
                ####plt.plot(tau3_time,chi2_diffTest,'r*--')
                ###base2 = file_input+"_CHI2.png"
                ###base_address = '/home/zach/project_wavelet/CHI2/'
                ###plt.savefig(base_address+base2,format = 'png')
                ###plt.clf()
                """
                whr  = where(diff(sign(mrg - prb))!= 0)
                sz_whr = size(whr)
                if (sz_whr > 0):
                    y1_1 = prb[whr[0]]
                    y1_2 = prb[whr[0]+1]
                    y2_1 = mrg[whr[0]]
                    y2_2 = mrg[whr[0]+1]
                    t_1 = tau3_time[whr[0]]
                    t_2 = tau3_time[whr[0]+1]
                    t_cross = ((t_2-t_1)*(y1_1-y2_1)-t_1*((y1_2-y1_1)-(y2_2-y2_1)))/((y2_2-y2_1)-(y1_2-y1_1))
                else:
                    t_cross = 'NA'
                """
                            
                wh_cv = where(CHI2 <=  chi2CriticalValues[:size(CHI2)])
                if (size(wh_cv) > 0):
                    ref = arange(size(wh_cv))
                    prnt_cvSize: str = inputFileName+':'+str(size(wh_cv))
                    text_file6.write("%s\n"%prnt_cvSize)
                    diff_ref = wh_cv - ref
                    wh_ref = where(diff_ref[0] !=  0)
                    if (size(wh_ref[0])  ==  0):
                        CHI2 = [CHI2[x] for x in wh_cv[0]]
                        #tau_time = [tau_time[y] for y in wh_cv[0]]
                        tau_time = tau_time[:max(wh_cv[0])+2]
                    else:
                        wh_indx = min(wh_ref[0])
                        CHI2 = CHI2[:wh_indx]
                        tau_time = tau_time[:wh_indx+1]
                        
                    CHI2 = array(CHI2)
                    tau_time = array(tau_time)
                    tau2_time = tau_time[1:]
                    tau3_time = 0.5*(tau2_time[1:]+tau2_time[:-1])
                    
                    chi2_diffTest = sqrt(abs(diff(CHI2)))
                    
                    if (chi2_diffTest.size > 0):
                        thrshld = 2.
                        sigma2_cl = where(chi2_diffTest >=  thrshld)
                        sigma2_cl2 = sigma2_cl[0]
                        if (sigma2_cl2.size  ==  0):
                            t_min = tau3_time[-1]
                            #prnt = file_input+'  : '+str(t_min)+'     D     '+str(chi2_diffTest[-1])
                            prnt: str = inputFileName+':'+str(t_min)+'  '+str(CHI2[-1]/(CHI2.size-1))+'  '+str(CHI2[-1])+'  '+str(CHI2.size-1)+' -1 '+str(chi2_diffTest[-1])
                            text_file1.write("%s\n"%prnt)
                        else:
                            sigma2_indx = sigma2_cl2[0]
                            if (sigma2_indx > 0):
                                y_1 = chi2_diffTest[sigma2_indx-1]
                                y_2 = chi2_diffTest[sigma2_indx]
                                t_1 = tau3_time[sigma2_indx-1]
                                t_2 = tau3_time[sigma2_indx]
                                t_min = (thrshld-y_2)/(y_2-y_1)*(t_2-t_1)+t_2
                                CHI2_0_dof = CHI2[sigma2_indx]/(sigma2_indx)
                                prnt: str = inputFileName+':'+str(t_min)+'  '+str(CHI2_0_dof)+'  '+str(CHI2[sigma2_indx])+'  '+str(sigma2_indx)+' 0 '+' 0 '
                                #prnt = file_input+'  : '+str(t_min)
                                text_file1.write("%s\n"%prnt)
                            else:
                                t_min = tau3_time[sigma2_indx]
                                prnt: str = inputFileName+':'+str(t_min)+'  '+str(CHI2[0])+'  '+' 0 '+'   ' +' 0 '+'  1  '+str(chi2_diffTest[0])
                                #prnt = file_input+'  : '+str(t_min)+'     U    '+str(chi2_diffTest[0])
                                text_file1.write("%s\n"%prnt)
                                
                        print(f't_min = {t_min}')
                    
    text_file1.close()    
    text_file2.close()
    text_file3.close()
    text_file4.close()
    text_file5.close()
    text_file6.close()

    return
