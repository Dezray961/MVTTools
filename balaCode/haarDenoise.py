"""
Haar wavelet denoising function. Vectorised implementation of the Haar wavelet denoising algorithm
found in https://github.com/sumanbala2210-USRA/GBM_MVT_paper. Original code by Suman Bala, modified
by Derek Pinkett

"""

from numpy import empty,log,sqrt,abs,median, ndarray, pad, diff, full

def haarDenoise(
        data: ndarray,
        error: ndarray = [],
        thresholdFactor: float = 1.0,
        estimateNoise: bool = False,
        soft: bool = False
        ):
    """
    Haar denoiser using a vectorised implementation of the Haar wavelet denoising algorithm. This function takes in a 1D array of data and applies Haar wavelet denoising to reduce noise while preserving important features in the signal.

    Args:
        data (ndarray): 1D array of input data to be denoised.
        error (ndarray, optional): 1D array of error values. Defaults to [].
        thresholdFactor (float, optional): Threshold factor for denoising. Defaults to 1.0.
        estimateNoise (bool, optional): Whether to estimate noise. Defaults to False.
        soft (bool, optional): Whether to use soft thresholding. Defaults to False.

    Returns:
        ndarray: Denoised data array.
    """

    data = data.astype('float32')
    lengthOfData: int = len(data)
    useError: bool = True

    if error is None or error.size == 0:
        useError = False
    else:
        # compute the variance of the error
        varianceError: ndarray = error.astype('float32') ** 2

    # determine the maximum decomposition level using the bit length of the data length minus one
    maxDecompositionLevel: int = (lengthOfData - 1).bit_length()
    # bitshifting to find the next power of two greater than or equal to the length of the data
    paddedSignalLength: int = 1 if (lengthOfData <= 0) else 1 << maxDecompositionLevel


    # extend to get 2^n length using symmetric reflection padding
    if paddedSignalLength > lengthOfData:
        # Calculate how much total padding is needed
        totalPaddingNeeded = paddedSignalLength - lengthOfData
        
        # Split the padding between left and right sides
        leftPaddingLength = totalPaddingNeeded // 2
        rightPaddingLength = totalPaddingNeeded - leftPaddingLength
        
        # pad the signal, 'reflect' mode mirrors the edge data accurately
        paddedSignalData = pad(
            data, 
            (leftPaddingLength, rightPaddingLength), 
            mode='reflect'
        )
        
        if useError:
            paddedNoiseVariance = pad(
                varianceError, 
                (leftPaddingLength, rightPaddingLength), 
                mode='reflect'
            )

    # initialize the fall-back noise estimate to 1.0
    estimatedNoiseSTD: float = 1.

    if (estimateNoise):
        # calculate the first difference of the signal to estimate noise
        signalFirstDifference: ndarray = abs(diff(paddedSignalData))

        # estimate the noise standard deviation using the median absolute deviation (MAD) method
        if (useError):
            squaredVarianceSum: ndarray = (
                paddedNoiseVariance[1:] ** 2 + paddedNoiseVariance[: -1] ** 2
                )
            propagatedNoiseSTD: ndarray = sqrt(0.5 * squaredVarianceSum)

            weightedDifferences: ndarray = signalFirstDifference / propagatedNoiseSTD
            estimatedNoiseSTD = 1.05 * median(weightedDifferences)
        else:
            estimatedNoiseSTD = 1.05 * median(signalFirstDifference)

    # compute and remove the baseline offset of the signal to center it around zero
    signalMeanOffset: float = paddedSignalData.mean()
    paddedSignalData  -=  signalMeanOffset

    # initialize arrays for the approximation coefficients and the reconstructed signal
    levelApproximationCoefficients: ndarray = empty(
        paddedSignalLength,
        dtype = 'float32')
    reconstructedSignal: ndarray = full(
        paddedSignalLength,
        signalMeanOffset,
        dtype = 'float32')
    
    # compute the cumulative sum of the padded signal data for the Haar transform
    paddedSignalData.cumsum(out = paddedSignalData)

    if (useError):
        # calculate the variance cumulative sum
        paddedNoiseVariance.cumsum(out = paddedNoiseVariance)
        levelVarianceCoefficients: ndarray = empty(
            paddedSignalLength,
            dtype = 'float32')

    # calculate the squared noise variance threshold for denoising
    squaredNoiseVarianceThreshold: float = (
        (2.0 * log(2.0)) 
        * ((thresholdFactor * estimatedNoiseSTD)
        * (thresholdFactor * estimatedNoiseSTD))
        )
    for decompositionLevelIndex in range(maxDecompositionLevel):
        # calculate scale using bitshifting
        waveletScaleBlockWidth: int = 1 << (maxDecompositionLevel - decompositionLevelIndex - 1)
        twoScale: int = 2 * waveletScaleBlockWidth

        # precalculate the inverse squared scale
        inverseSquaredScale: float = 1.0 / (twoScale * twoScale)

        # Forward tranform
        levelApproximationCoefficients[:-twoScale] = (
            2 * paddedSignalData[waveletScaleBlockWidth:-waveletScaleBlockWidth]
            - paddedSignalData[:-twoScale]
            - paddedSignalData[twoScale:]
        )
        levelApproximationCoefficients[-twoScale:-waveletScaleBlockWidth] = (
            2 * paddedSignalData[-waveletScaleBlockWidth:]
            - paddedSignalData[-twoScale:-waveletScaleBlockWidth]
            - paddedSignalData[:waveletScaleBlockWidth]
            - paddedSignalData[-1]
        )
        levelApproximationCoefficients[-waveletScaleBlockWidth:] = (
            2 * paddedSignalData[:waveletScaleBlockWidth]
            - paddedSignalData[-waveletScaleBlockWidth:]
            - paddedSignalData[waveletScaleBlockWidth:twoScale]
            + paddedSignalData[-1]
        )

        # level-dependent variance calculation
        if (useError):
            levelVarianceCoefficients[:-twoScale] = (
                varianceError[twoScale:] - varianceError[:-twoScale]
            )
            levelVarianceCoefficients[-twoScale:] = (
                varianceError[:twoScale] + varianceError[-1] - varianceError[-twoScale:]
            )
        else:
            # avoid allocating a whole array if the variance is uniform scalar
            levelVarianceCoefficients: int = twoScale

        # denoising thresholds
        thresholdLimit: float = (
            squaredNoiseVarianceThreshold * decompositionLevelIndex * levelVarianceCoefficients
        )
        noiseMask: ndarray = (
            (levelApproximationCoefficients * levelApproximationCoefficients) < thresholdLimit
        )

        # hard thresholding: set coefficients below the threshold to zero
        levelApproximationCoefficients[noiseMask] = 0

        if (soft):
            signalMask: ndarray = ~noiseMask
            survivingCoefficients: ndarray = levelApproximationCoefficients[signalMask]
            # prevent division by zero by ensuring non-zero coefficients
            squaredCoefficients: ndarray = (
                survivingCoefficients * survivingCoefficients
                + 1e-12
            )
            if (useError):
                varienceTerm: ndarray = (
                    thresholdLimit[signalMask] / squaredCoefficients
                )
            else:
                varienceTerm: ndarray = (
                    thresholdLimit / squaredCoefficients
                )
            # apply soft thresholding by scaling the surviving coefficients
            levelApproximationCoefficients[signalMask] *= (
                sqrt(1.0 - varienceTerm)
            )
        
        # set all noise coefficients to zero
        levelApproximationCoefficients[noiseMask] = 0.0


        # reconstruction via inverse transformation
        levelApproximationCoefficients.cumsum(out = levelApproximationCoefficients)

        # pre-calculate internal difference blocks to keep array operations vectorized
        diffMain: ndarray = (
            2 * levelApproximationCoefficients[waveletScaleBlockWidth:-waveletScaleBlockWidth]
            - levelApproximationCoefficients[:-twoScale]
            - levelApproximationCoefficients[twoScale:]
        )
        diffLeft: ndarray = (
            2 * levelApproximationCoefficients[-waveletScaleBlockWidth:]
            - levelApproximationCoefficients[-twoScale:-waveletScaleBlockWidth]
            - levelApproximationCoefficients[:waveletScaleBlockWidth]
            - levelApproximationCoefficients[-1]
        )
        diffRight: ndarray = (
            2 * levelApproximationCoefficients[:waveletScaleBlockWidth]
            - levelApproximationCoefficients[-waveletScaleBlockWidth:]
            - levelApproximationCoefficients[waveletScaleBlockWidth:twoScale]
            + levelApproximationCoefficients[-1]
        )

        # accumulate the results into the reconstructed signal
        reconstructedSignal[twoScale:] += (
            diffMain[::-1] * inverseSquaredScale
        )
        reconstructedSignal[waveletScaleBlockWidth:twoScale] += (
            diffLeft[::-1] * inverseSquaredScale
        )
        reconstructedSignal[:waveletScaleBlockWidth] += (
            diffRight[::-1] * inverseSquaredScale
        )

    # return the reconstructed signal, removing any padding that was added
    return reconstructedSignal[leftPaddingLength:leftPaddingLength+lengthOfData]
