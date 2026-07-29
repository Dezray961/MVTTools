from numpy import ndarray, log, mean, floor, log2, sign, cumsum, sqrt, maximum, where, asarray
from pandas import DataFrame
from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT, inverseMODWT
from tqdm import tqdm



def haarDenoise(
        data: DataFrame,
        thresholdMethod: str = "soft",
        backgroundNoiseModel: str = "poisson",
        thresholdScaleFactor: float = 0.5
        ) -> ndarray:
    """Denoises a light curve dataframe using the Haar wavelet transform and thresholding.

    Args:
        data (DataFrame): light curve dataframe to denoise. The dataframe must contain ``rate`` and ``totcounts`` columns.
        thresholdMethod (str, optional): thresholding method to use. Defaults to "soft". Options are "soft" and "hard". 
        backgroundNoiseModel (str, optional): background noise model to use. Defaults to "poisson". Options are "poisson" and "gaussian".
        thresholdScaleFactor (float, optional): factor to scale the threshold value. Defaults to 0.5. The rationale for this is that any noise that is left in should be small enough to not affect the MVT calculation. Whereas, if the threshold is too high, the MVT calculation will be affected by the loss of signal. This is a trade-off between noise and signal loss.

    Returns:
        ndarray: denoised rate values
    """
    def thresholdPoisson(
            level: int,
            maxLevel: int,
            meanRate: float
            ) -> float:
        """Calculates the threshold for the Haar denoising 

        Args:
            level (int): level of the wavelet transform
            maxLevel (int): maximum level of the wavelet transform
            meanRate (float): mean rate of the pre-burst data

        Returns:
            float: threshold for the given level
        """
        n_l: int = 2 ** (maxLevel - level)
        prefactor: float = 2 ** (-0.5 * (level + 2))
        firstLogTerm: float = log(n_l)
        safeMeanRate: float = max(meanRate, 1e-12)
        secondLogTerm: float = log(n_l * safeMeanRate)
        sqrtTerm: float = max((((4 * firstLogTerm) ** 2) + (8 * secondLogTerm)), 0) # max with 0 to avoid complex numbers
        thresholdValue: float = prefactor * (2 * firstLogTerm + sqrtTerm ** 0.5)
        return thresholdValue


    def thresholdGaussian(
            level: int,
            maxLevel: int
            ) -> float:
        """Calculates the threshold for the Haar denoising using the Gaussian noise model

        Args:
            level (int): level of the wavelet transform
            maxLevel (int): maximum level of the wavelet transform

        Returns:
            float: threshold for the given level
        """
        n_l: int = 2 ** (maxLevel - level)
        thresholdValue: float = sqrt(2 * log(n_l))
        return thresholdValue            


    def findPoissonRate(
            data: DataFrame,
            binSize: int
            ) -> float:
        """Finds the mean rate of the pre-burst data for a given bin size.

        Args:
            data (DataFrame): light curve data
            binSize (int): size of the bins in indices.

        Returns:
            float: mean rate of the pre-burst data
        """
        preBurstCounts: ndarray = data['totcounts'].to_numpy()
        timeInBin: float = 100e-6 * binSize
        # moving window sum
        cumulativeCounts: ndarray = cumsum(preBurstCounts)
        countsInBin: ndarray = cumulativeCounts[binSize:] - cumulativeCounts[:-binSize]
        meanRate: float = mean(countsInBin) / timeInBin
        return meanRate

    burstData: list[float] = data['rate'].to_list()

    maxLevel: int = int(floor(log2(len(burstData)))) - 1

    # get the Haar wavelet coefficients and scaling coefficients using the MODWT
    waveletCoeffs, scalingCoeffs = MODWT(burstData)

    # apply thresholding to the wavelet coefficients
    for level in range(maxLevel):
        # calculate the threshold for the current level
        meanRate: float = findPoissonRate(data, 2 ** level)
        
        if backgroundNoiseModel == "poisson":
            thresholdValue: float = thresholdPoisson(level, maxLevel, meanRate)
        elif backgroundNoiseModel == "gaussian":
            thresholdValue: float = thresholdGaussian(level, maxLevel)
        else:
            raise ValueError(f"Invalid background noise model: {backgroundNoiseModel}. Must be 'poisson' or 'gaussian'.")

        thresholdValue *= thresholdScaleFactor # scaling factor to adjust the threshold value

        # Convert coefficients to a NumPy array if they aren't already
        coeffs = asarray(waveletCoeffs[level])
        
        # apply the thresholding using vectorized NumPy operations
        if thresholdMethod == "soft":
            # Soft thresholding: sign(x) * max(0, |x| - threshold)
            abs_coeffs = abs(coeffs)
            waveletCoeffs[level] = sign(coeffs) * maximum(0, abs_coeffs - thresholdValue)
            
        elif thresholdMethod == "hard":
            # Hard thresholding: keep x if |x| >= threshold, else 0
            waveletCoeffs[level] = where(abs(coeffs) >= thresholdValue, coeffs, 0.0)
            
        else:
            raise ValueError(f"Invalid threshold method: {thresholdMethod}. Must be 'soft' or 'hard'.")

    # invert the transform to get the denoised data
    denoisedData: ndarray = inverseMODWT(
        waveletCoeffs,
        scalingCoeffs
        )
    return denoisedData


if __name__ == "__main__":
    def plotlightCurve(
        data: DataFrame
        ) -> None:
        import matplotlib.pyplot as plt

        # Set standard LaTeX font settings using built-in mathtext engine
        plt.rcParams.update({
            "text.usetex": False,            
            "mathtext.fontset": "cm",        # Force native Computer Modern LaTeX font
            "font.family": "serif",          # Match typical serif document layouts
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        })

        figWidthPTS = 246.0
        figWidthINCHES = figWidthPTS / 72.27
        heightScaleFactor = 0.8
        figHeightINCHES = figWidthINCHES * heightScaleFactor


        plt.figure(figsize=(figWidthINCHES, figHeightINCHES))
        plt.plot(
            data['time'],
            data['rate'],
            label='Original Data',
            color='lightgray'
            )
        plt.plot(
            data['time'],
            data['denoisedRate'],
            label='Denoised Data',
            color='Black',
            linewidth=1.0
            )
        plt.xlabel('Time (s)')
        plt.ylabel('Rate (counts/s)')
        plt.tight_layout()
        plt.xlim(0,60)
        plt.ylim(-8, 17)

        plt.show()


    data: LightCurveData = LightCurveData("GRB080319B", "swiftBAT")
    data.burstData['denoisedRate'] = haarDenoise(data.burstData)
    data.preBurstData['denoisedRate'] = haarDenoise(data.preBurstData)
    plotlightCurve(data.burstData)