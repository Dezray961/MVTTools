import numpy as np
from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT, inverseMODWT



def haarDenoise(
        data: LightCurveData,
        thresholdMethod: str = "soft"
        ) -> LightCurveData:
    """Denoises the light curve data using the Haar wavelet transform and thresholding. The denoised data is added to the LightCurveData object as a new column.

    Args:
        data (LightCurveData): light curve data to denoise
        thresholdMethod (str, optional): thresholding method to use. Defaults to "soft". Options are "soft" and "hard". 

    Returns:
        LightCurveData: light curve data with denoised data added as a new column
    """
    def threshold(
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
        firstLogTerm: float = np.log(n_l)
        secondLogTerm: float = np.log(n_l * meanRate)
        sqrtTerm: float = max((((4 * firstLogTerm) ** 2) + (8 * secondLogTerm)), 0)
        thresholdValue: float = prefactor * (2 * firstLogTerm + sqrtTerm ** 0.5)
        return thresholdValue


    def findPoissonRate(
            data: LightCurveData,
            binSize: int
            ) -> float:
        """Finds the mean rate of the pre-burst data for a given bin size.

        Args:
            data (LightCurveData): light curve data
            binSize (int): size of the bins in indices.

        Returns:
            float: mean rate of the pre-burst data
        """
        preBurstCounts: np.ndarray = data.burstData['totcounts'].to_numpy()
        timeInBin: float = 100e-6 * binSize
        countsInBin: np.ndarray = np.zeros(len(preBurstCounts) - binSize)
        for i in range(len(preBurstCounts) - binSize):
            countsInBin[i] = np.sum(preBurstCounts[i:i + binSize])
        meanRate: float = np.mean(countsInBin) / timeInBin
        return meanRate



    burstData: list[float] = data.burstData['rate'].to_list()
    burstError: list[float] = data.burstData['error'].to_list()

    maxLevel: int = int(np.floor(np.log2(len(burstData)))) - 1

    # get the Haar wavelet coefficients and scaling coefficients using the MODWT
    waveletCoeffs, scalingCoeffs = MODWT(burstData)

    # apply thresholding to the wavelet coefficients
    for level in range(maxLevel):
        # calculate the threshold for the current level
        meanRate: float = findPoissonRate(data, 2 ** level)
        thresholdValue: float = threshold(level, maxLevel, meanRate)

        # apply the thresholding to the level's wavelet coefficients
        if thresholdMethod == "soft":
            for i in range(len(waveletCoeffs[level])):
                if abs(waveletCoeffs[level][i]) < thresholdValue:
                    waveletCoeffs[level][i] = 0
                else:
                    waveletCoeffs[level][i] = np.sign(waveletCoeffs[level][i]) * (abs(waveletCoeffs[level][i]) - thresholdValue)
        elif thresholdMethod == "hard":
            for i in range(len(waveletCoeffs[level])):
                if abs(waveletCoeffs[level][i]) < thresholdValue:
                    waveletCoeffs[level][i] = 0
        else:
            raise ValueError(f"Invalid threshold method: {thresholdMethod}. Must be 'soft' or 'hard'.")


    # invert the transform to get the denoised data
    denoisedData: np.ndarray = inverseMODWT(
        waveletCoeffs,
        scalingCoeffs
        )
    data.burstData['denoisedRate'] = denoisedData


    return data


def plotlightCurve(
        data: LightCurveData
        ) -> None:
    """Plots the light curve data

    Args:
        data (LightCurveData): light curve data to plot
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.plot(
        data.burstData['time'],
        data.burstData['rate'],
        label='Original Data'
        )
    plt.plot(
        data.burstData['time'],
        data.burstData['denoisedRate'],
        label='Denoised Data',
        color='red',
        alpha=0.5
        )
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (counts/s)')
    plt.legend()
    plt.show()


def getData():
    return LightCurveData("GRB080319B")


data = getData()
haarDenoise(data)
plotlightCurve(data)


#    
#    
#
#
#
#    # wavelet transform of the data. THis needs to be replaced with some C++ code as it is too slow.
#    coeffs: list = pywt.swt(
#        paddedData,
#        'haar'
#        )
#    print(f"coeffs: {coeffs}")
#    return np.array(coeffs)
    # thresholding of the wavelet coefficients

    # inverse wavelet transform to get the denoised data



#if __name__ == "__main__":
#    grbName: str = "GRB080319B"
#    data: LightCurveData = LightCurveData(grbName)
#    denoisedData: np.ndarray = haarDenoise(
#        data
#        )

