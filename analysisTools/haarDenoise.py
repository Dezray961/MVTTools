from numpy import ndarray, zeros, log, mean, floor, log2, sign
from analysisTools.importLightCurve import LightCurveData
from analysisTools.pyramidsDWTs import MODWT, inverseMODWT



def haarDenoise(
        data: LightCurveData,
        thresholdMethod: str = "soft"
        ) -> None:
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
        firstLogTerm: float = log(n_l)
        secondLogTerm: float = log(n_l * meanRate)
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
        preBurstCounts: ndarray = data.burstData['totcounts'].to_numpy()
        timeInBin: float = 100e-6 * binSize
        countsInBin: ndarray = zeros(len(preBurstCounts) - binSize)
        for i in range(len(preBurstCounts) - binSize):
            countsInBin[i] = sum(preBurstCounts[i:i + binSize])
        meanRate: float = mean(countsInBin) / timeInBin
        return meanRate



    burstData: list[float] = data.burstData['rate'].to_list()
    burstError: list[float] = data.burstData['error'].to_list()

    maxLevel: int = int(floor(log2(len(burstData)))) - 1

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
                    waveletCoeffs[level][i] = sign(waveletCoeffs[level][i]) * (abs(waveletCoeffs[level][i]) - thresholdValue)
        elif thresholdMethod == "hard":
            for i in range(len(waveletCoeffs[level])):
                if abs(waveletCoeffs[level][i]) < thresholdValue:
                    waveletCoeffs[level][i] = 0
        else:
            raise ValueError(f"Invalid threshold method: {thresholdMethod}. Must be 'soft' or 'hard'.")


    # invert the transform to get the denoised data
    denoisedData: ndarray = inverseMODWT(
        waveletCoeffs,
        scalingCoeffs
        )
    data.burstData['denoisedRate'] = denoisedData


if __name__ == "__main__":
    def plotlightCurve(
        data: LightCurveData
        ) -> None:
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
    data: LightCurveData = LightCurveData("GRB080319B")
    haarDenoise(data)
    plotlightCurve(data)