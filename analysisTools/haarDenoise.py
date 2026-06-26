import pywt
import numpy as np
from analysisTools.importLightCurve import LightCurveData


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
        prefactor: float = (2 ** (-0.5 * (level + 2)))
        firstLogTerm: float = np.log(n_l)
        secondLogTerm: float = np.max([np.log(n_l * meanRate), 0])
        thresholdValue: float = prefactor * (2 * firstLogTerm + (((4 * firstLogTerm) ** 2) + (8 * secondLogTerm)) ** 0.5)
        return thresholdValue
    

    burstData: list[float] = data.burstData['rate'].to_list()
    burstError: list[float] = data.burstData['error'].to_list()

    # find the mean rate of the pre-burst data
    maxLevel: int = int(np.floor(np.log2(len(burstData))))
    # find the wavelet coefficients of the pre-burst data. This produces a (J - L + 1) x N matrix
    # this does not need to be the translation invariant wavelet transform, so it uses the DWT from pywavelets

    haarTransformCoefficients: np.ndarray = pywt.wavedec(
        burstData,
        'haar',
        mode = 'zero'
        )

    # find L - this may be used later, currently the wavelet transform is done over the whole space
    L: int = len(haarTransformCoefficients) - maxLevel - 1
    # find the threshold for each level of the wavelet transform
    print(haarTransformCoefficients[2][2])
    for level in range(1, maxLevel - L):
        meanRate: float = np.abs(haarTransformCoefficients[level][0] / ((2 ** level) * 100e-6))
        thresholdValue: float = threshold(
            level,
            maxLevel,
            meanRate
            )
        for i, coefficient in enumerate(haarTransformCoefficients[level][1:]):
            if thresholdMethod == "soft":
                if np.abs(coefficient) < 3 * thresholdValue:
                    haarTransformCoefficients[level][i + 1] = 0
                else:
                    haarTransformCoefficients[level][i + 1] = np.sign(coefficient) * (np.abs(coefficient) - 3 * thresholdValue)
            elif thresholdMethod == "hard":
                if np.abs(coefficient) < thresholdValue:
                    haarTransformCoefficients[level][i + 1] = 0
            else:
                raise ValueError(f"Invalid thresholding method: {thresholdMethod}. Must be 'soft' or 'hard'.")
    print(haarTransformCoefficients[2][2])
    # inverse wavelet transform to get the denoised data
    denoisedData: np.ndarray = pywt.waverec(
        haarTransformCoefficients,
        'haar',
        mode = 'zero'
        )


    data.burstData['denoisedRate'] = denoisedData[:len(data.burstData)]

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
        color='red'
        )
    plt.xlabel('Time (s)')
    plt.ylabel('Rate (counts/s)')
    plt.legend()
    plt.show()





#    # pad the data to the next power of 2
#    n: int = int(np.ceil(np.log2(len(data))))
#    nextPowerOf2: int = int(2 ** n)
#    padSize: int = (nextPowerOf2 - len(data)) // 2
#    paddedData: np.ndarray = np.pad(data, (padSize, padSize), mode='constant')
#    # check the length of the padded data is the correct length
#    difference: int = len(paddedData) - nextPowerOf2
#    match difference:
#        case _ if difference < 0:
#            # pad the data to the next power of 2
#            paddedData: np.ndarray = np.pad(paddedData, (0, -difference), mode='constant')
#        case _ if difference > 0:
#            # truncate the data to the next power of 2
#            paddedData: np.ndarray = paddedData[:nextPowerOf2]
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

