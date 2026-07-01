"""
4. Window the light curve data to a specific time range
5. Find the undecimated Haar transform of each window
6. Find the MVT of each window
7. MVT(t)
"""

from analysisTools.importLightCurve import LightCurveData
from analysisTools.haarDenoise import haarDenoise
from analysisTools.rebinLightCurve import rebinLightCurve
from analysisTools.pyramidsDWTs import MODWT, inverseMODWT
from pandas import DataFrame
from tqdm import tqdm


def analysisPipe(
        GRBName: str,
        timeWindowSize: float = 0.3,
        energyRange: str = "15-350"
    ) -> None:
    """This function is the main analysis pipeline for a given GRB. It imports the light curve data, denoises the data, rebins the data to a constant SNR, windows the data to a specific time range, finds the undecimated Haar transform of each window, and finds the MVT of each window."""
    
    def getWindowIndices(
            data: DataFrame,
            timeWindowSize: float
            ) -> list:
        """This function returns the indices of the light curve data that fall within a specific time window size.

        Args:
            data (DataFrame): the denoised light curve data
            timeWindowSize (float): the size of the time window to consider

        Returns:
            list: a list of indices corresponding to the time values that fall within the specified time window size
        """
        # get the timeInBins column from the data
        timeInBins = data['timeInBin'].to_numpy()

        # initialize variables to keep track of the window indices and the accumulated time
        windowIndices: list[tuple[int, int]] = []
        accumulatedTime: float = 0.0
        startIndex: int = 0
        endIndex: int = 0

        # get the number of time bins
        n = len(timeInBins)

        # iterate through the time bins and accumulate the time until it exceeds the time window size
        while startIndex < n:
            while accumulatedTime < timeWindowSize and endIndex < n:
                accumulatedTime += timeInBins[endIndex]
                endIndex += 1
            # add the window indices to the list and reset the accumulated time and start index
            if accumulatedTime >= timeWindowSize:
                windowIndices.append((startIndex, endIndex - 1))
                accumulatedTime -= timeInBins[startIndex]
                startIndex += 1
            else:
                break
        
        return windowIndices
    

    # Do I want this to also import?

    # import the light curve data for the given GRB name and energy range
    data: LightCurveData = LightCurveData(GRBName)
    # denoise the light curve data using the Haar wavelet transform
    haarDenoise(data)
    # rebin the light curve data to a constant SNR
    rebinLightCurve(data,
                    instrument = "swiftBAT",
                    snrThreshold = 5.0
                    )
    # window the light curve data to a specific time range
    windowIndices: list[tuple[int, int]] = getWindowIndices(
        data.rebinnedData,
        timeWindowSize
        )
    # find the undecimated Haar transform of each window
    for startIndex, endIndex in tqdm(windowIndices, desc="Processing Windows", unit="window"):
        windowData: DataFrame = data.rebinnedData.iloc[startIndex:endIndex + 1]
        # perform the undecimated Haar transform on the windowed data
        waveletCoefficients, scaleCoefficients = MODWT(windowData['rate'].to_numpy())
        # (this is a placeholder for the actual implementation)
        # haarTransform(windowData)

        # find the MVT of each window

        # MVT(t)
    print(waveletCoefficients, scaleCoefficients)


if __name__ == "__main__":
    grbName: str = "GRB080319B"
    analysisPipe(
        GRBName = grbName,
        energyRange = "15-350"
        )