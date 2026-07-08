"""
4. Window the light curve data to a specific time range
5. Find the undecimated Haar transform of each window
6. Find the MVT of each window
7. MVT(t)
"""

from analysisTools.importLightCurve import LightCurveData
from analysisTools.rebinLightCurve import rebinLightCurve
from analysisTools.pyramidsDWTs import MODWT
from analysisTools.windowData import getWindowIndices
from analysisTools.parametricMCUncertainty import MonteCarloUncertainty
from pandas import DataFrame
from tqdm import tqdm


def analysisPipe(
        GRBName: str,
        timeWindowSize: float = 0.3,
        energyRange: str = "15-350",
        numSimulations: int = 1000,
        highRAMSystem: bool = True
    ) -> None:
    """This function is the main analysis pipeline for a given GRB. It imports the light curve data, denoises the data, rebins the data to a constant SNR, windows the data to a specific time range, finds the undecimated Haar transform of each window, and finds the MVT of each window."""
    

    # Do I want this to also import?

    # import the light curve data for the given GRB name and energy range
    data: LightCurveData = LightCurveData(GRBName)
    # denoise the light curve data using the Haar wavelet transform
    MonteCarloUncertainty(
        data = data,
        denoisedDataArgs = {
            "thresholdMethod": "hard",
            "thresholdScaleFactor": 0.5
        },
        numSimulations = numSimulations,
        highRAMSystem = highRAMSystem
    )
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
    from loadConfig import importConfiguration
    config = importConfiguration()

    grbName: str = "GRB080319B"
    analysisPipe(
        GRBName = grbName,
        numSimulations = 500,
        highRAMSystem = False
        )