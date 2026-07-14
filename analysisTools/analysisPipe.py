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
        source: str
    ) -> None:
    """This function is the main analysis pipeline for a given GRB. It imports the light curve data, denoises the data, rebins the data to a constant SNR, windows the data to a specific time range, finds the undecimated Haar transform of each window, and finds the MVT of each window."""
    

    # Do I want this to also import?

    # import the light curve data for the given GRB name and energy range
    data: LightCurveData = LightCurveData(GRBName)



if __name__ == "__main__":
    from loadConfig import importConfiguration
    config = importConfiguration()

    grbName: str = "GRB080319B"
    analysisPipe(
        GRBName = grbName,
        source = "swift"
        )