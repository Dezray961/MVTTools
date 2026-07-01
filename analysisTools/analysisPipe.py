"""
3. Rebin the light curve data to constant SNR
4. Window the light curve data to a specific time range
5. Find the undecimated Haar transform of each window
6. Find the MVT of each window
7. MVT(t)
"""

from analysisTools.importLightCurve import LightCurveData
from analysisTools.haarDenoise import haarDenoise
from analysisTools.rebinLightCurve import rebinLightCurve


def analysisPipe(
        GRBName: str,
        energyRange: str = "15-350"
    ) -> None:
    """This function is the main analysis pipeline for a given GRB. It imports the light curve data, denoises the data, rebins the data to a constant SNR, windows the data to a specific time range, finds the undecimated Haar transform of each window, and finds the MVT of each window."""
    
    # Do I want this to also import?

    # import the light curve data for the given GRB name and energy range
    data: LightCurveData = LightCurveData(GRBName)
    # denoise the light curve data using the Haar wavelet transform
    haarDenoise(data)
    # rebin the light curve data to a constant SNR
    rebinLightCurve(data, instrument = "swiftBAT", snrThreshold = 5.0)
    # window the light curve data to a specific time range

    # find the undecimated Haar transform of each window

    # find the MVT of each window

    # MVT(t)


if __name__ == "__main__":
    grbName: str = "GRB080319B"
    analysisPipe(
        GRBName = grbName,
        energyRange = "15-350"
        )