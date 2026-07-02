from analysisTools.importLightCurve import LightCurveData
from analysisTools.haarDenoise import haarDenoise



class MonteCarloUncertainty:
    """
    A class to calculate the uncertainty of the denoised light curve using a parametric Monte Carlo method. The uncertainty is calculated by generating synthetic light curves based on the original light curve data and its propagated error, and then calculating the denoised values for each synthetic light curve. The standard deviation of the denoised values from the synthetic light curves is used as the uncertainty of the denoised light curve.

    Attributes:
        data (LightCurveData): The light curve data for which to calculate the uncertainty.
        numSimulations (int): The number of synthetic light curves to generate for the Monte Carlo simulation.
    """

    def __init__(
            self,
            data: LightCurveData,
            numSimulations: int = 1000
        ) -> None:
        """
        Initializes the MonteCarloUncertainty class with the given light curve data and number of simulations.

        Args:
            data (LightCurveData): The light curve data for which to calculate the uncertainty.
            numSimulations (int, optional): The number of synthetic light curves to generate for the Monte Carlo simulation. Defaults to 1000.
        """
        self.data: LightCurveData = data
        self.numSimulations: int = numSimulations
    

    