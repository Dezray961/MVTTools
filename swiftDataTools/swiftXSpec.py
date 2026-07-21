"""Needs to be run in a python 3.10 environment with xspec installed."""
import sys, os, xspec

def spectralFit(
        phaFile: str,
        model: str
        ) -> dict:
    """
    Perform spectral fitting using XSPEC.

    Parameters:
    phaFile (str): The path to the PHA file containing the spectral data.

    Returns:
    dict: A dictionary containing the fit results and statistics.
    """
    # Load the data
    print(f"Loading PHA file: {phaFile}")
    xspec.AllData.clear()
    xspec.AllData(phaFile)

    print(f"Setting model: {model}")
    # Set the model
    xspec.Model(model)

    print("Ignoring energy bins outside the range 15-150 keV")

    # Swift BAT has a minimum energy of 15 keV and over 150 keV the mask is transparent.
    # remove bins outside the "good" energy range of Swift (< 15 keV and > 150 keV)
    xspec.AllData.ignore("0.0-14.0, 150.0-**")

    print("Performing the fit...")
    # Perform the fit
    xspec.Fit.query = "yes"  # Automatically answer 'yes' to any prompts during fitting
    xspec.Fit.perform()

    print("Fit completed.")

    print(xspec.AllModels(1))

#    # Retrieve fit results
#    fitResults = {
#        'parameters': [param.values for param in xspec.AllModels(1).parameters],
#        'errors': [param.error for param in xspec.AllModels(1).parameters],
#        'chiSquare': xspec.Fit.statistic,
#        'degreesOfFreedom': xspec.Fit.dof
#    }
#
#    return fitResults


def getCleanArgs():
    # detect if the script is being run in an interactive environment (like Jupyter Notebook)
    isJupyter = 'ipykernel' in sys.modules or any("-f" in arg for arg in sys.argv)

    if isJupyter:
        # If running in Jupyter, return a default set of arguments for testing
        return ["swiftXSpec.py", "outputSpectrum.pha", "grbm", "test"]
    else:
        # If running from the command line, return the actual command-line arguments
        return sys.argv





def main():
    args = getCleanArgs()
    if len(args) < 3:
        print("Usage: python swiftXSpec.py <phaFile> <model>")
        sys.exit(1)

    phaFile = args[1]
    model = args[2]

    if not os.path.isfile(phaFile):
        print(f"Error: PHA file '{phaFile}' does not exist.")
        sys.exit(1)

    fitResults = spectralFit(phaFile, model)
    print("Fit Results:")
    print(f"Parameters: {fitResults['parameters']}")
    print(f"Errors: {fitResults['errors']}")
    print(f"Chi-Square: {fitResults['chiSquare']}")
    print(f"Degrees of Freedom: {fitResults['degreesOfFreedom']}")


if __name__ == "__main__":
    sys.argv = getCleanArgs()  # Clean the arguments before proceeding
    if sys.argv[3] == "test":
        # contextlib.chdir is python 3.11+, so we change the directory manually for the test case
        originalDir = os.getcwd()
        try:
            dataPath = "data/reproc/00306757000/bat/event"
            os.chdir(dataPath)
            sys.argv[1] = "outputSpectrum.pha"  # Replace the test argument with a default PHA file for testing
            main()
        finally:
            os.chdir(originalDir)
    else:
        main()