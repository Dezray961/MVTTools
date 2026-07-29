"""
Process a single PHA file with a re-parameterised Band function and output results to a CSV file.
Takes input from environment variables for slice ID, PHA file, RSP file, and output CSV file.
The script is designed to be run in a parallel processing environment where each slice of data is
processed independently, and results are appended to a shared CSV file in a thread-safe manner.
It MUST be run in a Python 3.10 environment with XSPEC installed and configured correctly. I recommend
using a dedicated conda environment for this.
"""

import os, csv, xspec

if __name__ == "__main__":
    # Pull parameters directly from the OS environment variables
    sliceId = int(os.environ["SLICE_ID"])
    phaFile = os.environ["PHA_FILE"]
    rspFile = os.environ["RSP_FILE"]
    outputCsv = os.environ["OUTPUT_CSV"]

    # Suppress plots and heavy log output
    xspec.Plot.device = "none"
    xspec.Xset.chatter = 1
    xspec.Xset.logChatter = 1

    try:
        # Load data and apply BAT constraints
        currentSpec = xspec.Spectrum(phaFile)
        currentSpec.response = rspFile
        currentSpec.ignore("0.0-15.0 150.0-**")
        currentSpec.systematic = 0.02

        # Setup re-parameterised Band function
        bandModel = xspec.Model("bndrep")
        bandModel.alpha.values = -1.0
        bandModel.beta.values = -2.3
        bandModel.beta.frozen = True
        bandModel.Epeak.values = 70.0

        # Fit and calculate 90% error bounds on Epeak
        xspec.Fit.nIterations = 100
        xspec.Fit.perform()
        xspec.Fit.error("3")

        # Collect metrics
        reducedChiSq = xspec.Fit.statistic / xspec.Fit.dof if xspec.Fit.dof > 0 else 0.0
        epeakVal = bandModel.Epeak.values
        epeakLow = bandModel.Epeak.error
        epeakHigh = bandModel.Epeak.error
        statusText = "Success"

    except Exception as fittingError:
        reducedChiSq, epeakVal, epeakLow, epeakHigh = 0.0, 0.0, 0.0, 0.0
        statusText = f"Failed: {str(fittingError)}"

    # Thread-safe appending to CSV
    rowResult = {
        "sliceId": sliceId, "status": statusText, "reducedChiSq": round(reducedChiSq, 2),
        "epeak": round(epeakVal, 2), "epeakLow": round(epeakLow, 2), "epeakHigh": round(epeakHigh, 2)
    }

    fileNeedsHeader = not os.path.exists(outputCsv)
    with open(outputCsv, "a", newline="") as csvFileObject:
        dictWriter = csv.DictWriter(csvFileObject, fieldnames=rowResult.keys())
        if fileNeedsHeader:
            dictWriter.writeheader()
        dictWriter.writerows([rowResult])
