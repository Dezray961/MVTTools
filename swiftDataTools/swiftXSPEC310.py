"""
Process a PHA file and a RSP file with a re-parameterised Band function and output results to a CSV file.
Takes input from environment variables for slice ID, target directory, and output CSV file.
The script is designed to be run in a parallel processing environment where each slice of data is
processed independently, and results are appended to a shared CSV file in a thread-safe manner.
It MUST be run in a Python 3.10 environment with XSPEC installed and configured correctly. I recommend
using a dedicated conda environment for this.
"""

import os, sys, csv, xspec


if __name__ == "__main__":
    sliceId: int = int(os.environ["SLICE_ID"])
    targetDir: str = os.environ["TARGET_DIR"]

    # convert the incoming output file path to an absolute path
    outputCsv: str = os.path.abspath(os.environ["OUTPUT_CSV"])

    phaFile: str = "outputSpectrum.pha"
    rspFile: str = "outputResponse.rsp"

    # global XSPEC environment flags
    xspec.Plot.device = "none"
    xspec.Xset.chatter = 0
    xspec.Xset.logChatter = 0
    xspec.Fit.query = "yes"

    try:
        os.chdir(targetDir)
        xspec.AllData.clear()
        xspec.AllModels.clear()

        if not os.path.exists(phaFile) or not os.path.exists(rspFile):
            print(f"CRITICAL: Required files missing in {targetDir}", flush=True)
            sys.exit(1)

        # load spectrum and apply criteria
        currentSpec: xspec.Spectrum = xspec.Spectrum(dataFile=phaFile, respFile=rspFile)
        currentSpec.ignore("0.0-15.0 150.0-**")
        xspec.AllModels.systematic = 0.02

        # configure grbm model using 1-based index notation
        bandModel: xspec.Model = xspec.Model("grbm")
        bandModel(1).values = -1.0      # alpha
        bandModel(2).values = -2.3      # beta
        bandModel(2).frozen = True      # Freeze beta index
        bandModel(3).values = [50.0, 1.0, 1.0, 1.0, 1000.0, 1500.0]  # Enforce hard upper limit on E0

        # fit
        xspec.Fit.nIterations = 100
        xspec.Fit.perform()

        # calculate error bounds
        xspec.Fit.error("1 3")

        # gather metrics
        reducedChiSq: float = xspec.Fit.statistic / xspec.Fit.dof if xspec.Fit.dof > 0 else 0.0
        
        alphaVal: float = bandModel(1).values[0]
        alphaLow, alphaHigh, _ = bandModel(1).error
        
        temVal: float = bandModel(3).values[0]
        temLow, temHigh, _ = bandModel(3).error

        # Epeak = E0 * (2 + alpha)
        epeakVal: float = temVal * (2.0 + alphaVal)
        epeakLow: float = temLow * (2.0 + alphaLow)
        epeakHigh: float = temHigh * (2.0 + alphaHigh)
        statusText: str = "Success"

    except Exception as fittingError:
        reducedChiSq, epeakVal, epeakLow, epeakHigh = 0.0, 0.0, 0.0, 0.0
        statusText: str = f"Failed: {str(fittingError)}"
        print(f"FITTING EXCEPTION: {statusText}", flush=True)

    # save to row
    rowResult = {
        "sliceId": sliceId, 
        "status": statusText, 
        "reducedChiSq": reducedChiSq,
        "epeak": epeakVal, 
        "epeakLow": epeakLow, 
        "epeakHigh": epeakHigh
    }

    fileNeedsHeader = not os.path.exists(outputCsv)
    with open(outputCsv, "a", newline="") as csvFileObject:
        dictWriter = csv.DictWriter(csvFileObject, fieldnames=rowResult.keys())
        if fileNeedsHeader:
            dictWriter.writeheader()
        dictWriter.writerows([rowResult])
        
    print("WORKER COMPLETE: Row appended.", flush=True)