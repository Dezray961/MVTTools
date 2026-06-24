"""
This needs to be split into a terminal running script with a class that can be imported into other
scripts. It also needs to use 'with chdir()' to change the working directory rather than using 
os as it will automatically change back to the original working directory when the block is exited.
Further, the SwiftBAT tools should be in their own directory.
"""

import swiftDataTools.swiftBATCatalogueGRB as catalogue
import swiftDataTools.swiftBATDataFetcher as fetcher
import swiftDataTools.swiftBATInitalProcessing as processing
from shellTools.runShell import ShellRunner
from contextlib import chdir
from pathlib import Path
import pandas as pd
import shutil
import os

# function to process the data into a usable format for analysis
def processSwiftBATData(
        GRBName: str,
        SNRThreshold: float,
        energyRange: str,
        timeBinSize: float = 0.0,
        download: bool = False,
        deleteOriginal: bool = False,
        fullDataSet: bool = False
        ) -> None:
    """Processes the data for a given GRB name into a usable format for analysis. This includes downloading the data from the swift.ac.uk archive, unzipping the data if it is compressed, and using the HEASoft tools to process the data into light curves for different energy ranges.

    Args:
        GRBName (str): The name of the GRB to be processed. Fomatted as "GRBYYMMDDX", where X is the letter assigned to the GRB. 
        SNRThreshold (float): The signal-to-noise ratio threshold to be used in the light curves bins. This will determine the size of the time bins in the light curves, with higher SNR thresholds resulting in larger time bins. If timeBinSize is provided, this will be ignored and the time bins will be of the specified size instead.
        energyRange (str): The energy range to be used in the light curves. This should be in the format "min-max", where min and max are the minimum and maximum energies in keV. The BAT sensor has a energy range of 15-350 keV.
        timeBinSize (float): The size of the time bins, in seconds, to be used in the light curves.
        download (bool): Whether to download the data from the swift.ac.uk archive. If False, it is assumed that the data has already been downloaded and is available in the "data/reproc" directory.
        deleteOriginal (bool): Whether to delete the original data file after processing. If True, the original data file will be deleted after the processed data has been saved.

    Returns:
        None: This function does not return anything, but it will save the processed data to a file that can be used for analysis. The processed data will include light curves for the specified energy range. The file will be saved in the same directory as the original data, with the same name but with "_processed" appended to the end of the file name.
    """

    # function to convert the light curve data into a CSV file
    def convertLightCurveToCSV(
            filename: str,
            GRBName: str
            ) -> None:
        # function to run fdump
        def runFdump(
                shell: ShellRunner,
                filename: str
                ) -> None:
            # run the command to convert the light curve data into a CSV file
            fdumpCommand: str = f'fdump {filename}.lc outfile={filename}.txt prhead=no clobber=yes columns="*" rows="-"'
            shellOutput: list[str] = shell.runShellCommand(fdumpCommand)
            for line in shellOutput:
                print(line)



        # function to read in the data, find the gap lines and return the data blocks as pandas DataFrames
        def linesToDataFrames(
                directory: str,
                filename: str
                ) -> tuple[pd.DataFrame, pd.DataFrame]:
            print(os.getcwd())
            with open(f'{directory}/{filename}.txt', 'r') as outputFile:
                lines: list[str] = outputFile.readlines()
            # find the gap lines
            gapLines: list[int] = [i for i, line in enumerate(lines) if not line.strip()]
            # The first two are blank, the data rows are between the 2nd and 3rd, and the 3rd and 4th
            block1Lines: list[str] = lines[gapLines[1]:gapLines[2]]
            block2Lines: list[str] = lines[gapLines[2]:gapLines[3]]
            # clean the data and convert it into pandas DataFrames
            cleanedBlock1Lines: pd.DataFrame = linesToDataFrame(block1Lines)
            cleanedBlock2Lines: pd.DataFrame = linesToDataFrame(block2Lines)
            return cleanedBlock1Lines, cleanedBlock2Lines


        # function to clean the data and convert it into a pandas DataFrame
        def linesToDataFrame(lines: list[str]) -> pd.DataFrame:
            data = [line.strip().split() for line in lines if line.strip()]
            # change the column names to lower case and remove the units
            data[0] = [col.split('(')[0].lower() for col in data[0]]
            data.pop(1)
            # seperate the index column into a separate list
            indexList = [row[0] for row in data[1:]]
            for i, list in enumerate(data):
                if i == 0:
                    continue
                list.pop(0)
            
                    
            return pd.DataFrame(data[1:], columns=data[0], index=indexList)


        # split the filename to get the directory
        directory: str = filename.rsplit("/", 1)[0]
        # open a shell
        with chdir(directory):
            shell: ShellRunner = ShellRunner()
            # run the command to convert the pre-burst light curve data into a CSV file
            runFdump(shell, f"outputPB")
            # run the command to convert the burst light curve data into a CSV file
            runFdump(shell, f"output")
            shell.closeShell()
            
        
        # create a new file path for the CSV file
        Path('data/processed').mkdir(parents=True, exist_ok=True)
        csvFilePath: str = f"data/processed/{GRBName}LC.csv"


        #############################################################################
        # change this to save the pre-burst data to CSV as well. It no lnoger       #
        # needs to calculate the standard deviation of the pre-burst count rate,    #
        #############################################################################

        # burst
        # read the output.txt file and write the data to the new CSV file
        burstData1, burstData2 = linesToDataFrames(directory, "output")

        # combine the two blocks of data side by side, using the index to align the rows
        mergedData: pd.DataFrame = pd.concat([burstData1, burstData2], axis=1)

        # write the merged data to a CSV file
        mergedData.to_csv(csvFilePath, index=False)

        # pre-burst
        # read the pre-burst output.txt file and calculate the standard deviation of the count rate
        preBurstData1, preBurstData2 = linesToDataFrames(directory, "outputPB")

        # combine the two blocks of pre-burst data side by side, using the index to align the rows
        preBurstMergedData: pd.DataFrame = pd.concat([preBurstData1, preBurstData2], axis=1)

        # truncate the pre-burst data to the same length as the burst data
        preBurstMergedData = preBurstMergedData.iloc[:len(mergedData)]

        # write the merged pre-burst data to a CSV file
        preBurstMergedData.to_csv(f"data/processed/{GRBName}PB.csv", index=False)


    # generate the wget statement to download the data for the given GRB name
    wgetStatement: str = fetcher.generateWgetStatement(GRBName)

    # download the data using the generated wget statement
    if download:
        print(f"Downloading data for {GRBName}...")
        fetcher.downloadSwiftBATData(wgetStatement)
        print("Download complete.")

    # generate the file path for the downloaded data
    data, _ = catalogue.importData('summary_general.csv')
    observationID: str = catalogue.getObservationID(GRBName, data)
    dataFilePath: str = f"data/reproc/{observationID}/bat/event/sw{observationID}bevshsp_uf.evt.gz"

    # process the data using the HEASoft tools to generate light curves for the specified energy range and time bin size
    print(f"Processing data for {GRBName}...")
    processing.processSwiftBATData(
        filename = dataFilePath,
        SNRThreshold = SNRThreshold,
        energyBins = energyRange,
        GRBName = GRBName,
        timeDeliniation = timeBinSize,
        fullDataSet = fullDataSet
    )
    print("Processing complete.")

    # convert the processed data into a CSV file for analysis
    print(f"Converting processed data for {GRBName} into CSV format...")
    convertLightCurveToCSV(dataFilePath, GRBName)


    print("Conversion complete.")

    # delete the original data file if specified
    if deleteOriginal:
        print(f"Deleting original data file for {GRBName}...")
        shutil.rmtree(f'data/reproc/{observationID}')
        print("Deletion complete.")

if __name__ == "__main__":
    grbName: str = "GRB080319B"
    processSwiftBATData(
        GRBName = grbName,
        SNRThreshold = 5.0,
        energyRange = "15-350",
        download = False,
        deleteOriginal = False,
        #timeBinSize = 100e-6,
        fullDataSet = False
        )
