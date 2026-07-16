from os import system, environ, getcwd, listdir, path as osPath
from heasoftpy import Config, fkeyprint, bateconvert, batmaskwtevt, batbinevt, fdump
from contextlib import chdir
from swiftDataTools.swiftBATCatalogueGRB import getObservationID, importData, getCoordinates, getStartStopTime
from pathlib import Path
from shutil import move
from astropy.io import fits

class ProcessSwiftData:
    """Class to process Swift BAT data for a given GRB. This class handles the processing of Swift BAT data for a given GRB, including checking and applying gain correction, checking and applying mask weighting, extracting light curves for different time periods, and converting the light curve data into CSV files. It has the ability to process a custom time range, however this will require the user to call methods from outside the class."""
    def __init__(
            self,
            GRBName: str,
            energyBins: str = "15-350"
            )-> None:
        """Constructor for `ProcessSwiftData` class. Runs a complete pipeline to convert data from Swift BAT raw data to a standardised CSV file. Requires the data to exist in `/data/reproc/{observationID}/bat/`

        Args:
            GRBName (str): Name of the GRB to process
            energyBins (str, optional): Range of energy bins to use. Should be a string in the format "min-max", with the units in keV. Defaults to "15-350" This is the maximum energy range of the BAT sensor, see the BAT Data Analysis Guide for more information. https://swift.gsfc.nasa.gov/analysis/bat_swguide_v6_3.pdf

        Raises:
            EnvironmentError: if the HEADAS environment variable is not set. This is required for the HEASoft tools to function properly. The user should ensure that they have initialized HEASoft before running this code.
            FileNotFoundError: if the required data files are not found in the specified directory.
            ValueError: if the provided GRB name is invalid or not found in the dataset.
            ValueError: if the provided energy bins are invalid or not supported.

        Returns:
            None: Output is saved to a CSV file in the `data/processed/{GRBName}` directory.
        """
        # Check if the HEADAS environment variable is set. If it is not set, print an error message and exit the program. This is important because the HEASoft tools require the HEADAS environment variable to be set in order to function properly. If the variable is not set, the program will not be able to find the necessary tools and will fail to run. By checking for the variable at the beginning of the program, we can ensure that the user is aware of the issue and can take steps to fix it before proceeding with the data processing.
        self.__headasPath = environ.get("HEADAS")
        if self.__headasPath:
            environ["PFILES"] = f"{environ['HOME']}/pfiles;{self.__headasPath}/syspfiles"
        else:
            raise EnvironmentError("Error: HEADAS environment variable not found. Did you initialize HEASoft?")
        Config.allow_failure = False

        self.GRBName = GRBName
        self.energyBins = energyBins
        self.__data, self.__columnNames = importData("swiftDataTools/summary_general.csv")
        self.__triggerID: str = getObservationID(
            GRBName,
            self.__data,
            isTrigID=False
            )

        # create the output directory if it does not exist
        Path(f"data/processed/{self.GRBName}").mkdir(parents=True, exist_ok=True)

        # process the data using the HEASoft tools
        with chdir(f"data/reproc/{self.__triggerID}/bat/event"):
            print(getcwd())
            print(f"Processing data for {GRBName}...")
            self.__eventFilename: str = f"sw{self.__triggerID}bevshsp_uf.evt.gz"
            
            # check if the gain correction has been applied to the event file.
            if not self.__checkGain():
                # correct the gain if it has not been applied.
                self.__correctGain()
            
            # check if the mask has been applied to the event file. 
            if not self.__checkMask():
                # apply the mask if it has not been applied.
                self.__applyMask()

            # get the start and stop times of the burst from the summary_general.csv file.
            self.__startTime, self.__stopTime, _ = getStartStopTime(
                self.GRBName,
                self.__data
            )

            # find the burst duration
            self.__burstDuration: float = self.__stopTime - self.__startTime
            print(f"Start time: {self.__startTime} seconds")
            print(f"Stop time: {self.__stopTime} seconds")
            print(f"Burst duration: {self.__burstDuration} seconds")

            # extract the light curves
            for period in range(2):
                self.extractLightCurve(period)

            # TODO
            # load the light curve data from the output files using astropy.io.fits

            # generate the spectrum for the burst period

            # revise the spectrum using batupdatephakw and batphasyserr

            # generate the response matrix using batdrmgen

            # might as well calculate the Epeak while the spectrum is being generated.
            # analyze the spectrum using XSPEC to find the Epeak

            # calculate the effective area of the detector (sum the elements of the response matrix)

            # convert the bins in the light curve to photon counts using the effective area photons = counts / effective area

            # potentially recalculate the light curve and uncertaities using monte carlo simulations?
            # this step is done by Bala et al 2026 on Fermi data.
            # to do this assume that the light curve is the poisson mean and get the background rate from
            # the pre-burst and post-burst light curves. See if parametricMCUncertainty.py can be used
            # for this.







    def __unzipFile(
            self,
            filename: str
        )-> None:
        """Unzips a file using the gunzip command. This is a wrapper function for the gunzip command, which is used to unzip files that have been compressed using the gzip algorithm. The function takes in the name of the file to be unzipped, and uses the gunzip command to unzip the file. The function does not return anything, but it will print out any output from the gunzip command to the console.

        Args:
            filename (str): The name of the file to be unzipped. This should be the name of the file that has been compressed using gzip.
        """
        print("Unzipping file: ", filename)
        system(f"gunzip {filename}")


    def __zipFile(
            self,
            filename: str
        )-> None:
        """Zips a file using the gzip command. This is a wrapper function for the gzip command, which is used to compress files using the gzip algorithm. The function takes in the name of the file to be zipped, and uses the gzip command to compress the file. The function does not return anything, but it will print out any output from the gzip command to the console.

        Args:
            filename (str): The name of the file to be zipped. This should be the name of the file that is to be compressed using gzip.
        """
        print("Zipping file: ", filename)
        system(f"gzip {filename} -v")


    def __checkGain(
            self
        )-> bool:
        """Checks the gain has been applied to BAT data. Uses `fkeyprint` to check the GAINAPP and GAINMETH keywords in the event file header.

        Args:
            filename (str): The name of the file to be checked. This should be the name of the event file. Assumes that the file is in the current working directory and is still compressed.

        Returns:
            bool: True if the gain correction has been applied, False if the gain correction has not been applied
        """
        print("Checking gain correction")
        output = fkeyprint(
            infile = self.__eventFilename,
            keynam = "GAIN"
            )
        print(output.stdout)
        bools:list[bool] = []
        for line in output.stdout.splitlines():
            if line.startswith("GAINAPP"):
                if "T" in line or "Gain correction has been applied" in line:
                    bools.append(True)
                else:
                    bools.append(False)
            if line.startswith("GAINMETH"):
                if "FIXEDDAC" in line:
                    bools.append(True)
                else:
                    bools.append(False)
        if all(bools):
            print("Gain correction already applied")
        
        return all(bools)


    def __correctGain(
            self
        )-> None:
        """Corrects the gain of BAT data.

        Currently needs to find the calibration file from the HEASARC FTP site.
        /swift/data/trend/YYYY_MM/bat/bgainoffs
        see page 31-32 of the BAT Data Analysis Guide https://swift.gsfc.nasa.gov/analysis/bat_swguide_v6_3.pdf
        Args:
            filename (str): The name of the file to be corrected. This should be the name of the event file. Assumes that the file is in the current working directory and is still compressed.

        Returns:
            None: This function does not return anything, but it will correct the gain of the BAT data in the specified file.
        """
        print("Gain correction has not been applied.")
        print("Applying gain correction to file: ", self.__eventFilename)
        # find the calibration file in the hk directory. 
        self.__calibrationFile: str = f"../hk/{self.__eventFilename.replace('bevshsp_uf.evt.gz', 'bcbo01deg00ab.fits.gz')}"
        if not osPath.exists(self.__calibrationFile):
            ### download the calibration file - impliment this later
            raise FileNotFoundError(f"Calibration file not found: {self.__calibrationFile}")


        # unzip the file if it is compressed. Needed for the bateconvert command to work. 
        if self.__eventFilename.endswith(".gz"):
            self.__unzipFile(self.__eventFilename)
            self.__eventFilename = self.__eventFilename.replace(".gz", "")

        # run the bateconvert command to correct the gain. 
        output = bateconvert(
            infile = self.__eventFilename,
            calfile = self.__calibrationFile,
            residfile = "CALDB",
            pulserfile = "CALDB",
            fitpulserfile = "CALDB",
            outfile = "NONE",
            calmode = "INDEF",
        )
        print(output.stdout)


        # zip the file back up
        self.__zipFile(self.__eventFilename)
        self.__eventFilename = self.__eventFilename + ".gz"


    def __checkMask(
            self
        )-> bool:
        
        def checkMaskWeighting(
            self
        )-> bool:
            """Checks the mask has been applied to BAT data. Uses the output from `fkeyprint` to check the MASKAPP and MASKMETH keywords in the event file header.

            Args:
                filename (str): The name of the file to be checked. This should be the name of the event file. Assumes that the file is in the current working directory and is still compressed.

            Returns:
                bool: True if the mask has been applied, False if the mask has not been applied
            """
            print("Checking mask weighting")
            output = fkeyprint(
                infile = self.__eventFilename,
                keynam = "BAT_"
            )
            print(output.stdout)
            bools:list[bool] = []
            for line in output.stdout.splitlines():
                if line.startswith("BAT_RA"):
                    splitLine: str = line.split('=')[1].strip().split('/')[0].strip()
                    try:
                        # do not use this value for the ra and dec as if the mask is incorrectly applied, these values will be incorrect.
                        float(splitLine)
                        bools.append(True)
                    except:
                        bools.append(False)
                if line.startswith("BAT_DEC"):
                    splitLine: str = line.split('=')[1].strip().split('/')[0].strip()
                    try:
                        float(splitLine)
                        bools.append(True)
                    except:
                        bools.append(False)
            if all(bools):
                print("Mask weighting already applied")
            
            return all(bools)


        def checkMaskVersion(
            self
        )-> bool:
            """Checks the `batmaskwtevt` version that has been applied to the data. Uses the output from `fkeyprint` to check the MASKVER keyword in the event file header.

            Args:
                output (FKeyPrintOutput): _description_

            Returns:
                bool: _description_
            """
            print("Checking batmaskwtevt version")
            output = fkeyprint(
                infile = self.__eventFilename,
                keynam = "BATCREAT"
            )
            print(output.stdout)
            bools:list[bool] = []
            for line in output.stdout.splitlines():
                if line.startswith("BATCREAT"):
                    if "batmaskwtevt" in line:
                        version: float = float(line.split('batmaskwtevt')[1].strip().split(' ')[0].strip('\''))
                        if version >= 1.16:
                            bools.append(True)
                            print(f"batmaskwtevt version is {version}>= 1.16") 
                        else:
                            bools.append(False)
                            print(f"batmaskwtevt version is {version}< 1.16")
                    else:
                        bools.append(False)
            return all(bools)


        bools: list[bool] = []
        bools.append(checkMaskWeighting(self))
        bools.append(checkMaskVersion(self))
        return all(bools)


    def __applyMask(
            self
        )-> None:
        """Applies the mask to BAT data. Uses `batmaskwtevt` to apply the mask to the event file.

        Args:
            filename (str): The name of the file to be masked. This should be the name of the event file. Assumes that the file is in the current working directory and is still compressed.
        """
        print("Mask has not been applied.")
        print("Applying mask to file: ", self.__eventFilename)

        # find the attitude file in the aux directory.
        self.__attitudeFile: str = f"../aux/{self.__eventFilename.replace('bevshsp_uf.evt.gz', 'sat.fits.gz')}"
        if not osPath.exists(self.__attitudeFile):
            ### download the attitude file - impliment this later
            raise FileNotFoundError(f"Attitude file not found: {self.__attitudeFile}")
            
        # find the quality map file in the hk directory.
        self.__qualityMapFile: str = f"../hk/{self.__eventFilename.replace('bevshsp_uf.evt.gz', 'bdqcb.hk.gz')}"
        if not osPath.exists(self.__qualityMapFile):
            ### download the quality map file - impliment this later
            raise FileNotFoundError(f"Quality map file not found: {self.__qualityMapFile}")

        # unzip the file if it is compressed. Needed for the batmaskwtevt command to work. 
        if self.__eventFilename.endswith(".gz"):
            self.__unzipFile(self.__eventFilename)
            self.__eventFilename = self.__eventFilename.replace(".gz", "")

        # find the right ascension and declination of the GRB from the summary_general.csv file.
        self.__rightAscension, self.__declination = getCoordinates(self.__triggerID, self.__data)

        # run the batmaskwtevt command to apply the mask.
        output = batmaskwtevt(
            infile = self.__eventFilename,
            attitude = self.__attitudeFile,
            ra = self.__rightAscension,
            dec = self.__declination,
            detmask = self.__qualityMapFile,
            rebalance = "YES",
            corrections = "default",
            auxfile = f"/local/data/gcn5b/craigm/{self.__triggerID}bevtr.fits",
            clobber = "YES"
        )
        print(output.stdout)

        # zip the file back up
        self.__zipFile(self.__eventFilename)
        self.__eventFilename = self.__eventFilename + ".gz"


    def extractLightCurve(
            self,
            period: int,
            customTimeRange: tuple[float, float] = None
        )-> None:
        """Extracts the lightcurve using `batbinevt`.

        Args:
            period (int): Time period to extract the light curve for. 0 = pre-burst, 1 = burst, 2 = post-burst, 3 = custom. If custom is selected, the user must provide `customTimeRange`
            customTimeRange (tuple[float, float], optional): Custom time range to extract a custom light curve for. Defaults to None.

        Raises:
            ValueError: If the period is not 0, 1, 2, or 3, or if the period is 3 and no custom time range is provided.
        """
        # determine the start and stop times for the light curve extraction
        match period:
            case 0: # pre-burst
                fileName: str = "outputPreBurst.lc"
                startTime: float = self.__startTime - self.__burstDuration
                stopTime: float = self.__startTime
                print("Extracting pre-burst uniform light curve")
            case 1: # burst
                fileName: str = "outputBurst.lc"
                startTime: float = self.__startTime
                stopTime: float = self.__stopTime
                print("Extracting burst uniform light curve")
            case 2: # post-burst
                fileName: str = "outputPostBurst.lc"
                startTime: float = self.__stopTime
                stopTime: float = self.__stopTime + self.__burstDuration
                print("Extracting post-burst uniform light curve")
            case 3: # custom
                if customTimeRange is None:
                    raise ValueError("Custom time range must be provided for period 3")
                fileName: str = "outputCustom.lc"
                startTime: float = customTimeRange[0]
                stopTime: float = customTimeRange[1]
                print("Extracting custom uniform light curve")
            case _:
                raise ValueError("Invalid period. Must be 0 (pre-burst), 1 (burst), 2 (post-burst), or 3 (custom).")
        
        # run batbinevt to create the light curve (uniform bins)
        output = batbinevt(
            infile = self.__eventFilename,
            outfile = fileName,
            outtype = "LC",
            timedel = 100e-6,
            timebinalg = "u",
            energybins = self.energyBins,
            detmask = f"../hk/sw{self.__triggerID}bdqcb.hk.gz",
            tstart = startTime,
            tstop = stopTime,
            clobber = "YES",
            outunits = "COUNTS"
        )
        print(output.stdout)






if __name__ == "__main__":
    GRBName: str = "GRB080319B"
    data: ProcessSwiftData = ProcessSwiftData(GRBName)


