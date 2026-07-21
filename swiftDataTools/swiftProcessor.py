from os import system, environ, getcwd, makedirs
from os.path import exists, join, abspath
from heasoftpy import Config, bateconvert, batmaskwtevt, batbinevt, batupdatephakw, batphasyserr, batdrmgen, HSPTask
from contextlib import chdir
from swiftDataTools.swiftBATCatalogueGRB import getObservationID, importData, getCoordinates, getStartStopTime
from pathlib import Path
from astropy.io import fits
import numpy as np
import subprocess

class ProcessSwiftData:
    """Class to process Swift BAT data for a given GRB. This class handles the processing of Swift BAT data for a given GRB, including checking and applying gain correction, checking and applying mask weighting, extracting light curves for different time periods, and converting the light curve data into CSV files. It has the ability to process a custom time range, however this will require the user to call methods from outside the class."""
    def __init__(
            self,
            GRBName: str
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
            # setup local separate writeable parameter directory
            local_pfiles = abspath("./pfiles")
            makedirs(local_pfiles, exist_ok=True)
            environ["PFILES"] = f"{local_pfiles};{self.__headasPath}/syspfiles"
            
            # inject HEASoft binaries directly into Python's active execution path
            environ["PATH"] = f"{join(self.__headasPath, 'bin')}:{environ.get('PATH', '')}"
            environ["LD_LIBRARY_PATH"] = f"{join(self.__headasPath, 'lib')}:{environ.get('LD_LIBRARY_PATH', '')}"
        else:
            raise EnvironmentError("Error: HEADAS environment variable not found. Did you initialize HEASoft?")
        Config.allow_failure = False

        # set the CALDB environment variables
        self.__caldbPath = environ.get("CALDB")
        if self.__caldbPath:
            environ["CALDBCONFIG"] = join(self.__caldbPath, 'software', 'tools', 'caldb.config')
            environ["CALDBALIAS"] = join(self.__caldbPath, 'software', 'tools', 'alias_config.fits')

        self.GRBName = GRBName
        self.energyBins = config.preProcessingConfig.swiftBATConfig.download.energyRange
        self.__data, self.__columnNames = importData("swiftDataTools/summary_general.csv")
        self.__triggerID: str = getObservationID(
            GRBName,
            self.__data,
            isTrigID=False
            )

        # find the right ascension and declination of the GRB from the summary_general.csv file.
        self.__rightAscension, self.__declination = getCoordinates(self.GRBName, self.__data)

        # create the output directory if it does not exist
        Path(f"data/processed/{self.GRBName}").mkdir(parents=True, exist_ok=True)

        # process the data using the HEASoft tools
        with chdir(f"data/reproc/{self.__triggerID}/bat/event"):
            print(f"Processing data for {GRBName}...")
            self.__eventFilename: str = f"sw{self.__triggerID}bevshsp_uf.evt.gz"
            
            # read the event FITS file using astropy.io.fits. 
            self.__eventFITS: fits.HDUList = fits.open(self.__eventFilename)

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

            # get the pre and post bust midpoints 
            self.__preBurstMidpoint: float = (
                self.__startTime - self.__burstDuration
                )
            self.__postBurstMidpoint: float = (
                self.__stopTime + self.__burstDuration
                )


            # extract the light curves
            for period in range(3):
                self.__extractLightCurve(period)

            # TODO
            # load the light curve data from the output files using astropy.io.fits
            self.__preBurstLightCurve: fits.HDUList = fits.open("outputPreBurst.lc")
            self.__burstLightCurve: fits.HDUList = fits.open("outputBurst.lc")
            self.__postBurstLightCurve: fits.HDUList = fits.open("outputPostBurst.lc")

            # generate the spectrum for the burst period
            self.__generateSpectrum()

            # generate the response matrix using batdrmgen
            self.__generateResponseMatrix()

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
        # get the GAINAPP and GAINMETH keywords from the event file header
        GAINAPP: bool = self.__eventFITS[1].header['GAINAPP']
        GAINMETH: str = self.__eventFITS[1].header['GAINMETH']

        if GAINAPP and GAINMETH == "FIXEDDAC":
            print("Gain correction already applied")
            return True
        else:
            print("Gain correction has not been applied.")
            return False


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
        if not exists(self.__calibrationFile):
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
            # get the right ascension and declination from the event file header
            rightAscension: float = self.__eventFITS[1].header['BAT_RA']
            declination: float = self.__eventFITS[1].header['BAT_DEC']

            # load the absolute tolerance for position matching from the config file
            absoluteTolerance: float = float(config.preProcessingConfig.swiftBATConfig.processing.positionTolerance)

            # check if the right ascension and declination are within the absolute tolerance of the expected values
            rightAscensionMatch: bool = np.isclose(
                rightAscension,
                self.__rightAscension,
                atol=absoluteTolerance
                )
            declinationMatch: bool = np.isclose(
                declination,
                self.__declination,
                atol=absoluteTolerance
                )
            
            if rightAscensionMatch and declinationMatch:
                print("Right ascension and declination match expected values, mask weighting already applied")
                return True
            else:
                print("Right ascension and declination do not match expected values, mask weighting has not been applied")
                return False


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
            # get the BATCREAT keyword from the event file header
            BATCREAT: str = self.__eventFITS[1].header['BATCREAT']
            version: float = float(BATCREAT.split(' ')[1].strip())
            if version >= 1.16:
                print(f"batmaskwtevt version is {version} >= 1.16")
                return True
            else:
                print(f"batmaskwtevt version is {version} < 1.16")
                return False


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
        if not exists(self.__attitudeFile):
            ### download the attitude file - impliment this later
            raise FileNotFoundError(f"Attitude file not found: {self.__attitudeFile}")
            
        # find the quality map file in the hk directory.
        self.__qualityMapFile: str = f"../hk/{self.__eventFilename.replace('bevshsp_uf.evt.gz', 'bdqcb.hk.gz')}"
        if not exists(self.__qualityMapFile):
            ### download the quality map file - impliment this later
            raise FileNotFoundError(f"Quality map file not found: {self.__qualityMapFile}")

        # unzip the file if it is compressed. Needed for the batmaskwtevt command to work. 
        if self.__eventFilename.endswith(".gz"):
            self.__unzipFile(self.__eventFilename)
            self.__eventFilename = self.__eventFilename.replace(".gz", "")

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


    def __extractLightCurve(
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
                startTime: float = self.__preBurstMidpoint - 1.0
                stopTime: float = self.__preBurstMidpoint + 1.0
                print("Extracting pre-burst uniform light curve")
            case 1: # burst
                fileName: str = "outputBurst.lc"
                if config.preProcessingConfig.swiftBATConfig.processing.fullBurst:
                    startTime: float = self.__startTime - 2.0
                    stopTime: float = self.__stopTime + 2.0
                    print("Extracting full burst uniform light curve")
                else:
                    startTime: float = self.__startTime - 2.0
                    stopTime: float = (
                        self.__startTime
                        + config.preProcessingConfig.swiftBATConfig.processing.sliceDuration - 2.0
                        )
                    print("Extracting burst uniform light curve")
            case 2: # post-burst
                fileName: str = "outputPostBurst.lc"
                startTime: float = self.__postBurstMidpoint - 1.0
                stopTime: float = self.__postBurstMidpoint + 1.0
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
            timedel = config.preProcessingConfig.swiftBATConfig.processing.initialBinSize,
            timebinalg = "u",
            energybins = self.energyBins,
            detmask = f"../hk/sw{self.__triggerID}bdqcb.hk.gz",
            tstart = startTime,
            tstop = stopTime,
            clobber = "YES",
            outunits = "COUNTS"
        )
        print(output.stdout)


    def __subprocessRunCommand(
            self,
            commandString: str
        )-> None:
        """Runs a command string in a bash shell using subprocess.run. This is a wrapper function for subprocess.run, which is used to run commands in a bash shell. The function takes in a command string, and uses subprocess.run to run the command in a bash shell. The function does not return anything, but it will print out any output from the command to the console.

        Args:
            commandString (str): The command string to be run. This should be a valid bash command string.
        """
        output = subprocess.run(
            commandString,
            shell=True,
            executable="/bin/bash",
            capture_output=True,
            text=True
        )
        print(f"STDOUT Log:\n{output.stdout}")
        if output.returncode != 0:
            print(f"Command failed with code {output.returncode}!")
            print(f"Error Log:\n{output.stderr}")


    def __generateSpectrum(
            self
        )-> None:
        print("Generating spectrum for burst period")
        # run batbinevt to create the spectrum for the burst period
        output = batbinevt(
            infile = self.__eventFilename,
            outfile = "outputSpectrum.pha",
            outtype = "PHA",
            timedel = 0.0,
            timebinalg = "u",
            tstart = self.__startTime,
            tstop = self.__stopTime,
            energybins = 'CALDB:80',
            outunits = "RATE",
            detmask = f"../hk/sw{self.__triggerID}bdqcb.hk.gz",
            clobber = "YES"
        )
        print(output.stdout)

        print("Applying corrections to the spectrum")
        # run batupdatephakw and batphasyserr to apply corrections to the spectrum
        # currently this is done using subprocess.run to run the commands in a bash shell. 
        # This is because the heasoftpy wrapper for batupdatephakw is throwing errors about the 
        # CALDB environment variable not being set, even though it is set??? TODO
        commandString = (
            f"source {self.__headasPath}/headas-init.sh && "
            f"batupdatephakw outputSpectrum.pha sw{self.__triggerID}bevtr.fits.gz clobber=YES"
        )
        self.__subprocessRunCommand(commandString)  

        # Construct the bash command string, injecting your HEASoft and CALDB setup
        commandString = (
            f"source {self.__headasPath}/headas-init.sh && "
            f"batphasyserr outputSpectrum.pha CALDB clobber=YES"
        )
        self.__subprocessRunCommand(commandString)  


    def __generateResponseMatrix(
            self
        )-> None:
        print("Generating response matrix for burst period")
        # run batdrmgen to create the response matrix for the burst period
        output = batdrmgen(
            infile = "outputSpectrum.pha",
            outfile = "outputResponse.rsp",
            hkfile = 'NONE',
            clobber = "YES"
        )
        print(output.stdout)


if __name__ == "__main__":
    from loadConfig import importConfiguration
    config = importConfiguration("config.yaml")
    GRBName: str = "GRB080319B"
    data: ProcessSwiftData = ProcessSwiftData(GRBName)


