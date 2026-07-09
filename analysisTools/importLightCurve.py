from pandas import read_csv, to_numeric, DataFrame
from numpy import ndarray, zeros_like, diff, argmax


# class to hold the data for a GRB
class LightCurveData:
    """GRBData class. This class holds the data for a GRB, including the name of the GRB, the path to the directory containing the CSV files, and the data itself as a pandas DataFrames.

        Args:
            GRBName (str): Name of the GRB, this should follow the standard naming convention for GRBs, e.g. "GRB080319B"
        
        Attributes:
            name (str): Name of the GRB, this should follow the standard naming convention for GRBs, e.g. "GRB080319B"
            preBurstData (DataFrame): DataFrame containing the pre-burst light curve data for the GRB. This data is read in from the CSV file "data/processed/{GRBName}/{GRBName}PreBurstLC.csv"
            burstData (DataFrame): DataFrame containing the burst light curve data for the GRB. This data is read in from the CSV file "data/processed/{GRBName}/{GRBName}LC.csv"
        
        ## Notes:
            The csvToDataFrame() method reads in the CSV file and stores the data in a pandas DataFrame. Pandas is a suboptimal choice for this task as it can be computationally expensive. However, it is quick to convert columns into numeric data. Try not to do any operations on the DataFrame itself; rather convert the columns to numpy arrays and operate on those instead, then put the results back into the DataFrame if needed.
    """
    def __init__(
            self,
            GRBName: str,
            source: str
        ) -> None:
        """Initializes the GRBData class. This method reads in the CSV files for the given GRB name and stores the data in pandas DataFrames.

        Args:
            GRBName (str): Name of the GRB, this should follow the standard naming convention for GRBs, e.g. "GRB080319B"
            source (str): The source of the data e.g. "Swift", "Fermi", "SVOM".
        """
        self.name: str = GRBName
        self.source: str = source
        self.__preBurstCSV: str = f"data/processed/{self.name}/{self.name}PreBurstLC.csv"
        self.__burstCSV: str = f"data/processed/{self.name}/{self.name}BurstLC.csv"
        print("Importing the pre-burst light curve data...")
        self.preBurstData: DataFrame = self.__generateDataFrame(self.__preBurstCSV)
        print("Importing the burst light curve data...")
        self.burstData: DataFrame = self.__generateDataFrame(self.__burstCSV)
        self.__truncateAtSlewPoint()
        print("Data import complete.")

    def __generateDataFrame(
            self,
            CSVFilePath: str
        ) -> DataFrame:
        # read in the CSV file and store the data as a pandas DataFrame
        data: DataFrame = read_csv(CSVFilePath)

        # convert all the data to numeric, coercing errors to NaN
        for column in data.columns:
            data[column] = to_numeric(data[column], errors='coerce')

        # correct the time column to be relative to the trigger time
        triggerTime: float = data['time'].to_numpy()[0]
        data['time'] = data['time'] - triggerTime

        # get the time in each bin
        timeArray = data['time'].to_numpy()
        timeInBin: ndarray = zeros_like(timeArray)
        timeInBin[:-1] = diff(timeArray)
        data['timeInBin'] = timeInBin

        return data


    def __truncateAtSlewPoint(
            self
        ) -> None:
        """Truncates the light curve data at the slew point. While slewing the BAT data is unreliable and is not included in the light curve data. This method finds the slew point in the light curve data and truncates the data at that point."""
        # find the slew point in the burst data
        maxBinSize: float = self.burstData['timeInBin'].max()
        if maxBinSize > 100e-6:
            slewPointIndex: int = argmax(self.burstData['timeInBin'].to_numpy())
            # truncate the data after the slew point
            self.burstData = self.burstData.iloc[:slewPointIndex]
            self.preBurstData = self.preBurstData.iloc[:slewPointIndex]


if __name__ == "__main__":
    # print the working directory
    import matplotlib.pyplot as plt
    grbName: str = "GRB080319B"
    data: LightCurveData = LightCurveData(grbName, "Swift")

    # plot the data
    fig, axs = plt.subplots(2, 1)
    axs[0].plot(
        data.burstData['time'],
        data.burstData['rate'],
        color='blue',
        label='Lightcurve'
        )
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Count Rate (counts/s)')
    axs[0].legend()
    axs[1].plot(
        data.preBurstData['time'],
        data.preBurstData['rate'],
        color='silver',
        label='Pre-Burst Lightcurve'
        )
    axs[1].set_xlabel('Time (s)')
    axs[1].set_ylabel('Count Rate (counts/s)')
    axs[1].legend()
    fig.show()
