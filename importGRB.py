import pandas as pd
import numpy as np
from scipy.signal import lfilter


# class to hold the data for a GRB
class GRBData:
    """GRBData class. This class holds the data for a GRB, including the name of the GRB, the path to the CSV file containing the data, and the data itself as a pandas DataFrame.

        Args:
            name (str): Name of the GRB, this should follow the standard naming convention for GRBs, e.g. "GRB080319B"
            CSVfilePath (str): Path to the CSV file containing the GRB data. The CSV file should have the following columns: time, rate, error, and any other columns that may be relevant for the analysis. The time column should be in seconds and can either be relative to the trigger time of the GRB or absolute. The rate column should be in counts per second. The error column should be the error on the rate measurement in counts per second.
        
        ## Notes:
            The csvToDataFrame() method reads in the CSV file and stores the data in a pandas DataFrame. Pandas is a suboptimal choice for this task as it can be computationally expensive. However, it is quick to convert columns into numeric data. Try not to do any operations on the DataFrame itself; rather convert the columns to numpy arrays and operate on those instead, then put the results back into the DataFrame if needed.
    """
    def __init__(
            self,
            name: str,
            CSVfilePath: str
            ) -> None:
        self.name: str = name
        self.CSVfilePath: str = CSVfilePath
        self.csvToDataFrame()
        self.correctTime()
        self.timeInBin()
        self.logRate()

    # method to read in a CSV file and store the data
    def csvToDataFrame(
            self
            ) -> None:
        # This method reads in the CSV file and stores the data in a pandas DataFrame. Pandas is
        # a suboptimal choice for this task as it can be computationally expensive. However, it is
        # quick to convert columns into numeric data. Try not to do any operations on the DataFrame
        # itself; rather convert the columns to numpy arrays and operate on those instead, then put
        # the results back into the DataFrame if needed.
        self.data: pd.DataFrame = pd.read_csv(self.CSVfilePath)
        # set the index to be the first column
        self.data.set_index(self.data.columns[0], inplace=True)
        #remove the first column name
        self.data.index.name = None
        # convert all the data to numeric, coercing errors to NaN

        for column in self.data.columns:
            self.data[column] = pd.to_numeric(self.data[column], errors='coerce')


    # method to correct the time array to be relative to the trigger time
    def correctTime(
            self
            ) -> None:
        triggerTime: float = self.data['time'].to_numpy()[0]
        self.data['time'] = self.data['time'] - triggerTime


    # method to get the time in each bin
    def timeInBin(
            self
            ) -> None:
        timeInBin: np.ndarray = np.zeros_like(self.data['time'])
        for i in range(len(self.data['time'])):
            if i == len(self.data['time']) - 1:
                timeInBin[i] = 0
            else:
                timeInBin[i] = self.data['time'].to_numpy()[i + 1] - self.data['time'].to_numpy()[i]
        self.data['timeInBin'] = timeInBin
    

    # method to get the log of the rate
    def logRate(
            self
            ) -> None:
        rate: np.ndarray = self.data['rate'].to_numpy()
        logRate: np.ndarray = np.log(rate)
        self.data['logRate'] = logRate


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    grbName: str = "GRB080319B"
    csvFilePath: str = f"data/processed/{grbName}LC.csv"
    grbData: GRBData = GRBData(grbName, csvFilePath)
    print(grbData.data.head())
    fig, ax = plt.subplots()
    ax.plot(grbData.data['time'], grbData.data['rate'], color='silver', label='Lightcurve')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Count Rate (counts/s)')
    ax.set_xlim(-0.5, 62)
    ax.set_ylim(0, max(grbData.data['rate'].where(grbData.data['time'] < 60).dropna()) * 1.05)
    ax.legend()