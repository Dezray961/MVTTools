from pandas import DataFrame

def getWindowIndices(
        data: DataFrame,
        timeWindowSize: float
        ) -> list:
    """This function returns the indices of the light curve data that fall within a specific time window size.

    Args:
        data (DataFrame): the denoised light curve data
        timeWindowSize (float): the size of the time window to consider

    Returns:
        list: a list of indices corresponding to the time values that fall within the specified time window size
    """
    # get the timeInBins column from the data
    timeInBins = data['timeInBin'].to_numpy()

    # initialize variables to keep track of the window indices and the accumulated time
    windowIndices: list[tuple[int, int]] = []
    accumulatedTime: float = 0.0
    startIndex: int = 0
    endIndex: int = 0

    # get the number of time bins
    n = len(timeInBins)

    # iterate through the time bins and accumulate the time until it exceeds the time window size
    while startIndex < n:
        while accumulatedTime < timeWindowSize and endIndex < n:
            accumulatedTime += timeInBins[endIndex]
            endIndex += 1
        # add the window indices to the list and reset the accumulated time and start index
        if accumulatedTime >= timeWindowSize:
            windowIndices.append((startIndex, endIndex - 1))
            accumulatedTime -= timeInBins[startIndex]
            startIndex += 1
        else:
            break
    
    return windowIndices