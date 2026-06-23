"""
for a given data set, this code will find the MVT using the following process:
1. Take a class object with the data set as input. This should have the following attributes:
    - time: array of time values
    - deltaTime: array of time intervals
    - rate: array of rate values
    - deltaRate: array of rate uncertainties
2. Denoise the data using Haar wavelet transform
3. Bin the data into constant S/N bins
4. Find the Haar coefficients for the data using the un-decimated Haar wavelet transform
5. Find the MVT for the data using a statistical analysis of the Haar coefficients
    This is will be a 2σ confidence level deviation from a straight line fit to the start of the Haar coefficients. The MVT is the time bin at which this deviation occurs.
6. Store the MVT in the class object
"""


