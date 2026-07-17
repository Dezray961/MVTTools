from astropy.io import fits


filePath = "data/reproc/00306757000/bat/event/outputBurst.lc"

data = fits.open(filePath)

print(data.info())

print(type(data[0]))
