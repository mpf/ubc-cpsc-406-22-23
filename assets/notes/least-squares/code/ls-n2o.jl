# This file was generated, do not modify it. # hide
fname = "mlo_N2O_All.dat"
if !isfile(fname)
	download("https://gml.noaa.gov/aftp/data/hats/n2o/insituGCs/CATS/hourly/mlo_N2O_All.dat", fname)
end