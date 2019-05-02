# SeedlingModel
Code for seedling model described in Moran et. al 2019, Ecosphere

Raw data files are:
- PlotInfo2.txt (Names, elevation, and location info for each site - UTM coordinates rounded to nearest 500)
- tree2.csv (Tree size and mortality data)
- treeyears.txt (Which years trees in each site were measured)
- quadrat_precise.csv (locations of seedling quadrats within site)
- TaggedSdl.txt (Tagged seedling observations)
- climate_data_long.txt (climate data for all years for all sites from CA BCM model)

DataOrganization.r is the R file that processes this raw data into the input file TagSdl_Data5.RData

G14a2.r is the R file for the best-fit growth model
S15m1.r is the R file for the best-fit survival model

FutureProjTrans.r (R code that creates the transition rates under future climate figure)
FutureProjMort.r (R code that creates the survival rates under future climate figure)
