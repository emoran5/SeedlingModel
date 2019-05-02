# SeedlingModel
Code for seedling model described in Moran et. al 2019, Ecosphere

Raw data files are:
- PlotInfo2.txt (Names, elevation, and location info for each site - UTM coordinates rounded to nearest 500)
- tree2.csv (Tree size and mortality data)
- treeyears.txt (Which years trees in each site were measured)
- quadrat_precise.csv (locations of seedling quadrats within site)
- TaggedSdl.txt (Tagged seedling observations)

DataOrganization.r is the R file that processes this raw data into the input file TagSdl_Data5.RData

G14a2.r is the R file for the best-fit growth model
S15m1.r is the R file for the best-fit survival model
