library(neonUtilities)
outdir <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw"
outpath <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw/rawTempStacked"
# Download soil temperature data
zipsByProduct(dpID = "DP1.00041.001", site = "all", 
             # enddate = "2015-12", 
              avg = 30, check.size = T, 
              savepath = outdir)
out <- stackByTable(filepath = paste0(outdir, "/filesToStack00041"), savepath = outpath, folder = TRUE, 
                    nCores = 16, saveUnzippedFiles = FALSE)

  