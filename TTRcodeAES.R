#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# This script provides a portable demonstration of how we fitted the plant growth
# models in this study. It assumes that the user is already a proficient user of R
# An example data set sd_Combretum apiculatum.Rdata is provided, this includes
# forcing data from CHELSA 2.1 and occurrence data as described in the manuscript.
#
# Disclaimer: The script is to be interpreted as documentation of how we did 
# the analyses and should not be interpreted as software that will produce 
# coherent error messages should you do something wrong.
#
# For this to work you require the following files in your R session's 
# working directory:
# TTRcodeAES.R (this file)
# include.c
# include.R
# sd_Combretum apiculatum.Rdata
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


library(DEoptim) # to fit the model
library(ROCR)    # to evaluate the model

#------------------------------------------------------------------------------#
# Source a file that includes all helper functions and will compile the library.
# Requires an Rtools installation!
#------------------------------------------------------------------------------#
source("include.R")

#------------------------------------------------------------------------------#
# prepare the input data
#------------------------------------------------------------------------------#

#-- Assume you have a data frame with appropriate data, e.g.
load("sd_Combretum apiculatum.Rdata") 
head(sd) # take a peak at what the data should look like
#-- See comments associated with prodcdata.chelsa for more information on data
#-- formatting. The function make.sd.chelsa documents how to make
#-- the sd data files, but it is not portable: to use make.sd.chelsa you would
#-- need to download the appropriate input data files.  See the methods section
#-- of the manuscript for information on the data sources.

#-- Create a list with the input data and add the photosynthesis forcing
photo = "C3"
if(photo=="C3") c.data.A<-prodcdata.chelsa(sd,C3=TRUE,CA=338/10) # 338ppm = CO2 concentration 1980
if(photo=="C4") c.data.A<-prodcdata.chelsa(sd,C3=FALSE,CA=338/10) # 338ppm = CO2 concentration 1980

#-- Create another list that structures the data for DEoptim
MyData <- make.MyData(c.data.A)

#------------------------------------------------------------------------------#
# estimate the parameters using DEoptim
# this will take a few minutes
# read the DEoptim help for more information
#------------------------------------------------------------------------------#
deoptim_pop = 50
deoptim_iter = 250
deoptim_trace = 5 # if you want to monitor progress, else FALSE

SolDEoptim<-DEoptim(ModelDE, upper=MyData$up, lower=MyData$lo,             
            control=DEoptim.control(NP = deoptim_pop,
                                    itermax = deoptim_iter,
                                    trace = deoptim_trace,
                                    CR=0.9, F=0.8, c=0.750,
                                    initialpop = NULL,
                                    strategy = 2,
                                    parallelType=0),
            Data=MyData)

#------------------------------------------------------------------------------#            
# evaluate the model using a wrapper function to ROCR
# by default the get_eval function writes a ROC plot to the working directory
#------------------------------------------------------------------------------#
eval <- get_eval(SolDEoptim$optim$bestmem,MyData)
print(eval)
