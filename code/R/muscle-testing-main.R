## CONSISTENCY OF MUSCLE TEST RESULTS: DATA ANALYSIS
# 
# This script performs all analyses, in particular Bayesian inference on 
# rank data using Kendall's tau distance.  
# 
# For details about the data-coding scheme, see the README file
# (CMTR_readme.txt).  This file and others may be found at our ... 
#               OSF project page:      https://osf.io/8d4wy/
#               OSF pre-registration:  https://osf.io/wymjr/

source('muscle-testing-helpers.R')

data <- read.mt.data("https://osf.io/4rp7c/download")

print((a <- run.mt(muscle.data = data,
                   confident.only = TRUE, 
                   believers.only = TRUE))$bayes.factor)

print((a <- run.mt(muscle.data = data,
                   confident.only = FALSE, 
                   believers.only = TRUE))$bayes.factor)

print((a <- run.mt(muscle.data = data,
                   confident.only = TRUE, 
                   believers.only = FALSE))$bayes.factor)

print((a <- run.mt(muscle.data = data,
                   confident.only = FALSE, 
                   believers.only = FALSE)$bayes.factor))

