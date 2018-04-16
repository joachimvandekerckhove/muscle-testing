####  Helper function for the muscle testing project
# A collection of helper functions for the muscle testing paper.
# (c)2018 CC-BY joachim vandekerckhove

### Function: kendalls.tau()
# Computes an unnormalized Kendall's tau-a coefficient
kendalls.tau <- function(list1, list2) {
# If there are any NAs, return NA
if (any(is.na(list1+list2))) {
  return(NA)
}

# If no NA entries, compute Kendall's tau
n <- 1:length(list2)

rank1 <- sort(x = list1, index.return = TRUE)$ix
rank2 <- sort(x = list2, index.return = TRUE)$ix

ii <- outer(n, n, function(x, y){x})
jj <- outer(n, n, function(x, y){y})

sign1 <- rank1[jj[jj > ii]] > rank1[ii[jj > ii]]
sign2 <- rank2[jj[jj > ii]] > rank2[ii[jj > ii]]

return(sum(sign1 != sign2))
}

### Function: expected.tau()
# Computes tau for every possible outcome and tabulates the frequencies
expected.tau <-function(concordance, length) {
  # All the possible response sequences assuming some level of agreement
  domain <- combinat::permn(length-concordance)
  N <- length(domain)
  
  # Compute tau for all pairs
  arr <- sapply(domain, kendalls.tau, domain[[1]])
  
  # Tabulate the taus
  freqs  <- table(arr)
  labels <- as.numeric(labels(freqs)[[1]])
  values <- as.numeric(freqs)
  
  # Convert to histogram
  out <- array(0, length * (length-1)/2 + 1)
  out[labels+1] <- values
  
  # Convert to probability function and return
  return(out / sum(out))
}

### Function: tau.matrix()
# Compute probability function for each possible level of concordance
tau.matrix <- function(length) {
  return(t(sapply(0:(length-1),expected.tau,length)))
}

### Function: norm()
# A normalization function
norm <- function(x) { 
  return(x / sum(x))
}

### Function: read.mt.data()
# This function will read the .csv data file from OSF and format into a list
read.mt.data <- function(data.url = "https://osf.io/4rp7c/download") {
  # Read .csv from the url and add column names
  a <- read.csv(file   = url(data.url), 
                header = FALSE,
                col.names = c("session","tester","participant",
                              "suppl.1","suppl.2","suppl.3","suppl.4","suppl.5",
                              "label.1","label.2","label.3","label.4","label.5",
                              "confidence","belief"))
  
  # Collate the "label" and "supplements" columns into matrices
  a$labels      <- cbind(a$label.1, a$label.2, a$label.3, a$label.4, a$label.5)
  a$supplements <- cbind(a$suppl.1, a$suppl.2, a$suppl.3, a$suppl.4, a$suppl.5)
  
  # Trim the redundant columns and return
  return(a[,-grep("\\.", colnames(a))])
}

### Function: process.data()
# Extract the relevant subset of the muscle testing data
process.data <- function(muscle.data,
                         tester.subset  = c(1:5, 7, 9),
                         data.type      = 'supplements',
                         confident.only = FALSE, 
                         believers.only = FALSE) {
  # Use the confidence and belief columns to select all rankings or a subset of rankings, based on inputs
  use <- switch (confident.only + 2 * believers.only + 1, 
                 !logical(length(muscle.data$belief)),         # all rankings
                 muscle.data$confidence==1,                    # from confident testers only
                 muscle.data$belief==1,                        # from participants who believe only
                 muscle.data$confidence & muscle.data$belief,  # from confident testers and participants who believe only 
                 stop("process.data# Invalid input(s) check 'confident_only' and 'believers_only' inputs.")
  )
  
  # Choose which rankings to analyze based on datatype input
  rankings <- switch(data.type, 
                     supplements = muscle.data$supplements[use,], 
                     labels      = muscle.data$labels[use,], 
                     stop("process.data# Invalid input. Check 'data.type' input.")
  )
  
  ## Subset the data
  # session      <- muscle.data$session[use] #(currently unused)
  tester       <- muscle.data$tester[use]
  participant  <- muscle.data$participant[use]
  
  # The list of unique participant numbers
  participantList <- unique(muscle.data$participant[use])
  
  ## Compute all pairwise comparisons to generate the output vector (tau)
  # Initialize the output variable
  taus <- integer()
  
  # Loop over unique participants
  for (p in 1:length(participantList)) {
    # Find the rows containing rankdata for the current participant (p)
    # for all the muscle testers (listed in the input variable mt) this 
    # gives a set of rankings
    foundranks <- rankings[muscle.data$participant[use]==participantList[p] & 
                             is.element(muscle.data$tester[use],tester.subset),]
    
    # Number of rankings found
    M <- dim(foundranks)[1]
    
    # Preallocate a MxM variable here
    taus.tmp <- array(NA, c(M,M))
    
    # Loop over all M rankings
    for (m1 in 1:M) {
      for (m2 in 1:M) {
        if (m2 > m1) {
          # Calculate tau between these pairs of rankings
          taus.tmp[m1,m2] <- kendalls.tau(foundranks[m1,], foundranks[m2,])
        }
        # (Note: If this combination is not reached (e.g., (2,1)), this 
        # remains NA. It is also NA if any of the rankings are NA.)
      }
    }
    
    # Find all the taus that are not NA and add them to the output variable
    taus <- c(taus, taus.tmp[!is.na(taus.tmp)])
  }
  return(taus)
}

### Function: run.mt()
# This function runs the analysis

run.mt <- function(muscle.data,
                   tester.subset  = c(1:5, 7, 9),
                   data.type      = 'supplements',
                   confident.only = TRUE, 
                   believers.only = TRUE,
                   sigma.prior    = function(x){ seq(1,1,along.with = x) },
                   lambda.prior   = function(x){ seq(1,1,along.with = x) }){
  # First, process the data
  data                <-  process.data(muscle.data    = muscle.data,
                                       tester.subset  = tester.subset,
                                       data.type      = data.type,
                                       confident.only = confident.only, 
                                       believers.only = believers.only)
  
  # The number of items to be sorted: 
  items               <-  5
  
  # Possible values of agreement (sigma), based on the number of items:
  sigma               <-  0:(items-1)
  
  # A list of possible lapse rates (lambda):
  lambda              <-  seq(0, 1/2, by = 1/length(data))
  ln                  <-  length(lambda)
  
  # Copy the lookup table for each lambda:
  sigma.table         <-  tau.matrix(length = items)
  
  # Make one matrix for the lapse lookup table:
  lapse.table         <-  t(
    matrix(
      rep(expected.tau(concordance = 0, length = items), each = items), 
      ncol = items,
      byrow = TRUE
    )
  )
  
  # Loop over lambda values to compute mixture probability:
  likelihood          <-  array(NA, c(items, items*(items-1)/2+1, ln))
  for (ctr in 1:ln) {
    likelihood[,,ctr] <-  lambda[ctr] * lapse.table + (1 - lambda[ctr]) * sigma.table
  }
  
  # Compute priors
  prior.sigma         <-  norm(sigma.prior(sigma))
  prior.lambda        <-  norm(lambda.prior(lambda))
  prior.joint         <-  outer(prior.sigma, prior.lambda)
  
  ## Compute inferential quantities
  # Apply lookup table to obtain summed log-likelihood:
  sum.log.likelihood  <-  apply(log(likelihood)[,data + 1,], c(1,3), sum)
  
  # Compute joint log-posterior:
  log.posterior       <-  sum.log.likelihood + log(prior.joint)
  
  # Get joint and marginal posteriors
  posterior.joint     <-  norm(exp(log.posterior - median(log.posterior)))
  posterior.lambda    <-  apply(posterior.joint, 2, sum)
  posterior.sigma     <-  apply(posterior.joint, 1, sum)
  
  ## Some quantities for decision-making
  posterior.ratio     <-  posterior.sigma[1] / sum(posterior.sigma[-1])
  prior.ratio         <-  prior.sigma[1] / sum(prior.sigma[-1])
  
  bayes.factor        <-  posterior.ratio / prior.ratio
  
  return(
    list(
      joint.posterior   =  posterior.joint  ,
      lambda.posterior  =  posterior.lambda ,
      sigma.posterior   =  posterior.sigma  ,
      bayes.factor      =  bayes.factor     )
  )
}
