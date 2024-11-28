# This R script implements a basic greedy model for inferring
# fatal incidence rates from Covid deaths in English hospitals using a 
# modified Pearson statistic as the loss function. Where the loss is between
# real death data and simulated death data. The simulated death data is then 
# used to infer infection dates. The algorithm fits by randomly shifting individual 
# deaths by a random amount and checking if the new infection/death days reduce the loss.

# Finally, we check our results by performing an uncertainty analysis with
# error quantiles created by bootstrapping the real dataset with fake death data.
# The fake death data is sampled from a poisson distr for each day with rate 
# parameter = true death data for that day.


# Load data from the "engcov.txt" file, assuming the first row contains column headers
covid_data <- read.table("engcov.txt", 
                         header = TRUE,   # The file has column headers
                         nrow = 150)      # Read the first 150 rows


modified_pearson <- function(real_death_days, simulated_death_days) {
  # Description: This function calculates the modified Pearson loss between 
  #              two vectors. The loss is computed by taking the squared difference 
  #.              then, normalised by the second vector or 1.
  #              This function ensures no division by zero by replacing any zero 
  #              death counts in the denominator with 1.
  #              
  # Inputs:
  #   - real_death_days: Vector of non-negative integers.
  #   - simulated_death_days: Vector of non-negative integers.
  # Outputs:
  #   - Sum of normalized squared differences between real and simulated deaths
  #     (i.e., the modified Pearson loss)
  #
  difference_squared <- (real_death_days - simulated_death_days)^2
  denominator <- pmax(1, simulated_death_days)  # Avoid division by zero
  return(sum(difference_squared / denominator))
}


generate_sequence <- function(step_count) {
  # Description: This function generates a sequence of shift values for deconv
  #               based on the input step count. It returns smaller shift ranges
  #               as the step count increases
  # Inputs:
  #   - step_count: Integer representing the number of steps
  # Outputs:
  #   - Vector of integers with step values based on the step count
  if (step_count < 25) {
    return(c(-8, -4, -2, -1, 1, 2, 4, 8))
  } else if (step_count < 50) {
    return(c(-4, -2, -1, 1, 2, 4))
  } else {
    return(c(-2, -1, 1, 2))
  }
}


deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL) {
  # Function: deconv
  # Description: This function generates predicted infection days given data
  #               on deaths per day for a virus. It minimises the loss function
  #               modified_pearson using a greedy algorithm that randomly
  #               shifts infection/death days to see if they reduce the loss.
  #.              The function includes a bootstrapping mode with the bs 
  #.              indicator variable that runs the algorithm on "simulated"
  #.              true death data based on the real data.
  #              
  # Inputs:
  #   - t: Vector of positive integers representing potential days people can die on
  #   - deaths: Vector of non-negative integers representing deaths per day in t
  #   - n.rep: integer, number of epochs to run the model fit for
  #   - bs: Boolean: bootstrapping mode indicator variable
  #         bs=TRUE: fit model to fake death data generated from true death data (t/deaths)
  #         bs=FALSE: fit model to true death data (t/deaths)
  #   - t0: vector or NULL: variable for parsing pre-defined infection days
  #         t0=NULL: generate vector of infection days 
  #         t0=vector: use pre-defined infection days for model fit
  # Outputs:
  #   - list(P, inft, t0),
  #       - P: vector of Pearson loss values after each epoch
  #       - inft: Matrix(310,n.rep), containing infection days for each iteration
  #               of the algorithm
  #       - t0: Vector of final infection days after n.rep epochs of the algorithm
  
  
  # Set up probability distribution for infection-to-death duration
  days_survived <- 1:80                      # Range of survival days
  meanlog <- 3.152                           # Mean of log-normal distribution
  sdlog <- 0.451                             # Standard deviation of log-normal distribution
  probs <- dlnorm(days_survived, meanlog, sdlog)
  normalised_probs <- probs / sum(probs)     # Normalise probabilities to sum to 1
  
  # Set up initial death days based on bs (bootstrapping) parameter
  if (bs == FALSE) {
    # Without bootstrapping: replicate t according to deaths counts
    death_day <- rep(t, times = deaths)
  } else {
    # Perform the bootstrapping with poisson data
    
    
    # pre define the matrix and vector for storing each bootstrap result
    bs_data <-  matrix(0, nrow = 310, ncol = n.rep)
    P <- numeric(n.rep)     
    for (epoch in 1:n.rep){
      if(epoch %% 10 ==0 ){
        print(paste("Bootstrap: ",epoch,"/",n.rep))
      }
      # Bootstrapping: sample new death data using a Poisson distribution
      death_sample <- rpois(length(deaths), lambda = deaths)
      death_day <- rep(t, times = death_sample)
      # fit the model using one step
      poisson_fit <- deconv(1:310,tabulate(death_day, nbins = 310), n.rep = 1, bs = FALSE, t0 = t0)
      # record the resulting t0/p values in each column/ element 
      # of the pre-defined matrix/vector
      bs_data[, epoch] <- tabulate(poisson_fit$t0, nbins = 310)
      P[epoch] <- poisson_fit$P[[1]]
      
    }
    
    # return all the bootstrapped t0 values
    return(list(P=P,inft=bs_data,t0=bs_data[, epoch]))
    
  }
  
  # Initialize t0 if not provided; represents initial infection days
  if (is.null(t0)) {
    null_indicator <- TRUE
    duration <- sample(days_survived, replace = TRUE, prob = normalised_probs, size = length(death_day))
    t0 <- death_day - duration
  } else {
    null_indicator <-FALSE
  }
  
  # Initialize result containers for Pearson values and incidence matrix
  P <- numeric(n.rep)                     # Stores Pearson values for each epoch
  inft <- matrix(0, nrow = 310, ncol = n.rep) # Stores incidence trajectory
  
  # Tabulate real death data for comparison with simulations
  death_day_310 <- tabulate(death_day, nbins = 310)
  
  # Main iteration loop over specified number of repetitions (epochs)
  durations <- sample(days_survived, replace = TRUE, prob = normalised_probs, size = length(t0)*n.rep)
  dim(durations)<-c(n.rep,length(t0))
  for (epoch in 1:n.rep) {
    
    # Generate new infection-to-death durations and simulated death days
    duration<-durations[epoch,]
    simulated_deaths <- t0 + duration
    simulated_deaths_310 <- tabulate(simulated_deaths, nbins = 310)
    
    # Compute Pearson statistic to evaluate the fit of simulated deaths to real deaths
    Pearson <- modified_pearson(death_day_310, simulated_deaths_310)
    
    # Generate sequence of possible shifts based on epoch count
    shift_amount <- generate_sequence(epoch)
    
    # Generate random sample and shifts for t0 adjustments
    random_sample<-sample(1:length(t0), size = length(t0), replace = FALSE)
    random_shifts<-sample(shift_amount, size = length(t0), replace = TRUE)

    # Adjust infection days in t0 to improve fit
    for (x in random_sample) {
      current_death_day <- simulated_deaths[x]
      adjusted_death_day_310 <- simulated_deaths_310
      
      # Adjust counts for moving infection time
      # Explanation: to avoid retabulating the simulated/death day vector,
      # we edit the existing vector in place.
      adjusted_death_day_310[current_death_day] <- adjusted_death_day_310[current_death_day] - 1
      adjusted_death_day_310[current_death_day + random_shifts[x]] <- adjusted_death_day_310[current_death_day + random_shifts[x]] + 1
      
      # Compute the partial Pearson statistic for this shift and accept if improvement is achieved
      
      # EXPLANATION: Only two values in the Pearson statistic are changed when
      # we move each infection day. These are the P value for the death days whose 
      # count is reduced and increased by one when we randomly shift the infection day.
      # Therefore, only the changed parts of the sum need to be recalculated.
      # These values can then be compared instead of recomputing the entire P
      # value. 
      # This adjustment to the pearson calculation causes a performance speedup
      
      orig_max1 <- max(c(1,simulated_deaths_310[current_death_day]))
      orig_p1 <-(death_day_310[current_death_day]- simulated_deaths_310[current_death_day])^2/orig_max1
      
      orig_max2 <- max(c(1,simulated_deaths_310[current_death_day+ random_shifts[x]]))
      orig_p2 <-(death_day_310[current_death_day+ random_shifts[x]]- simulated_deaths_310[current_death_day+ random_shifts[x]])^2/orig_max2
      orig_p <- orig_p1+orig_p2
      
      
      changed_max1 <- max(c(1,adjusted_death_day_310[current_death_day]))
      changed_p1 <-(death_day_310[current_death_day]- adjusted_death_day_310[current_death_day])^2/changed_max1
      
      changed_max2 = max(c(1,adjusted_death_day_310[current_death_day+ random_shifts[x]]))
      changed_p2 <-(death_day_310[current_death_day+ random_shifts[x]]- adjusted_death_day_310[current_death_day+ random_shifts[x]])^2/changed_max2
      
      changed_p <-changed_p1 +changed_p2
      if (length(changed_p < orig_p)==0||length(changed_p < orig_p)>1) {
        # out of range proposed change error
        next
      }
        if (changed_p < orig_p) {
          # Update t0,simulated deaths and Pearson loss with the improved values
          # if the new infection day reduces the loss.
          # NOTE: adjusting both t0 and simulated deaths keeps the duration of 
          # infection fixed, as desired.
          t0[x] <- max(t0[x] + random_shifts[x], 1)
          simulated_deaths[x] <- simulated_deaths[x] + random_shifts[x]
          simulated_deaths_310 <- adjusted_death_day_310
          Pearson <- Pearson - orig_p + changed_p
        }
      
    }
    
    # Store the Pearson statistic for this epoch
    P[epoch] <- Pearson
    
    
    # Record current infection estimates into inft
    inft[, epoch] <- tabulate(t0, nbins = 310)
    
    # Plot real vs. simulated data every 10 epochs for validation
    if(epoch %% 10 ==0 && bs==FALSE && null_indicator==TRUE){
    print(paste("epoch", epoch))
    print(paste("end P:", Pearson))
    plot(1:310, death_day_310, type = "l", col = "red", lwd = 2, ylim = c(0, 1700), 
         xlab = "Days", ylab = "Counts", main = paste("Epoch", epoch, "- Simulation State"))
    lines(1:310, simulated_deaths_310, col = "blue", lwd = 2,lty=5)
    lines(1:310, inft[,epoch], col = "black", lwd = 2)
    abline(v=84, col="purple")
    legend("topright", 
           legend = c("Estimated Incidence","Real Deaths","Simulated Deaths","First day of UK Lockdown"),
           lty = c(1,1,5,1),                 
           pch = c(NA,NA,NA,NA),
           col=c("black","red","blue","purple"))  
    }
  }
  
  # Return results containing loss values, incidence trajectory 
  # ,final t0 values and the simulated death days from the fitted model
  return(list(P = P, inft = inft, t0 = t0,sim_deaths=simulated_deaths_310))
}



# Run deconv function without bootstrapping to obtain initial results
returns <- deconv(covid_data$julian, covid_data$nhs, n.rep = 100, bs = FALSE, t0 = NULL)
# Run deconv with bootstrapping to estimate uncertainty
uncertainty_runs <- 100
bootstrap_t0 <- deconv(covid_data$julian, covid_data$nhs, n.rep = uncertainty_runs, bs = TRUE, t0 = returns$t0)

# get values for 95th percentiles of deaths on each of the 310 days
percentiles <- apply(bootstrap_t0$inft,1,quantile,probs=c(0.025,0.975))
lower_percentile = percentiles[1,] 
upper_percentile = percentiles[2,]

# Create final plot of 
# real/simulated deaths
death_day <- rep(covid_data$julian, times = covid_data$nhs)
plot(1:310,  tabulate(death_day, nbins = 310), type = "l", col = "red", lwd = 2, ylim = c(0, 1700), 
     xlab = "Days after 1st Jan 2020", ylab = "Counts", main = paste("Final Simulation State"))
lines(1:310, returns$sim_deaths, col = "blue", lwd = 2,lty=5)

# plot the accompanying confidence interval
polygon(x = c(1:310, rev(1:310)),
        y = c(upper_percentile, 
              rev(lower_percentile)),
        col =  adjustcolor("green", alpha.f = 0.5), border = NA)
# incidence, day of lockdown
lines(1:310, tabulate(returns$t0,nbins=310), col = "black", lwd = 2)
abline(v=84, col="purple")

legend("topright", 
       legend = c("Estimated Incidence","Real Deaths","Simulated Deaths","First day of UK Lockdown","95% Uncertainty Interval"),
       lty = c(1,1,5,1,1),                 
       pch = c(NA,NA,NA,NA,NA),
       col=c("black","red","blue","purple","green"))  


# Plot Pearson value fit plot
Sys.sleep(5)
plot(1:length(returns$P), returns$P, type = "l", col = "red", lwd = 2,
     xlab = "Epochs", ylab = "Modified Pearson Statistic", main = "Modified Pearson Statistic Fitness vs Epoch")
legend("topright", 
       legend = c("P(D,D^S)"),
       lty = c(1),                 
       pch = c(NA),
       col=c("red"))



# Experiment with changing shift values in generate_sequence function


# generate_sequence <- function(step_count) {
#  if (step_count < 10) {
#    return(c(-8, -4, -2, -1, 1, 2, 4, 8))
#  } else if (step_count < 30) {
#    return(c(-4, -2, -1, 1, 2, 4))
#  } else {
#    return(c(-2, -1, 1, 2))
#  }
#}


# Run deconv with the modified generate_sequence function to compare convergence
# returns_dif_step_lower <- deconv(covid_data$julian, covid_data$nhs, n.rep = 100, bs = FALSE, t0 = NULL)

# generate_sequence <- function(step_count) {
#  if (step_count < 35) {
#    return(c(-8, -4, -2, -1, 1, 2, 4, 8))
#  } else if (step_count < 65) {
#    return(c(-4, -2, -1, 1, 2, 4))
#  } else {
#    return(c(-2, -1, 1, 2))
#  }
# }

# returns_dif_step_higher <- deconv(covid_data$julian, covid_data$nhs, n.rep = 100, bs = FALSE, t0 = NULL)

# Overlay the Pearson convergence for the modified step sequence
# plot(1:length(returns$P), returns$P, type = "l", col = "red", lwd = 2,
     # xlab = "Epochs", ylab = "Pearson Statistic", main = "Pearson Statistic Fitness vs Epoch")
# lines(1:length(returns_dif_step_higher$P), returns_dif_step_higher$P, col ="blue", lwd = 2)
# lines(1:length(returns_dif_step_lower$P), returns_dif_step_lower$P, col = "yellow", lwd = 2)

## Observations
# After testing some different steps to change the random sequence by (and running it multiple times) we concluded that the 25,50 was the best choice
# as the pearson converged the fastest for the changing points 25,50. It performed worse noticeably for changing the steps on
# 10,30 potentially suggesting that the simulated death days were not yet as close to the actual ones.  The 35,65 performed similarly to the 25.50

