# HZAM-Sym_release_1.0.R
# Hybrid Zone Assortative Mating model, in sympatry with ecological differentiation.
# To model the effects of assortative mating, low hybrid fitness, and ecological differentiation on species coexistence.

# This software is a major modification and extension of the original HZAM code
# that was originally released by Darren Irwin, November 2019, University of British Columbia, Biodiversity Research Centre and Dept. of Zoology.
# That version was presented in the following paper, and is available at the associated Dryad archive listed below:

# Data produced using the original HZAM script is presented in the following paper: 
# Irwin, D.E. 2020. Assortative mating in hybrid zones is remarkably ineffective in promoting speciation. American Naturalist 195: E150-E167. https://doi.org/10.1086/708529

# Irwin, Darren (2019), Assortative mating in hybrid zones is remarkably ineffective in promoting speciation, Dryad, Dataset, https://doi.org/10.5061/dryad.k98sf7m30

# The present script, called HZAM-Sym, is a major modification of HZAM with substantial changes:
# This version models a non-spatial (i.e., purely sympatric) model of interactions 
# of two populations that have just come into contact, e.g. in a single lake or other habitat patch. 
# This means: 
# no dispersal curve (random uniform distribution of offspring "positions").
# no decline of competition over distance.
# HZAM-Sym incorporates ecological differentiation between the populations, and ecological competition.
# HZAM-Sym also has a choice of using either choice-based or "group-based" assortative mating. 

# HZAM-Sym is presented in the following manuscript:

# Irwin, D., and D. Schluter. 2021. Hybridization and the coexistence of species. bioRxiv 2021.04.04.438369; doi: https://doi.org/10.1101/2021.04.04.438369

# You are welcome to use and modify this script. If you do, please cite the manuscript above.

# Please note there are some options / capabilities in the code below that were not used in the simulations
# presented in the paper, but are included as options for future work.

# This file contains two main parts:
# Part 1: the code for running HZAM-Sym 
# Part 2: the code for generating the figures in the manuscript (from previously run and saved simulations)
# Part 3: additional code for generating other possibly useful figures from HZAM-Sym data

# Questions: Email Darren at irwin@zoology.ubc.ca


#### Part 1: The HZAM-Sym code ----

# Starting setup ----

# If these packages are not installed, run these lines to install:
# install.packages("viridis")
# install.packages("mgcv")
# install.packages("boot")
# install.packages("dplyr")
# install.packages("zoo")
####

# load packages:
library(viridis)  # for producing graphs using color-blind friendly palette
library(mgcv)  # needed for gam command
library(boot)  # needed for the inv.logit command
library(dplyr, warn.conflicts = F, quietly = T)  # for count, group_by, and summarize commands
library(zoo)  # for the rollapply command

# set the working directory (change as appropriate):
setwd("/Users/darrenirwin/Dropbox/Darren\'s\ current\ work/HZAM-Sym_project_and_paper")  

# set up colors that are good for color blind (found here: https://venngage.com/blog/color-blind-friendly-palette/)
colorBlindBlackPlasma41 <- c("black", plasma(40))
# to see these colours: 
# pie(rep(1, 41), col = colorBlindBlackPlasma41)

# Main script ----

# To run multiple replicates of the whole set, use this:
replications <- 1:5   # or just 1 for 1 replicate, or something like (2:5) to add replicates after 1 is done

set_name <- "setA"  # provide a name for this set of runs, which will be in the filenames

for (individual_replicate in replications) {
  
  # option to run a bunch of simulations consecutively
  run_set_name <- paste0(set_name,"_rep", individual_replicate)
  
  # maximum fitness of maximal heterozygote compared to pure forms
  hybrid_fitness_set <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 1)   # for just one run, just put one number in this and next line

  pref_ratio_set <- c(0, 0.001, 0.00333333, 0.01, 0.0333333, 0.1, 0.333333, 1)   # ratio in mating pref between hetero- and homospecific
  # Note that when pref_ratio is specified as 0 or 1 above, the code below changes those to numbers that allow the math to work (very close to 0 or 1) 
  
  mating_method <- "choice"  # options: "choice" means female preference--the main method in the paper; 
                                        # "grouping" means mating pools, without preference per se.
  
  if (mating_method=="grouping") {  # the settings here only needed if "grouping" mating is used  (I think choice is more realistic)
    mating_grouping_movement_ratio <- NULL  # when mating_method = "grouping", the fraction of species A that join mating pool of species B (and according to normal curve in between)
    mating_grouping_movement_SD <- NULL  # this will be set in code further down, based on mating_grouping_movement_ratio
    # To use a set of values for mating_grouping_movement_ratio that results in equivalent rate of production of F1s in the first generation 
    # as the pref_ratio produces in the choice-based mating, use the following for the corresponding pref-ratio: 
    mating_grouping_movement_ratio_set <- c(0,                # pref_ratio = 0
                                            0.0000000000075,  # pref_ratio = 0.001
                                            0.0000000027,     # pref_ratio = 0.00333333
                                            0.00000038,       # pref_ratio = 0.01
                                            0.000035,         # pref_ratio = 0.0333333
                                            0.00145,          # pref_ratio = 0.1
                                            0.039,            # pref_ratio = 0.333333
                                            0.8)              # pref_ratio = 1
  }
  
  track_fitness <- TRUE  # Set to T to keep track of realized fitness (offspring production) of the different trait values (i.e., HI categories) over time.
                         # Set to F to not track realized fitness, and get a boost in speed.
  
  # make sure the random number seed is not set, and instead is based on system time:
  set.seed(NULL)
  
  # Choose parameter states
  K <- 500     #500 # EVEN NUMBER; carrying capacity (on main resource alpha) of entire range (for two sexes combined), regardless of species
  K_half <- K/2
  K_B <- 500    #500 # EVEN NUMBER; carrying capacity (on resource beta) of entire range (for two sexes combined), regardless of species
  K_B_half <- K_B / 2
  N_0 <- K   # starting N of species 0
  N_0_half <- N_0/2
  N_1 <- K_B   # starting N of species 1
  N_1_half <- N_1/2
  
  ecolDiff <- 1
  
  individual_useResource_method <- 1   # 1 is linear gradient between species (the method used in the paper); 2 is Gaussian decline
  if (individual_useResource_method == 1) {   # The 4 values below can be fractional, and don't have to add to one for each species
    competAbility_useResourceA_species0 <- 0.5 + ecolDiff/2    # 1   
    competAbility_useResourceB_species0 <- 1 - competAbility_useResourceA_species0
    competAbility_useResourceA_species1 <- 0.5 - ecolDiff/2   # 0
    competAbility_useResourceB_species1 <- 1 - competAbility_useResourceA_species1
  } else if (individual_useResource_method == 2) {
    # If using Gaussian decline in resource use of each species (requires two species each on different resource peaks:
    # species0 at HI=0 using resource A, and species1 at HI=1 using resource B):
    width_resource_use_Gaussian <- 0.3
    competAbility_useResourceA_species0 <- 1
    competAbility_useResourceB_species0 <- exp(-((1 - 0)^2) / (2 * (width_resource_use_Gaussian^2)))
    competAbility_useResourceA_species1 <- exp(-((0 - 1)^2) / (2 * (width_resource_use_Gaussian^2)))
    competAbility_useResourceB_species1 <- 1
    if (FALSE) {  # option to make a graph of the fitness landscape for intermediate individuals
      quartz()
      spaced_HI <- round(seq(from = 0, to = 1, length.out = 1001), digits=3)  # the round is needed to ensure only 3 decimal places, avoiding error later
      ind_useResourceA_curve <- exp(-((spaced_HI - 0)^2) / (2 * (width_resource_use_Gaussian^2)))
      plot(spaced_HI, ind_useResourceA_curve, ylim=c(0,1.2))
      ind_useResourceB_curve <- exp(-((spaced_HI - 1)^2) / (2 * (width_resource_use_Gaussian^2)))
      points(spaced_HI, ind_useResourceB_curve)
      combined_useResource_curve <- ind_useResourceA_curve + ind_useResourceB_curve   # This is the underlying relationship between HI and total resource competAbility
      points(spaced_HI, combined_useResource_curve)
    }
  }
  
  max_generations <- 1000;  # the time at which the simulation will stop 
  
  do_plot <- T # whether to do plot during sims
  plot_MTL <- F  # whether to plot Mating Trait Loci
  plot_UDL <- T  # whether to plot UnderDominant Loci (these were the loci shown / modeled in the paper, and they are also Mating Trait Loci)
  plot_NL <- F  # whether to plot Neutral Loci  
  plot_density <- F   # T (or F) to show (not show) density plot
  plot_int <- 1 # plotting intervals, in generations  (also will fit cline at these intervals, if that is chosen below)
  make_movie <- F   # set whether to make a movie (T or F)
  movie_gen <- seq(5,250, by=5)  # applies only if making a movie
  movie_name <- "movie_temp"  # applies only if making a movie
  
  sympatry <- T   # T means space does not matter at all--this was used in the HZAM-sym paper; 
  # If T, then "location" is used purely as a random value that determines which male is encountered by female
  # If F, which is not recommended yet, then space does matter (in terms of dispersal and local competition, so the following lines matter)
  if (sympatry==F) {   # Note that the non-synpatric case has not been tested together with the ecological differentiation incorporated into HZAM-Sym
                       # and is not recommended. It is included here for future capability only. Take home message: keep sympatry set to T.
    meandispersal <- 10    # SD of dispersal curve  (not used by HZAM-Sym as presented in paper)
    sigma_comp <- 10 # the standard deviation (in units of space) of the density dependence effect  (not used by HZAM-Sym as presented in paper) 
  } 
  
  fit_to_cline <- F  # choose whether to fit a model at end (Not recommended for HZAM-Sym)
  fit_to_cline_during <- F  # choose whether to fit a cline model at each time plotted (Not recommended for HZAM-Sym)
  cline <- "UDL and NL"  # choose UDL or MTL or "UDL and NL" (for two clines)
  width_method <- "eightieth"   # choose "tangent" or "eightieth" percentile in middle  ("tangent" may not work well)
  
  survival_fitness_method <- 1  # option 1: underdominance only; option 2: epistasis (model of Barton & Gale 1993, Fig. 2-2) 
  beta <- 1  # Applies only if option 2 above: The beta value in the epistatic fitness equation: w(x) = 1 - s(4x[1-x])^beta
  
  range_limit_left <- 0   # in HZAM-Sym, non-spatial, this should always be 0
  range_limit_right <- 1  # in HZAM-Sym, non-spatial, this should always be 1
  
  per_reject_cost <- 0  # proportion reduction in fitness due to search cost, per rejected male (compounded). (not used in HZAM-Sym paper)
  
  growth_rate <- 1.05  # The variable R in the paper, this is the average maximum expected number of offspring per individual, when pop size far below K
  
  # toggles regarding whether to save and close a histogram of trait values at the end of each simulation:
  plot_histogram_at_end <- TRUE
  close_histogram <- TRUE   # To close the histogram after saving it
 
  # define starting ranges (in the case of HZAM-sym, both pops should cover the full "range" of zero to one)
  pop1_range_limit_left <- 0
  pop1_range_limit_right <- 1    #0.48
  pop2_range_limit_left <-  0    #0.52
  pop2_range_limit_right <- 1
  
  beginning_columns <- 1  # number of initial columns in matrix containing data such as location
  male_trait_loci <- 3  # number of loci determining male trait (and female trait too if assort_mating=1)
  assort_mating <- 1   # set for 1 if male trait and female trait are determined by same loci, as done in the paper (zero otherwise) 
  female_trait_loci <- 3  # if assort_mating = 1, then will automatically be equal to male_trait_loci, as they are the same loci; otherwise, this is number of loci determining female 'trait' (her ideal male trait)
  underdominant_loci <- 3  # number of loci in which hybrids selected against
  same_loci_MTL_UDL <- TRUE   # set to TRUE if the UDL loci are the same as the MTL loci; set to false if different loci
  neutral_loci <- 3  # number of neutral loci (used for neutral measure of hybrid index; not used in the HZAM-sym paper)
   
  mating_trait_loci_dominant <- F # set to T to make mating loci encode trait in completeley dominant / recessive way, Otherwise additive.
  
  # set up HI (Hybrid Index, equivalent to Trait value when no dominance) categories based on male_trait_loci
  HI_categories <- 2*male_trait_loci + 1
  HI_category_values <- round(seq(0, 1, by=1/(HI_categories-1)), digits=5)
  HI_category_breaks <- rep(NA, HI_categories-1)  # create empty vector
  for (i in 1:(HI_categories-1)) {
    HI_category_breaks[i] <- mean(HI_category_values[i:(i+1)]) 
  }
  category_labels <- round(HI_category_values, 2)
  low_category_lower_bound <- (0 - HI_category_breaks[1]) # This and next line used for the lower and upper bounds of the adjusted mating trait, for the "grouping" mating option
  high_category_upper_bound <- (1 + HI_category_breaks[1])
  
  # set colours for graphing (just for graphs shown during simulations, not those in paper)
  UDL_colour <- "seagreen"
  MTL_colour <- "purple"
  UDLMTL_colour <- "green4"  # tried a lot of greens, and green4 seems best
  NL_colour <- "gray70"
  #set colours of dots on graphs:
  if (same_loci_MTL_UDL == T) {
    UDL_colour <- UDLMTL_colour
    MTL_colour <- UDLMTL_colour
  }
  
  if (sympatry==F) {   # This option was not used in paper, and not recommended to use sympatry=F. Not yet fully tested, included here for future work.
    # set up locations every 0.001 across range, and calculate density of use on resource A when K individuals using resource A are perfectly spaced:
    spaced_locations <- round(seq(from = range_limit_left, to = range_limit_right, length.out = 1001), digits=3)  # the round is needed to ensure only 3 decimal places, avoiding error later
    ind_locations_if_even_at_K <- seq(from = range_limit_left, to = range_limit_right, length.out = K)
    get_density_if_even_at_K <- function(focal_location) {
      return(sum(exp(-((ind_locations_if_even_at_K - focal_location)^2)/(2*(sigma_comp^2)))))
    } 
    ideal_K_densities_at_spaced_locations <- sapply(spaced_locations, get_density_if_even_at_K)
    
    # set up locations every 0.001 across range, and calculate density of use on resource B when K individuals using resource B are perfectly spaced:
    ind_locations_if_even_at_K_B <- seq(from = range_limit_left, to = range_limit_right, length.out = K_B)
    get_density_if_even_at_K_B <- function(focal_location) {
      return(sum(exp(-((ind_locations_if_even_at_K_B - focal_location)^2)/(2*(sigma_comp^2)))))
    } 
    ideal_K_B_densities_at_spaced_locations <- sapply(spaced_locations, get_density_if_even_at_K_B)
  }
  
  ## Loop through the different simulation sets:
  
  for (hybrid_fitness_case in 1:length(hybrid_fitness_set)) {
    for (pref_ratio_case in 1:length(pref_ratio_set)) {
      
      hybrid_fitness <- hybrid_fitness_set[hybrid_fitness_case]
      pref_ratio <- pref_ratio_set[pref_ratio_case]
      if (mating_method=="grouping") {
        mating_grouping_movement_ratio <- mating_grouping_movement_ratio_set[pref_ratio_case]
      }

      run_name <- paste0("HZAM_animation_run",run_set_name,"_growthrate",growth_rate,"_ecolDiff",ecolDiff,"_K",K,"_UDLMTL",male_trait_loci,"_gen",max_generations,"_hybridfitness",hybrid_fitness,"_prefratio",pref_ratio)
      
      if (pref_ratio == 1) {  # this bit needed for calculation of pref_SD below, to avoid error when pref_ratio equals 1 or 0
        pref_ratio_for_math <- 1 - 10^(-15)  # very very close to 1
      } else if (pref_ratio == 0) {
        pref_ratio_for_math <- 10^(-30)  # very very close to 0
      } else {
        pref_ratio_for_math <- pref_ratio
      }
      pref_SD <- sqrt( -1 / (2 * log(pref_ratio_for_math)))     # width of female acceptance curve for male trait
      
      if (mating_method == "grouping") {
        if (mating_grouping_movement_ratio == 1) {  # this bit needed for calculation of mating_grouping_movement_SD below, to avoid error when mating_grouping_movement_ratio equals 1 or 0
          mating_grouping_movement_ratio_for_math <- 1 - 10^(-15)  # very very close to 1
        } else if (mating_grouping_movement_ratio == 0) {
          mating_grouping_movement_ratio_for_math <- 10^(-30)  # very very close to 0
        } else {
          mating_grouping_movement_ratio_for_math <- mating_grouping_movement_ratio
        }
        mating_grouping_movement_SD <- sqrt( -1 / (2 * log(mating_grouping_movement_ratio_for_math)))     # width of mating_grouping_movement curve, when "grouping" method of mating used
      }
      
      s_per_locus <- 1 - hybrid_fitness^(1/underdominant_loci)  # loss in fitness due to each heterozygous locus
      
      # Set up column structure for genotype matrix, based on number of loci of each type.
      # each locus has two alleles, so takes up two columns
      # for clarity later, calculate start and end columns for each locus type
      MTL_col_start <- beginning_columns + 1
      MTL_col_end <- beginning_columns + 2*male_trait_loci
      if (assort_mating == 1) {
        FTL_col_start <- MTL_col_start
        FTL_col_end <- MTL_col_end
      }
      if (assort_mating == 0) {
        FTL_col_start = MTL_col_end + 1
        FTL_col_end = MTL_col_end + 2*female_trait_loci
      }
      if (same_loci_MTL_UDL == TRUE) {
        UDL_col_start <- FTL_col_start
        UDL_col_end <- FTL_col_end
      }
      if (same_loci_MTL_UDL == FALSE) {
        UDL_col_start <- FTL_col_end + 1
        UDL_col_end <- FTL_col_end + 2*underdominant_loci
      }
      NL_col_start <- UDL_col_end + 1
      NL_col_end <- UDL_col_end + 2*neutral_loci
      total_loci <- (NL_col_end - beginning_columns)/2
      
      # Set up starting N for each species and total
      pop1_starting_N <- round(N_0_half * (pop1_range_limit_right - pop1_range_limit_left))  # N_0_half is for each sex
      pop2_starting_N <- round(N_1_half * (pop2_range_limit_right - pop2_range_limit_left))  # N_1_half is for each sex
      starting_N <- pop1_starting_N + pop2_starting_N  #for each sex
      
      # generate the population of females
      pop_matrix_F = matrix(-9, nrow=starting_N, ncol=1+2*total_loci) # create matrix to store population locations and genotypes; columns in this order: location, genotype columns
      # generate starting values for pop1:
      pop_matrix_F[1:pop1_starting_N,1] <- runif(n=pop1_starting_N, min=pop1_range_limit_left, max=pop1_range_limit_right)  # assigns random locations for pop1
      pop_matrix_F[1:pop1_starting_N,2:(1+2*total_loci)] <- 0   # assigns genotypes of pop1
      # generate starting values for pop2:
      pop_matrix_F[(1+pop1_starting_N):starting_N,1] <- runif(n=pop2_starting_N, min=pop2_range_limit_left, max=pop2_range_limit_right)  # assigns random locations for pop2
      pop_matrix_F[(1+pop1_starting_N):starting_N,2:(1+2*total_loci)] <- 1   # assigns genotypes of pop2
      
      # generate the population of males
      pop_matrix_M = matrix(-9, nrow=starting_N, ncol=1+2*total_loci) # create matrix to store population locations and genotypes; columns in this order: location, genotype columns
      # generate starting values for pop1:
      pop_matrix_M[1:pop1_starting_N,1] <- runif(n=pop1_starting_N, min=pop1_range_limit_left, max=pop1_range_limit_right)  # assigns random locations for pop1
      pop_matrix_M[1:pop1_starting_N,2:(1+2*total_loci)] <- 0   # assigns genotypes of pop1
      # generate starting values for pop2:
      pop_matrix_M[(1+pop1_starting_N):starting_N,1] <- runif(n=pop2_starting_N, min=pop2_range_limit_left, max=pop2_range_limit_right)  # assigns random locations for pop2
      pop_matrix_M[(1+pop1_starting_N):starting_N,2:(1+2*total_loci)] <- 1   # assigns genotypes of pop2
      
      # set up plot, if set to (recommended to use only use do_plot=T and plot_density=F):
      if (plot_density) {
        quartz(title=paste0("HZAM; pref_ratio=", pref_ratio,"; hybrid_fitness=", hybrid_fitness, sep=""), width=6, height=6, bg="white")
        dev_display <- dev.cur()
        par(oma=c(1,1,1,1))  # set outer margins
        zones <- matrix(c(2,1), ncol=1, byrow=TRUE)  # numbers in matrix give order of plotting
        layout(zones, widths=c(1,1), heights=c(1/3,2/3))  
      } else if (do_plot) {
        quartz(title=paste0("HZAM; pref_ratio=", pref_ratio,"; hybrid_fitness=", hybrid_fitness, sep=""), width=6, height=5, bg="white")
        dev_display <- dev.cur()
        plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(-0.05,1.05), yaxp=c(0,1,5), cex.axis=0.8, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
      }
      
      # if keeping track of realized fitness (offspring production) of the different trait values (i.e., HI categories):
      if (track_fitness) {
        mother_count_per_HI_by_gen <- matrix(nrow=length(HI_category_values), ncol=max_generations)
        offspring_count_per_motherHI_by_gen <- matrix(nrow=length(HI_category_values), ncol=max_generations)
        father_count_per_HI_by_gen <- matrix(nrow=length(HI_category_values), ncol=max_generations)
        offspring_count_per_fatherHI_by_gen <- matrix(nrow=length(HI_category_values), ncol=max_generations)
        matings_count_per_maleHI_by_gen <- matrix(nrow=length(HI_category_values), ncol=max_generations)
      }
      
      # set some initial values:
      time_vs_width = NULL
      extinction <- FALSE    # set this as the starting condition (population exists) 
      outcome <- NULL
      final_distribution <- NULL
      
      
      ###### Run an actual simulation, with cycles of mate choice, reproduction, dispersal
      
      for (time in 1:max_generations) {
        N <- dim(pop_matrix_F)[1] + dim(pop_matrix_M)[1]
        print(paste0("gen: ",time,"; PopSize: ",N))
        
        # determine male traits (for all males):
        if (mating_trait_loci_dominant == T) {  # the dominant case (not used in HZAM-sym paper)
          fraction_trait_loci_dominant <- function(x) {
            dominant_count <- 0
            # cycle through trait loci and count up the number of loci with a dominant (="1") allele
            for (locus_count in 1:male_trait_loci) {
              locus_columns <- (2*(locus_count-1))+(MTL_col_start:((MTL_col_start)+1))
              if (sum(x[locus_columns]) %in% c(1,2)) {
                dominant_count <- dominant_count + 1
              }
            }
            return(dominant_count / male_trait_loci)  
          } 
          male_traits <- apply(pop_matrix_M, 1, fraction_trait_loci_dominant) # the dominant case
        } else {  # the additive case (used in paper)
          male_traits <- apply(pop_matrix_M, 1, function(x) mean(x[MTL_col_start:MTL_col_end]))  
        }
        
        # determine female preferences (for all females):
        if (mating_trait_loci_dominant == T) {  # the dominant case
          fraction_trait_loci_dominant <- function(x) {
            dominant_count <- 0
            # cycle through trait loci and count up the number of loci with a dominant (="1") allele
            for (locus_count in 1:female_trait_loci) {
              locus_columns <- (2*(locus_count-1))+(FTL_col_start:((FTL_col_start)+1))
              if (sum(x[locus_columns]) %in% c(1,2)) {
                dominant_count <- dominant_count + 1
              }
            }
            return(dominant_count / female_trait_loci)  
          } 
          female_traits <- apply(pop_matrix_F, 1, fraction_trait_loci_dominant) # the dominant case
        } else {  # the additive case (used in paper)
          female_traits <- apply(pop_matrix_F, 1, function(x) mean(x[FTL_col_start:FTL_col_end]))  # the additive case
        }
        
        # New way of calculating density dependence, using similar method as Irwin 2002 AmNat,
        # but with each individual counting toward portions of resource A and/or B resource use:
        # Two options for how intermediate individuals use resources:
        # - According to linear gradient between use of species 0 and species 1 (used in the HZAM-sym paper)
        # - According to Gaussian curve centered on their trait value  (not used in HZAM-sym paper)
        if (individual_useResource_method == 1) {
          # trait of zero counts as resource use according to competAbility_useResourceA_species0 and competAbility_useResourceB_species0, 
          # trait of 1 counts as resource use according to competAbility_useResourceA_species1 and competAbility_useResourceB_species1,
          # and in between is linear gradient:
          ind_useResourceA_F <- competAbility_useResourceA_species1 + ((1-female_traits)*(competAbility_useResourceA_species0 - competAbility_useResourceA_species1))
          ind_useResourceB_F <- competAbility_useResourceB_species0 + (female_traits * (competAbility_useResourceB_species1 - competAbility_useResourceB_species0))
          ind_useResourceA_M <- competAbility_useResourceA_species1 + ((1-male_traits)*(competAbility_useResourceA_species0 - competAbility_useResourceA_species1))
          ind_useResourceB_M <- competAbility_useResourceB_species0 + (male_traits * (competAbility_useResourceB_species1 - competAbility_useResourceB_species0))
        } else if (individual_useResource_method == 2) {
          # Alternate way of counting resource use of individuals, according to Gaussian curve centered on their trait value:
          ind_useResourceA_F <- exp(-((female_traits - 0)^2) / (2 * (width_resource_use_Gaussian^2)))
          ind_useResourceB_F <- exp(-((female_traits - 1)^2) / (2 * (width_resource_use_Gaussian^2)))
          ind_useResourceA_M <- exp(-((male_traits - 0)^2) / (2 * (width_resource_use_Gaussian^2)))
          ind_useResourceB_M <- exp(-((male_traits - 1)^2) / (2 * (width_resource_use_Gaussian^2)))
        }
        
        ind_useResourceA <- c(ind_useResourceA_F, ind_useResourceA_M)
        ind_useResourceB <- c(ind_useResourceB_F, ind_useResourceB_M)
        
        if (sympatry==F) {  # do this only if space matters (not used in HZAM-sym paper)--this is not fully tested and use is not recommended--included here for future work
          ind_locations_real <- c(pop_matrix_F[,1] , pop_matrix_M[,1])
          get_density_useResourceA_real <- function(focal_location) {
            ind_density_values <- exp(-((ind_locations_real - focal_location)^2)/(2*(sigma_comp^2)))
            return(sum(ind_density_values * ind_useResourceA))
          } 
          real_density_useResourceA_at_spaced_locations <- sapply(spaced_locations, get_density_useResourceA_real)
          
          get_density_useResourceB_real <- function(focal_location) {
            ind_density_values <- exp(-((ind_locations_real - focal_location)^2)/(2*(sigma_comp^2)))
            return(sum(ind_density_values * ind_useResourceB))
          } 
          real_density_useResourceB_at_spaced_locations <- sapply(spaced_locations, get_density_useResourceB_real)
          
          # Liou & Price (1994) equation, where density dependence based on same K everywhere, but limited growth rate:
          local_growth_rates_ResourceA <- growth_rate*ideal_K_densities_at_spaced_locations / (ideal_K_densities_at_spaced_locations + ((real_density_useResourceA_at_spaced_locations)*(growth_rate - 1)))
          local_growth_rates_ResourceB <- growth_rate*ideal_K_B_densities_at_spaced_locations / (ideal_K_B_densities_at_spaced_locations + ((real_density_useResourceB_at_spaced_locations)*(growth_rate - 1)))
        } else if (sympatry==T) {  # if space does not matter, i.e. pure sympatry (used in HZAM-sym paper)
          global_useResourceA <- sum(ind_useResourceA)
          global_useResourceB <- sum(ind_useResourceB)
          global_growth_rate_ResourceA <- growth_rate*K / (K + ((global_useResourceA)*(growth_rate - 1)))
          global_growth_rate_ResourceB <- growth_rate*K / (K + ((global_useResourceB)*(growth_rate - 1)))
        }
        
        # cycle through the mothers, determining male mates and numbers of offspring (and their genotypes and locations)
        daughters_per_mother_list <- replicate(n=nrow(pop_matrix_F), matrix(nrow=0, ncol=1+2*total_loci))  # creates a list with identical empty matrices (with zero rows)
        sons_per_mother_list <- replicate(n=nrow(pop_matrix_F), matrix(nrow=0, ncol=1+2*total_loci))
        # Also keep track of offspring per father, to determine reproductive fitness of fathers:
        daughters_per_father_list <- replicate(n=nrow(pop_matrix_M), matrix(nrow=0, ncol=1+2*total_loci)) 
        sons_per_father_list <- replicate(n=nrow(pop_matrix_M), matrix(nrow=0, ncol=1+2*total_loci))
        
        # Record number of matings per male, to determine sexual selection due to HI class:
        matings_per_male <- rep(0, times=nrow(pop_matrix_M))  # this will record the number of matings per male
        matings_per_female <- rep(0, times=nrow(pop_matrix_F))  # this will record the number of matings per female
       
         # now cycle through all females and determine male mate, using either "choice" or "grouping" for mating method
        if (mating_method == "choice") { 
          for (mother in 1:nrow(pop_matrix_F))  {
            # Find a male mate
            # rank males in terms of "distance" from female (in case of sympatry, not real distance but just the difference in their random "location" values)
            male_dists <- abs(pop_matrix_F[mother,1] - pop_matrix_M[,1])  # calculates vector of "dist" of each male from focal female
            # pick closest male, determine his male signaling trait
            mate <- 0
            rejects <- 0
            father <- NULL  # this and above two lines just initialize variables
            while (mate == 0) {
              focal_male <- which(male_dists == min(male_dists))[1]  # returns row number of male that is "closest"
              # compare male trait with female's trait (preference), and determine
              # whether she accepts; note that match_strength is determined by a
              # Gaussian, with a maximum of 1 and minimum of zero.
              match_strength <- (exp(1) ^ ((-(male_traits[focal_male] - female_traits[mother])^2) / (2 * (pref_SD ^2))))
              if (runif(1) < match_strength) {
                # she accepts male, and they reproduce
                father <- focal_male
                matings_per_male[focal_male] <- matings_per_male[focal_male] + 1
                matings_per_female[mother] <- matings_per_female[mother] + 1
                mate <- 1
              } else {
                # she rejects male;
                # change that male's distance to a very large number (99) so he
                # won't be considered again (just in this bit)
                male_dists[focal_male] <- 99
                rejects <- rejects + 1  # count number of rejects, for imposition of fitness cost for search time
                if (all(male_dists == 99)) {  # if the female has rejected all males, then she doesn't breed
                  father <- NaN
                  break
                } 
              }
            }
            # Reproduce
            # determine fitness cost due to mate search time
            search_fitness <- (1-per_reject_cost) ^ rejects    # calculate proportion fitness lost (1-cost) due to mate search time (no cost in HZAM-sym paper)
            if (sympatry==F) {  # determine growth rate based on both resources A and B at mother's location:
              location_rounded <- round(pop_matrix_F[mother,1], digits=3)
              local_growth_dueToResourceA <- local_growth_rates_ResourceA[spaced_locations==location_rounded]
              local_growth_dueToResourceB <- local_growth_rates_ResourceB[spaced_locations==location_rounded]
              # weighted (based on mother's trait) average of growth rates on resources A and B: 
              growth_rate_of_focal_female <- (ind_useResourceA_F[mother] * local_growth_dueToResourceA) + (ind_useResourceB_F[mother] * local_growth_dueToResourceB) 
            } else if (sympatry==T) {  # determine growth rate based on global densities
              growth_rate_of_focal_female <- (ind_useResourceA_F[mother] * global_growth_rate_ResourceA) + (ind_useResourceB_F[mother] * global_growth_rate_ResourceB) 
            }
            #combine for total fitness:   
            reproductive_fitness <- 2 * growth_rate_of_focal_female * search_fitness  # the 2 is because only females, not males, produce offspring
            # now draw the number of offspring from a poisson distribution with a mean of total_fitness
            if (is.nan(father)) {  # if the female rejected all males, then she has no offspring
              offspring <- 0
            } else {  # else if the female accepted a male, then she draws offpsring number from a Poissson distribution with mean reproductive_fitness
              offspring <- rpois(1, reproductive_fitness) 
            }
            daughters <- matrix(nrow=0, ncol=1+2*total_loci) # initialize daughters and sons variables as empty matrices
            sons <- matrix(nrow=0, ncol=1+2*total_loci)
            if (offspring >= 1)  { # if offspring, generate their location and genotypes
              for (kid in 1:offspring) {
                kid_info <- rep(-9, 1+2*total_loci)
                if (sympatry==F) {  # if space matters, disperse the kid according to Gaussian dispersal curve, with sd=meandispersal
                  while (1 == 1)  {	
                    newloc <- pop_matrix_F[mother,1] + rnorm(1, mean=0, sd=meandispersal) 
                    if ((newloc <= range_limit_right) & (newloc >= range_limit_left))  {
                      break
                    }
                  }
                } else if (sympatry==T) {  # in sympatric case, kid has random "location" in range (from uniform random distribution)
                  newloc <- runif(1, min=range_limit_left, max=range_limit_right)
                }
                
                kid_info[1] <- newloc
                # generate genotypes; for each locus, first column for allele from mother, second for allele from father
                for (locus in 1:total_loci) {   
                  # choose allele from mother
                  kid_info[2+2*(locus-1)] <- pop_matrix_F[mother,1+2*(locus-1)+sample(2, 1)]
                  # choose allele from father
                  kid_info[3+2*(locus-1)] <- pop_matrix_M[father,1+2*(locus-1)+sample(2, 1)]
                }
                # determine sex of kid and add to table
                if (runif(1) > 0.5) {
                  daughters <- rbind(daughters, kid_info)
                } else {
                  sons <- rbind(sons, kid_info)
                }
              } 
              daughters_per_mother_list[[mother]] <- rbind(daughters_per_mother_list[[mother]], daughters)
              sons_per_mother_list[[mother]] <- rbind(sons_per_mother_list[[mother]], sons)
              daughters_per_father_list[[father]] <- rbind(daughters_per_father_list[[father]], daughters)  # this and next line can be used to track fitness of fathers
              sons_per_father_list[[father]] <- rbind(sons_per_father_list[[father]], sons)
            } 
          } 
          pop_matrix_daughters <- do.call("rbind", daughters_per_mother_list)
          pop_matrix_sons <- do.call("rbind", sons_per_mother_list)
        } 
        else if (mating_method == "grouping") {
          # assign females to mating groups
          female_traits_modified_for_group_assignment <- rep(-9, times=length(female_traits))
          for (female in 1:length(female_traits)) {
            success <- F
            value <- NULL
            while (success == F) {
              value <- female_traits[female] + rnorm(1, mean=0, sd=mating_grouping_movement_SD)
              if ((value > low_category_lower_bound) & (value < high_category_upper_bound)) {
                success <- T
              }
            }
            female_traits_modified_for_group_assignment[female] <- value
          }
          female_mating_group <- round(6 * female_traits_modified_for_group_assignment, digits = 0) + 1
          
          # assign males to mating groups
          male_traits_modified_for_group_assignment <- rep(-9, times=length(male_traits))
          for (male in 1:length(male_traits)) {
            success <- F
            value <- NULL
            while (success == F) {
              value <- male_traits[male] + rnorm(1, mean=0, sd=mating_grouping_movement_SD)
              if ((value > low_category_lower_bound) & (value < high_category_upper_bound)) {
                success <- T
              }
            }
            male_traits_modified_for_group_assignment[male] <- value
          }
          male_mating_group <- round(6 * male_traits_modified_for_group_assignment, digits = 0) + 1
          
          # Now, within each group, pair females and males. Each individual can mate only once. 
          # There will be some unpaired, but this should not be biased towards females or males,
          # or towards any particular values of HI.
          
          # Set up record of who has been mated:
          # female_mated <- rep(FALSE, times=length(male_mating_group))
          # male_mated <- rep(FALSE, times=length(male_mating_group))
          
          # determine number of pairings per mating group:
          pairings <- rep(-9, times=HI_categories)
          unpaired_females <- rep(-9, times=HI_categories)
          unpaired_males <- rep(-9, times=HI_categories)
          for (group in 1:(HI_categories)) {
            pairings[group] <- min(sum(female_mating_group == group), sum(male_mating_group == group))
            unpaired_females[group] <- sum(female_mating_group == group) - pairings[group]
            unpaired_males[group] <- sum(male_mating_group == group) - pairings[group]
          }
          #growth_rate_correction_due_to_unpaired_females <- (pairings + unpaired_females) / pairings
          growth_rate_correction_due_to_unpaired_individuals <- (2*pairings + unpaired_females + unpaired_males) / (2*pairings)
          total_pairings <- sum(pairings)
          
          # check if no pairings, and end the simulation if so
          if (total_pairings == 0) {
            extinction <- TRUE   # record an extinction of whole population (both "species")
            break    # break out of current loop (this simulation)
          }
          
          # for each mating group, pair up males and females randomly (once per individual)
          mate_order_females <- NULL
          mate_order_males <- NULL
          ind_female_rep_rate_correction <- NULL
          for (group in 1:(HI_categories)) {
            if (pairings[group] >= 1) {
              random_order_female_indices_in_mating_group <- sample(which(female_mating_group == group))
              random_order_male_indices_in_mating_group <- sample(which(male_mating_group == group))
              mate_order_females <- c(mate_order_females, random_order_female_indices_in_mating_group[1:pairings[group]])
              mate_order_males <- c(mate_order_males, random_order_male_indices_in_mating_group[1:pairings[group]])
              #ind_female_growth_rate_correction <- c(ind_female_growth_rate_correction, rep(growth_rate_correction_due_to_unpaired_females[group], times=pairings[group]))
              ind_female_rep_rate_correction <- c(ind_female_rep_rate_correction, rep(growth_rate_correction_due_to_unpaired_individuals[group], times=pairings[group]))
            }
          }
          
          matings_per_male[mate_order_males] <- 1
          matings_per_female[mate_order_females] <- 1
          
          # Reproduce according to mate pairings
          for (pair_num in 1:total_pairings) {
            mother <- mate_order_females[pair_num]
            father <- mate_order_males[pair_num]
            search_fitness <- 1    # a placeholder in case want to add some sort of search cost later (not clear how in case of group mating)
            if (sympatry==F) {  # determine growth rate based on both resources A and B at mother's location:
              location_rounded <- round(pop_matrix_F[mother,1], digits=3)
              local_growth_dueToResourceA <- local_growth_rates_ResourceA[spaced_locations==location_rounded]
              local_growth_dueToResourceB <- local_growth_rates_ResourceB[spaced_locations==location_rounded]
              # weighted (based on mother's trait) average of growth rates on resources A and B: 
              growth_rate_of_focal_female <- (ind_useResourceA_F[mother] * local_growth_dueToResourceA) + (ind_useResourceB_F[mother] * local_growth_dueToResourceB) 
            } else if (sympatry==T) {  # determine growth rate based on global densities
              growth_rate_of_focal_female <- (ind_useResourceA_F[mother] * global_growth_rate_ResourceA) + (ind_useResourceB_F[mother] * global_growth_rate_ResourceB) 
            }
            #combine for total fitness (the 2 is because only females, not males, produce offspring):   
            reproductive_fitness <- 2 * growth_rate_of_focal_female * search_fitness * ind_female_rep_rate_correction[pair_num] 
            # now draw the number of offspring from a poisson distribution with a mean of total_fitness
            offspring <- rpois(1, reproductive_fitness) 
            daughters <- matrix(nrow=0, ncol=1+2*total_loci)
            sons <- matrix(nrow=0, ncol=1+2*total_loci)
            if (offspring >= 1)  { # if offspring, generate their location and genotypes
              for (kid in 1:offspring) {
                kid_info <- rep(-9, 1+2*total_loci)
                if (sympatry==F) {  # disperse the kid according to Gaussian dispersal curve, with sd=meandispersal
                  while (1 == 1)  {	
                    newloc <- pop_matrix_F[mother,1] + rnorm(1, mean=0, sd=meandispersal)
                    if ((newloc <= range_limit_right) & (newloc >= range_limit_left))  {
                      break
                    }
                  }
                } else if (sympatry==T) {  # kid has random location in range (from uniform random distribution)
                  newloc <- runif(1, min=range_limit_left, max=range_limit_right)
                }
                
                kid_info[1] <- newloc
                # generate genotypes; for each locus, first column for allele from mother, second for allele from father
                for (locus in 1:total_loci) {   
                  # choose allele from mother
                  kid_info[2+2*(locus-1)] <- pop_matrix_F[mother,1+2*(locus-1)+sample(2, 1)]
                  # choose allele from father
                  kid_info[3+2*(locus-1)] <- pop_matrix_M[father,1+2*(locus-1)+sample(2, 1)]
                }
                # determine sex of kid and add to table
                if (runif(1) > 0.5) {
                  daughters <- rbind(daughters, kid_info)
                } else {
                  sons <- rbind(sons, kid_info)
                }
              } 
              daughters_per_mother_list[[mother]] <- rbind(daughters_per_mother_list[[mother]], daughters)
              sons_per_mother_list[[mother]] <- rbind(sons_per_mother_list[[mother]], sons)
              daughters_per_father_list[[father]] <- rbind(daughters_per_father_list[[father]], daughters)  # this and next line can be used to track fitness of fathers
              sons_per_father_list[[father]] <- rbind(sons_per_father_list[[father]], sons)
            } 
          } 
          pop_matrix_daughters <- do.call("rbind", daughters_per_mother_list)
          pop_matrix_sons <- do.call("rbind", sons_per_mother_list)
        }
          
        # check if either no daughters or no sons, and end the simulation if so
        if (nrow(pop_matrix_daughters)==0 | nrow(pop_matrix_sons)==0) {
          extinction <- TRUE   # record an extinction of whole population (both "species")
          break    # break out of current loop (this simulation)
        }
        
        if (track_fitness) {    # keep track of number of offspring per HI types of parent
          offspring_per_mother <- sapply(daughters_per_mother_list, nrow) + sapply(sons_per_mother_list, nrow)
          offspring_per_father <- sapply(daughters_per_father_list, nrow) + sapply(sons_per_father_list, nrow)
          # count mothers in each HI category, then offspring per category
          mother_traits_and_offspring <- data.frame(HI=round(female_traits, digits=5), offspring=offspring_per_mother)
          mother_traits_and_offspring$HI <- factor(mother_traits_and_offspring$HI, levels = HI_category_values)  # assign factor levels to HI
          mother_count_per_HI <- mother_traits_and_offspring %>%
            count(HI, .drop = FALSE)
          offspring_count_per_motherHI <- mother_traits_and_offspring %>%
            group_by(HI, .drop = FALSE) %>%
            summarize(total_offspring=sum(offspring))
          # save the counts of mothers and offspring per HI category
          mother_count_per_HI_by_gen[,time] <- mother_count_per_HI$n
          offspring_count_per_motherHI_by_gen[,time] <- offspring_count_per_motherHI$total_offspring
          # count fathers in each HI category, offspring per category, and matings per category
          father_traits_and_offspring <- data.frame(HI=round(male_traits, digits=5), offspring=offspring_per_father, matings=matings_per_male)
          father_traits_and_offspring$HI <- factor(father_traits_and_offspring$HI, levels = HI_category_values)  # assign factor levels to HI
          father_count_per_HI <- father_traits_and_offspring %>%
            count(HI, .drop = FALSE)
          offspring_count_per_fatherHI <- father_traits_and_offspring %>%
            group_by(HI, .drop = FALSE) %>%
            summarize(total_offspring=sum(offspring))
          matings_count_per_fatherHI <- father_traits_and_offspring %>%
            group_by(HI, .drop = FALSE) %>%
            summarize(total_matings=sum(matings))
          offspring_and_HIfather_counts <- data.frame(HI=father_count_per_HI$HI, father_count=father_count_per_HI$n,
                                                      total_offspring=offspring_count_per_fatherHI$total_offspring,
                                                      total_matings=matings_count_per_fatherHI$total_matings)
          # save the counts of fathers and offspring and matings per HI category
          father_count_per_HI_by_gen[,time] <- father_count_per_HI$n
          offspring_count_per_fatherHI_by_gen[,time] <- offspring_count_per_fatherHI$total_offspring
          matings_count_per_maleHI_by_gen[,time] <- matings_count_per_fatherHI$total_matings
        }
        
        # determine survival fitnesses of daughters due to heterozygosity at (option 1) heterozygosity, or (option 2) epistasis 
        if (survival_fitness_method == 1) {  # option 1 is survival probability based on heterozygosity (in HZAM-sym paper)
          underdominance_fitness_daughters <- rep(NaN, dim(pop_matrix_daughters)[1])  # initialize data structure
          for (daughter in 1:nrow(pop_matrix_daughters))  {
            heterozyg_loci <- 0
            for (locus in 1:underdominant_loci) {
              locus_col <- UDL_col_start + 2*(locus-1)
              if (mean(pop_matrix_daughters[daughter,locus_col:(locus_col+1)]) == 0.5 ) { 
                heterozyg_loci <- heterozyg_loci + 1
              }
            }
            underdominance_fitness_daughters[daughter] <- (1-s_per_locus) ^ heterozyg_loci
          }
          # determine whether each daughter survives to adulthood
          random_proportions_daughters <- runif(n=length(underdominance_fitness_daughters), min=0, max=1)
          daughters_survive <- underdominance_fitness_daughters > random_proportions_daughters
          
          # determine survival fitnesses of sons due to heterozygosity at underdominance loci
          underdominance_fitness_sons <- rep(NaN, dim(pop_matrix_sons)[1])
          for (son in 1:nrow(pop_matrix_sons))  {
            heterozyg_loci <- 0
            for (locus in 1:underdominant_loci) {
              locus_col <- UDL_col_start + 2*(locus-1)
              if (mean(pop_matrix_sons[son,locus_col:(locus_col+1)]) == 0.5 ) { 
                heterozyg_loci <- heterozyg_loci + 1
              }
            }
            underdominance_fitness_sons[son] <- (1-s_per_locus) ^ heterozyg_loci
          } 
          # determine whether each son survives to adulthood
          random_proportions_sons <- runif(n=length(underdominance_fitness_sons), min=0, max=1)
          sons_survive <- underdominance_fitness_sons > random_proportions_sons
          
        } else if (survival_fitness_method == 2) {  # option 2 is survival probability based on epistasis (not in HZAM-sym paper)
          HI_daughters <- rowMeans(pop_matrix_daughters[,UDL_col_start:UDL_col_end]) 
          epistasis_fitness_daughters <- 1 - (1-hybrid_fitness) * (4*HI_daughters*(1-HI_daughters))^beta
          random_proportions_daughters <- runif(n=length(epistasis_fitness_daughters), min=0, max=1)
          daughters_survive <- epistasis_fitness_daughters > random_proportions_daughters
          
          HI_sons <- rowMeans(pop_matrix_sons[,UDL_col_start:UDL_col_end])
          epistasis_fitness_sons <- 1 - (1-hybrid_fitness) * (4*HI_sons*(1-HI_sons))^beta
          random_proportions_sons <- runif(n=length(epistasis_fitness_sons), min=0, max=1)
          sons_survive <- epistasis_fitness_sons > random_proportions_sons  
        }
        
        # assign surviving offspring to new adult population
        pop_matrix_F <- pop_matrix_daughters[daughters_survive,,drop=F]  # the "drop=F" prevents R from collapsing a 1-row matrix to a vector
        pop_matrix_M <- pop_matrix_sons[sons_survive,,drop=F]
        
        # check if either no surviving females or no surviving males, and end the simulation if so
        if (nrow(pop_matrix_F)==0 | nrow(pop_matrix_M)==0) {
          extinction <- TRUE   # record an extinction of whole population (both "species")
          break    # break out of current loop (this one simulation)
        }
        
        # Option of fitting cline during the simulation. 
        # Not used and not meaningful in sympatric case, so this not used in HZAM-sym paper.
        # (Was used in the original HZAM paper in AmNat 2020, where space mattered)
        if (fit_to_cline_during & ((time %% plot_int) == 0)) {
          # choose the data to include
          
          x = c(pop_matrix_F[,1] , pop_matrix_M[,1])
          if (cline=="UDL") {
            y <- c(rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end]) , rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end]))
          }
          else if (cline=="MTL") {
            y <- c(rowMeans(pop_matrix_F[,MTL_col_start:MTL_col_end]) , rowMeans(pop_matrix_M[,MTL_col_start:MTL_col_end]))
          }
          else if (cline=="UDL and NL") {
            y <- c(rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end]) , rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end]))
            y_NL <- c(rowMeans(pop_matrix_F[,NL_col_start:NL_col_end]) , rowMeans(pop_matrix_M[,NL_col_start:NL_col_end]))
          }
          
          # fit UDL or MTL to cline
          mydata <- as.data.frame(cbind(x,y))
          z <- gam(y ~ s(x), data=  mydata, quasibinomial(link = "logit"),
                   method = "ML")
          x.spaced <- seq(0, 1, by=0.001)
          y.predicted <- inv.logit(as.vector(predict.gam(z, newdata=data.frame(x=x.spaced))))
          #thanks to Dolph Schluter for his help with above                      
          
          if (width_method=="tangent") {
            #find the x with maximum slope
            first.diff <- diff(y.predicted)
            index <- which.max(first.diff)
            x.1 <- x.spaced[index]
            x.2 <- x.spaced[index+1]
            y.1 <- y.predicted[index]
            y.2 <- y.predicted[index+1]
            max.slope <- (y.2 - y.1) / (x.2 - x.1)
            width <- 1/max.slope   # if asymptotes at 0 and 1
            print(paste0("width = ",round(width, digits=3)))
            intercept <- y.1 - max.slope*x.1
            X.high <- (1 - intercept) / max.slope
            X.low <- (0 - intercept) / max.slope
            
          }
          
          if (width_method=="eightieth") {
            # find place that fit goes below 0.05 on left of center
            left_width_margin <- max(x.spaced[y.predicted <= 0.1])
            # find place where fit goes above 0.95 on right of center
            right_width_margin <- min(x.spaced[y.predicted >= 0.9])
            dist_from_HI50percent <- abs(y.predicted-0.5)
            centre <- x.spaced[dist_from_HI50percent == min(dist_from_HI50percent)]
            width <- right_width_margin - left_width_margin
            print(paste0("width = ",round(width, digits=3)))
          }
          
          # this part is for the NL clinefit and width
          if (cline=="UDL and NL") {     
            mydata <- as.data.frame(cbind(x,y_NL))
            z <- gam(y_NL ~ s(x), data=  mydata, quasibinomial(link = "logit"),
                     method = "ML")
            NL.x.spaced <- seq(0, 1, by=0.001)
            NL.y.predicted <- inv.logit(as.vector(predict.gam(z, newdata=data.frame(x=x.spaced))))
            #thanks to Dolph Schluter for his help with above                      
            
            if (width_method=="tangent") {
              #find the x with maximum slope
              first.diff <- diff(NL.y.predicted)
              index <- which.max(first.diff)
              x.1 <- NL.x.spaced[index]
              x.2 <- NL.x.spaced[index+1]
              y.1 <- NL.y.predicted[index]
              y.2 <- NL.y.predicted[index+1]
              max.slope <- (y.2 - y.1) / (x.2 - x.1)
              width_NL <- 1/max.slope   # if asymptotes at 0 and 1
              print(paste0("width_NL = ",round(width_NL, digits=3)))
              intercept <- y.1 - max.slope*x.1
              NL.X.high <- (1 - intercept) / max.slope
              NL.X.low <- (0 - intercept) / max.slope
            }
            
            if (width_method=="eightieth") {
              # find place that fit goes below 0.05 on left of center
              NL.left_width_margin <- max(NL.x.spaced[NL.y.predicted <= 0.1])
              # find place where fit goes above 0.95 on right of center
              NL.right_width_margin <- min(NL.x.spaced[NL.y.predicted >= 0.9])
              dist_from_HI50percent <- abs(NL.y.predicted-0.5)
              NL.centre <- x.spaced[dist_from_HI50percent == min(dist_from_HI50percent)]
              width_NL <- NL.right_width_margin - NL.left_width_margin
              print(paste0("width_NL = ",round(width_NL, digits=3)))
            }
          }
          
          # add time and width to matrix
          if (cline=="UDL" | cline =="MTL") {
            time_vs_width <- rbind(time_vs_width, c(time, width, centre, left_width_margin, right_width_margin))
          } else if (cline=="UDL and NL") {
            time_vs_width <- rbind(time_vs_width, c(time, width, width_NL, centre, NL.centre, left_width_margin, NL.left_width_margin, right_width_margin, NL.right_width_margin))
          }
        }
        
        if (do_plot == 1) {  # will show plot of the simulation
          if ((time %% plot_int) == 0)  { # update figures at regular intervals (x%%y means remainder of x/y)    
            if (plot_density) {
              par(mar=c(3,3,1,1))  # specifies number of lines around plot (bottom, left, top right)
            }
            
            plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(-0.05,1.05), yaxp=c(0,1,5), cex.axis=0.8, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
            title(ylab="Hybrid index", line=2, cex.lab=1.2)
            if (sympatry==F) {
              title(xlab="Location", line=2, cex.lab=1.2)
            } else if (sympatry) {
              title(xlab="Random value (just to spread population visually)", line=2, cex.lab=1.2)
            }
            
            if (fit_to_cline_during) {
              if (width_method=="tangent") {
                # add fit to plot
                lines(c(X.high, X.low), c(1, 0), col=adjustcolor("grey", alpha.f = 0.5), lwd=4)
                lines(x.spaced, y.predicted, col=adjustcolor("blue", alpha.f = 0.75), lwd=4)
              }
              
              if (width_method=="eightieth") {
                polygon(x=c(left_width_margin, right_width_margin, right_width_margin, left_width_margin), y=c(-0.05, -0.05, 1.05, 1.05), col=adjustcolor("skyblue", alpha.f = 0.5), border=NA)
                lines(x.spaced, y.predicted, col=adjustcolor("blue", alpha.f = 0.75), lwd=4)
              }
              
              if (cline=="UDL and NL") {
                if (width_method=="tangent") {
                  # add fit to plot
                  lines(c(NL.X.high, NL.X.low), c(1, 0), col=adjustcolor("grey70", alpha.f = 0.5), lwd=4)
                  lines(NL.x.spaced, NL.y.predicted, col=adjustcolor("grey30", alpha.f = 0.75), lwd=4)
                }
                
                if (width_method=="eightieth") {
                  # add fit to plot
                  polygon(x=c(NL.left_width_margin, NL.right_width_margin, NL.right_width_margin, NL.left_width_margin), y=c(-0.05, -0.05, 1.05, 1.05), col=adjustcolor("grey70", alpha.f = 0.5), border=NA)
                  lines(NL.x.spaced, NL.y.predicted, col=adjustcolor("grey30", alpha.f = 0.75), lwd=4)
                }
              }
            }
            
            if (plot_NL) {
              # graph hybrid index based on NL loci
              x_M <- pop_matrix_M[,1]
              y_M <- rowMeans(pop_matrix_M[,NL_col_start:NL_col_end, drop=F])  # the "drop=F" prevents an error when only one individual
              colour <- adjustcolor(NL_colour, alpha.f = 0.5)
              points(x_M, jitter(y_M, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
              x_F <- pop_matrix_F[,1]
              y_F <- rowMeans(pop_matrix_F[,NL_col_start:NL_col_end, drop=F])
              colour <- adjustcolor(NL_colour, alpha.f = 0.5)
              points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
            }
            
            if (plot_UDL) {
              # graph hybrid index based on UDL loci
              x_M <- pop_matrix_M[,1]
              y_M <- rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end, drop=F])
              colour <- adjustcolor(UDL_colour, alpha.f = 0.5)
              points(x_M, jitter(y_M, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
              x_F <- pop_matrix_F[,1]
              y_F <- rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end, drop=F])
              colour <- adjustcolor(UDL_colour, alpha.f = 0.5)
              points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
            }
            
            if (plot_MTL) {
              # graph hybrid index based on MTL loci
              x_M <- pop_matrix_M[,1]
              y_M <- rowMeans(pop_matrix_M[,MTL_col_start:MTL_col_end])
              colour <- adjustcolor(MTL_colour, alpha.f = 0.5)
              points(x_M, jitter(y_M, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
              x_F <- pop_matrix_F[,1]
              y_F <- rowMeans(pop_matrix_F[,MTL_col_start:MTL_col_end])
              colour <- adjustcolor(MTL_colour, alpha.f = 0.5)
              points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
            }
            
            mtext(text=paste0("generation ",time), side=3, line=0.5, adj=0)
            
            
            if (plot_density) {
              # graph density
              par(mar=c(1,3,1,1))  # specifies number of lines around plot (bottom, left, top right)
              F_hist <- hist(x_F, breaks=seq(0, 1, by=1/demes), plot=FALSE)
              M_hist <- hist(x_M, breaks=seq(0, 1, by=1/demes), plot=FALSE)
              barplot(F_hist$counts + M_hist$counts, space=0, cex.axis=0.8, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
              title(ylab="Density ", line=2, cex.lab=1.2)
            }
            
            if (make_movie & (time %in% movie_gen)) {
              main_dev_number <- dev.cur()
              digit <- NULL
              if (time < 10) {
                digit <- "000"
              } else if (time < 100) {
                digit <- "00"
              } else if (time < 1000) {
                digit <- "0"
              }
              dev.copy(png, width = 6, height = 5, units = "in", res=200, filename=paste0("new_movie/",movie_name,digit,time,".png"))  # copies plot to png file
              dev.off()
              dev.set(main_dev_number)  # returns to main quartz screen
            }
          }
        }           
      }    
      
      #### The above finishes the one simulation.
      
      # Record results of the one simulation
      
      if (extinction==TRUE) {   # if whole population went extinct
        outcome <- "extinction"
        final_distribution <- NULL
        species0_proportion <- NULL
        species1_proportion <- NULL
        species0_proportion_NL <- NULL
        species1_proportion_NL <- NULL
        if (track_fitness) {
          save(outcome,
               mother_count_per_HI_by_gen, offspring_count_per_motherHI_by_gen,
               father_count_per_HI_by_gen, offspring_count_per_fatherHI_by_gen, matings_count_per_maleHI_by_gen,
               HI_category_values,
               file=paste0("simulation_data.",run_name,".RData"))
        } else {
          save(outcome, file=paste0("simulation_data.",run_name,".RData"))
        }
        
      } else if (extinction==FALSE) {   # if population did not go extinct
        # choose data to include in final output
        locus_type <- "UDL and NL"
        x = c(pop_matrix_F[,1] , pop_matrix_M[,1])
        if (locus_type=="UDL") {
          y <- c(rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end]) , rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end]))
        } else if (locus_type=="MTL") {
          y <- c(rowMeans(pop_matrix_F[,MTL_col_start:MTL_col_end]) , rowMeans(pop_matrix_M[,MTL_col_start:MTL_col_end]))
        } else if (locus_type=="UDL and NL") {
          y <- c(rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end]) , rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end]))
          y_NL <- c(rowMeans(pop_matrix_F[,NL_col_start:NL_col_end]) , rowMeans(pop_matrix_M[,NL_col_start:NL_col_end]))
        }
        
        # add calculation of outcomes over entire population over whole range
        # (started adding this on 7Jan2020)
        final_distribution <- data.frame(x, y, y_NL)
        species0_proportion <- sum(y == 0) / length(y)
        species1_proportion <- sum(y == 1) / length(y)
        if (species0_proportion >= 0.85 | species1_proportion >= 0.85) {
          outcome <- "one_species"
        } else if ((species0_proportion + species1_proportion > 0.85) & (species0_proportion >= 0.15) & (species1_proportion >= 0.15)) {
          outcome <- "two_species"
        } else {
          outcome <- "blended"
        } 
        species0_proportion_NL <- sum(y_NL == 0) / length(y_NL)
        species1_proportion_NL <- sum(y_NL == 1) / length(y_NL)
        
        if (track_fitness) {
          save(outcome, final_distribution, species0_proportion, species1_proportion, 
               species0_proportion_NL, species1_proportion_NL, 
               mother_count_per_HI_by_gen, offspring_count_per_motherHI_by_gen,
               father_count_per_HI_by_gen, offspring_count_per_fatherHI_by_gen, matings_count_per_maleHI_by_gen,
               HI_category_values,
               file=paste0("simulation_data.",run_name,".RData"))
        } else {
          save(outcome, final_distribution, species0_proportion, species1_proportion, species0_proportion_NL, species1_proportion_NL, file=paste0("simulation_data.",run_name,".RData"))
        }
      } 
      
      # Below is the final cline fitting and plotting:
      # Fit to cline (not used in HZAM-sym paper; not meaningful when pure sympatry)
      if (fit_to_cline) {
        # choose the data to include
        
        x = c(pop_matrix_F[,1] , pop_matrix_M[,1])
        if (cline=="UDL") {
          y <- c(rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end]) , rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end]))
        }
        else if (cline=="MTL") {
          y <- c(rowMeans(pop_matrix_F[,MTL_col_start:MTL_col_end]) , rowMeans(pop_matrix_M[,MTL_col_start:MTL_col_end]))
        }
        else if (cline=="UDL and NL") {
          y <- c(rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end]) , rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end]))
          y_NL <- c(rowMeans(pop_matrix_F[,NL_col_start:NL_col_end]) , rowMeans(pop_matrix_M[,NL_col_start:NL_col_end]))
        }
        
        # fit UDL or MTL to cline
        
        mydata <- as.data.frame(cbind(x,y))
        
        z <- gam(y ~ s(x), data=  mydata, quasibinomial(link = "logit"),
                 method = "ML")
        x.spaced <- seq(0, 1, by=0.001)
        y.predicted <- inv.logit(as.vector(predict.gam(z, newdata=data.frame(x=x.spaced))))
        
        if (width_method=="tangent") {
          #find the x with maximum slope
          first.diff <- diff(y.predicted)
          index <- which.max(first.diff)
          x.1 <- x.spaced[index]
          x.2 <- x.spaced[index+1]
          y.1 <- y.predicted[index]
          y.2 <- y.predicted[index+1]
          max.slope <- (y.2 - y.1) / (x.2 - x.1)
          width <- 1/max.slope   # if asymptotes at 0 and 1
          print(paste0("width = ",round(width, digits=3)))
          intercept <- y.1 - max.slope*x.1
          X.high <- (1 - intercept) / max.slope
          X.low <- (0 - intercept) / max.slope
          # add fit to plot
          lines(c(X.high, X.low), c(1, 0), col=adjustcolor("grey", alpha.f = 0.5), lwd=4)
          lines(x.spaced, y.predicted, col=adjustcolor("blue", alpha.f = 0.75), lwd=4)
        }
        
        if (width_method=="eightieth") {
          # find place that fit goes below 0.05 on left of center
          left_width_margin <- max(x.spaced[y.predicted <= 0.1])
          # find place where fit goes above 0.95 on right of center
          right_width_margin <- min(x.spaced[y.predicted >= 0.9])
          width <- right_width_margin - left_width_margin
          print(paste0("width = ",round(width, digits=3)))
          polygon(x=c(left_width_margin, right_width_margin, right_width_margin, left_width_margin), y=c(-0.05, -0.05, 1.05, 1.05), col=adjustcolor("skyblue", alpha.f = 0.5), border=NA)
          lines(x.spaced, y.predicted, col=adjustcolor("blue", alpha.f = 0.75), lwd=4)
          
          # add calculation of "bimodality" in centre of zone (where HI fit = 0.5)
          # will base on a half dispersal distance on either side of centre
          # (added this on 30March2019)
          dist_from_HI50percent <- abs(y.predicted-0.5)
          centre <- x.spaced[dist_from_HI50percent == min(dist_from_HI50percent)]
          upper_lim <- centre + meandispersal/2
          lower_lim <- centre - meandispersal/2
          selection <- (x > lower_lim) & (x < upper_lim)
          x_near_centre <- x[selection]
          y_near_centre <- y[selection]
          y_NL_near_centre <- y_NL[selection]
          bimodality <- (sum(y_near_centre==0) + sum(y_near_centre==1)) / length(y_near_centre)
          centre_distribution <- cbind(x_near_centre, y_near_centre, y_NL_near_centre)
        }
        
        # this part is for the NL clinefit and width
        if (cline=="UDL and NL") {     
          
          mydata <- as.data.frame(cbind(x,y_NL))
          
          z <- gam(y_NL ~ s(x), data=  mydata, quasibinomial(link = "logit"),
                   method = "ML")
          x.spaced <- seq(0, 1, by=0.001)
          y.predicted <- inv.logit(as.vector(predict.gam(z, newdata=data.frame(x=x.spaced))))
          #thank Dolph for his help with above                      
          
          if (width_method=="tangent") {
            #find the x with maximum slope
            first.diff <- diff(y.predicted)
            index <- which.max(first.diff)
            x.1 <- x.spaced[index]
            x.2 <- x.spaced[index+1]
            y.1 <- y.predicted[index]
            y.2 <- y.predicted[index+1]
            max.slope <- (y.2 - y.1) / (x.2 - x.1)
            width_NL <- 1/max.slope   # if asymptotes at 0 and 1
            print(paste0("width_NL = ",round(width_NL, digits=3)))
            
            intercept <- y.1 - max.slope*x.1
            
            X.high <- (1 - intercept) / max.slope
            X.low <- (0 - intercept) / max.slope
            # add fit to plot
            lines(c(X.high, X.low), c(1, 0), col=adjustcolor("grey70", alpha.f = 0.5), lwd=4)
            lines(x.spaced, y.predicted, col=adjustcolor("grey30", alpha.f = 0.75), lwd=4)
          }
          
          if (width_method=="eightieth") {
            # find place that fit goes below 0.05 on left of center
            left_width_margin <- max(x.spaced[y.predicted <= 0.1])
            # find place where fit goes above 0.95 on right of center
            right_width_margin <- min(x.spaced[y.predicted >= 0.9])
            width_NL <- right_width_margin - left_width_margin
            print(paste0("width_NL = ",round(width_NL, digits=3)))
            polygon(x=c(left_width_margin, right_width_margin, right_width_margin, left_width_margin), y=c(-0.05, -0.05, 1.05, 1.05), col=adjustcolor("grey70", alpha.f = 0.5), border=NA)
            lines(x.spaced, y.predicted, col=adjustcolor("grey30", alpha.f = 0.75), lwd=4)
          }
        }
        
        if (plot_NL) {
          # graph hybrid index based on NL loci
          x_M <- pop_matrix_M[,1]
          y_M <- rowMeans(pop_matrix_M[,NL_col_start:NL_col_end])
          colour <- adjustcolor(NL_colour, alpha.f = 0.5)
          points(x_M, jitter(y_M, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
          x_F <- pop_matrix_F[,1]
          y_F <- rowMeans(pop_matrix_F[,NL_col_start:NL_col_end])
          colour <- adjustcolor(NL_colour, alpha.f = 0.5)
          points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
        }
        
        if (plot_UDL) {
          # graph hybrid index based on UDL loci
          x_M <- pop_matrix_M[,1]
          y_M <- rowMeans(pop_matrix_M[,UDL_col_start:UDL_col_end])
          colour <- adjustcolor(UDL_colour, alpha.f = 0.5)
          points(x_M, jitter(y_M, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
          x_F <- pop_matrix_F[,1]
          y_F <- rowMeans(pop_matrix_F[,UDL_col_start:UDL_col_end])
          colour <- adjustcolor(UDL_colour, alpha.f = 0.5)
          points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
        }
        
        if (plot_MTL) {
          # graph hybrid index based on MTL loci
          x_M <- pop_matrix_M[,1]
          y_M <- rowMeans(pop_matrix_M[,MTL_col_start:MTL_col_end])
          colour <- adjustcolor(MTL_colour, alpha.f = 0.5)
          points(x_M, jitter(y_M, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
          x_F <- pop_matrix_F[,1]
          y_F <- rowMeans(pop_matrix_F[,MTL_col_start:MTL_col_end])
          colour <- adjustcolor(MTL_colour, alpha.f = 0.5)
          points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
        }
        
        # save image
        main_dev_number <- dev.cur()
        digit <- NULL
        if (time < 10) {
          digit <- "000"
        } else if (time < 100) {
          digit <- "00"
        } else if (time < 1000) {
          digit <- "0"
        }
        dev.copy(png, width = 6, height = 5, units = "in", res=200, filename=paste0("HZAM_animations/",run_name,".png"))  # copies plot to png file
        dev.off()
        dev.set(main_dev_number)  # returns to main quartz screen
      }
      
      # Make histogram of final HI distribution
      if (plot_histogram_at_end & outcome != "extinction") {
        quartz(width=5, height=3, title=paste0("HZAM: Histogram of HI distribution; pref_ratio=", pref_ratio,"; hybrid_fitness=", hybrid_fitness, sep=""), bg="white")
        # graph hybrid index based on MTL loci
        HI_F <- rowMeans(pop_matrix_F[,MTL_col_start:MTL_col_end])
        HI_M <- rowMeans(pop_matrix_M[,MTL_col_start:MTL_col_end])
        colour <- adjustcolor(MTL_colour, alpha.f = 0.5)
        # points(x_F, jitter(y_F, amount=0.015), col=colour, bg=colour, pch=16, cex=0.6)
        HI_categories <- 2*male_trait_loci + 1
        HI_category_values_for_hist <- seq(0, 1, by=1/(HI_categories-1))
        HI_category_breaks <- rep(NA, HI_categories-1)  # create empty vector
        for (i in 1:(HI_categories-1)) {
          HI_category_breaks[i] <- mean(HI_category_values_for_hist[i:(i+1)]) 
        }
        category_labels <- round(HI_category_values_for_hist, 2)
        F_hist <- hist(HI_F, breaks=c(-0.1, HI_category_breaks, 1.1), plot=FALSE)
        M_hist <- hist(HI_M, breaks=c(-0.1, HI_category_breaks, 1.1), plot=FALSE)
        barplot(F_hist$counts + M_hist$counts, space=0.1, names.arg=category_labels, col=colour, cex.axis=0.8, tcl=-0.25, xlab=NA, ylab=NA, mgp=c(3,0.5,0), las=1)
        title(ylab="Individuals", line=2, cex.lab=1)
        title(xlab="Hybrid Index (HI)", line=2, cex.lab=1)
        title(main="Distribution of HI at end of simulation", cex.lab=1)
        title(sub=paste0("pref_ratio=", round(pref_ratio,4),"; hybrid_fitness=", hybrid_fitness, sep=""), cex.lab=0.8)
        
        # save image
        main_dev_number <- dev.cur()
        dev.copy(png, width = 5, height = 3, units = "in", res=300, filename=paste0("final_histogram",run_name,".png"))  # copies plot to png file
        dev.off()
        dev.set(main_dev_number)  # returns to main quartz screen
        if (close_histogram) {
          dev.off()  # closes the quartz window with the histogram
        }
      }
    }
  }
}


# Code for making movie ----

# animation <- image_animate(img, fps = 1, loop = 0)
# print(animation)
# image_write(animation, path=paste0(movie_name,".gif"), format="gif")

#run ImageMagick in a terminal window
# convert *.png -delay 10 -loop 0 animation.gif

# converted this online into mp4 format, then brought into iMovie and slowed down (10% speed), then exported as MP4
# or use ezgif.com to change speed



#### Part 2: Figures in the manuscript ----

# for the below, the data files should be in a directory called this:
data_folder <- "simulation_data_HZAM_files"

#### Graphs of population over generations, by HI ----
# For Figs. 2 and 5

hybrid_fitness_set <- NULL
mate_pref_set <- NULL
run_number_set <- NULL
run_set_ID <- NULL
gens_to_show <- list(NULL)
make_graph_toggle <- NULL
plot_offspring_per_HI_toggle <- NULL
gen_window_set <- NULL

i <- 1   # Fig. 2A: no ecological diff, complete premating isolation, no postzygotic isolation --> one species (stochastic extinction)
hybrid_fitness_set[i] <- 1
mate_pref_set[i] <- 0
run_number_set[i] <- 1
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.05_ecolDiff0_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1,200,400,600,800,1000)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- F
gen_window_set[i] <- NULL

i <- 2    # Fig. 2B: yes ecol diff, complete premating isolation, no postzygotic isolation --> maintenance of two species
hybrid_fitness_set[i] <- 1
mate_pref_set[i] <- 0
run_number_set[i] <- 1
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1,200,400,600,800,1000)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- F
gen_window_set[i] <- NULL

i <- 3    # Fig. 2C: yes ecol diff, 10x assortative mating, no postzygotic isolation --> blending
hybrid_fitness_set[i] <- 1
mate_pref_set[i] <- 0.1
run_number_set[i] <- 1
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1,100,200,300,400)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- F
gen_window_set[i] <- NULL

i <- 4    # Fig. 2D: yes ecol diff, 10x assortative mating, 0.6 hybrid fitness --> blending
hybrid_fitness_set[i] <- 0.6
mate_pref_set[i] <- 0.1
run_number_set[i] <- 3
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1,50,100)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- F
gen_window_set[i] <- NULL

i <- 5    # Fig. 2E: yes ecol diff, 3X assortative mating, 0.6 hybrid fitness --> blending
hybrid_fitness_set[i] <- 0.6
mate_pref_set[i] <- 0.333333
run_number_set[i] <- 1
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1,10,20,30,40,50)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- F
gen_window_set[i] <- NULL

i <- 6   # Fig. 5A
hybrid_fitness_set[i] <- 1
mate_pref_set[i] <-  0.001
run_number_set[i] <- 1
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.2_ecolDiff1_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1, 200, 400, 600, 800, 1000)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- T
gen_window_set[i] <- 150

i <- 7   # Fig. 5B
hybrid_fitness_set[i] <- 1
mate_pref_set[i] <-  0.00333333
run_number_set[i] <- 5
run_set_ID[i] <- paste0("simulation_data.HZAM_animation_runA_rep",run_number_set[i],"_growthrate1.2_ecolDiff1_K500_UDLMTL3_gen1000")
gens_to_show[[i]] <- c(1, 100, 200, 300, 400)
make_graph_toggle[i] <- T
plot_offspring_per_HI_toggle[i] <- T
gen_window_set[i] <- 30

for (i in 1:length(run_set_ID)) {
  if (make_graph_toggle[i]) {
    max_gens_to_show <- max(gens_to_show[[i]])
    filename <- paste0(run_set_ID[i],"_hybridfitness",hybrid_fitness_set[i],"_prefratio",mate_pref_set[i],".RData")
    load(paste0(data_folder,"/",filename))  # load the data for a single run (contains final_distribution, bimodality, bimodality_NL)
    popN_per_HI_by_gen <- father_count_per_HI_by_gen[,1:max_gens_to_show] + mother_count_per_HI_by_gen[,1:max_gens_to_show]
    total_N_by_gen <- colSums(popN_per_HI_by_gen)
    
    total_generations <- ncol(popN_per_HI_by_gen)
    generation_values <- 1:total_generations
    HI_categories <- nrow(popN_per_HI_by_gen)
    HI_values <- c(0, (1:(HI_categories-1) / (HI_categories-1)))
    colours <- viridis(length(HI_values))
    
    if (plot_offspring_per_HI_toggle[i] == FALSE) {
      quartz(title=paste0("Graph ",i,": Population size per phenotypic category over time (hybridfitness",hybrid_fitness_set[i],"; prefratio",mate_pref_set[i],")"), width=5, height=4, bg="white")
      par(oma=c(3,3,3,1))  # set outer margins
      zones <- matrix(c(1,2), ncol=1, byrow=TRUE)  # numbers in matrix give order of plotting
      layout(zones, widths=c(1,1), heights=c(2/8,3/8))
      # in the first plotting area:
      par(mar=c(0,3,3,1))  # specifies number of father_count_per_HI_by_genlines around plot (bottom, left, top, right)
      plot(x=1:length(total_N_by_gen), y=total_N_by_gen, type="l", xlim = c(1, total_generations), ylim = c(0,1400), xaxs = "i", yaxs = "i",
           main="Outcome of contact of two species", cex.main=1.2, xaxt = "n", yaxt = "n")
      axis(1, at=gens_to_show[[i]], labels=F, tcl = -0.25)
      axis(2, at=c(0,500,1000), labels=T, cex.axis = 0.75, las=1, tcl = -0.25, mgp = c(3,0.5,0))
      # in the second plotting area:
      par(mar=c(3,3,1,1))  # specifies number of lines around plot (bottom, left, top, right)
      image(x=generation_values, y=HI_values, z= t(popN_per_HI_by_gen),
            col = gray.colors(1000, start = 0, end = 1, gamma = 0.1, alpha = NULL, rev = T),
            xaxt = "n", yaxt = "n")
      axis(1, at=gens_to_show[[i]], labels=T, cex.axis = 0.75, tcl = -0.25, mgp = c(3,0.5,0))
      axis(2, at=HI_values, labels=format(HI_values, digits=2), cex.axis = 0.75, line=0, las=1, tcl = -0.25, mgp = c(3,0.5,0))
      mtext("Generations", side=1, line=-1, outer=TRUE, cex=1.5,
            at=0.55)
      mtext("HI", side=2, line = -0.5, outer=TRUE, cex=1.5, las=1, 
            at=0.37)
      mtext("N", side=2, line = -0.5, outer=TRUE, cex=1.5, las=1, 
            at=0.68)
    } 
    
    if (plot_offspring_per_HI_toggle[i] == TRUE) {
      axis_label_size <- 1.1
      quartz(title=paste0("Graph ",i,": Population size and fitness per phenotypic category over time (hybridfitness",hybrid_fitness_set[i],"; prefratio",mate_pref_set[i],")"), width=5, height=5, bg="white")
      par(oma=c(3,6,3,1))  # set outer margins
      zones <- matrix(c(1,2,3), ncol=1, byrow=TRUE)  # numbers in matrix give order of plotting
      layout(zones, widths=c(1,1,1), heights=c(2/10,2/10,5/10))
      # in the first plotting area:
      par(mar=c(0,3,3,1))  # specifies number of lines around plot (bottom, left, top, right)
      plot(x=1:length(total_N_by_gen), y=total_N_by_gen, type="l", xlim = c(1, total_generations), ylim = c(0,1400), xaxs = "i", yaxs = "i",
           main="Outcome of contact of two species", cex.main=1.4, xaxt = "n", yaxt = "n")
      axis(1, at=gens_to_show[[i]], labels=F, tcl = -0.25)
      axis(2, at=c(0,500,1000), labels=T, cex.axis = axis_label_size, las=1, tcl = -0.25, mgp = c(3,0.5,0))
      # in the second plotting area:
      par(mar=c(0,3,1,1))  # specifies number of lines around plot (bottom, left, top, right)
      image(x=generation_values, y=HI_values, z= t(popN_per_HI_by_gen),
            col = gray.colors(1000, start = 0, end = 1, gamma = 0.1, alpha = NULL, rev = T),
            xaxt = "n", yaxt = "n")
      axis(1, at=gens_to_show[[i]], labels=F, cex.axis = axis_label_size, tcl = -0.25, mgp = c(3,0.5,0))
      axis(2, at=HI_values, labels=format(HI_values, digits=2), cex.axis = axis_label_size, line=0, las=1, tcl = -0.25, mgp = c(3,0.5,0))
      
      # in the third plotting area:
      par(mar=c(3,3,1,1))  # specifies number of lines around plot (bottom, left, top, right)
      # rolling mean graph for total reproductive fitness (mothers + fathers)
      # graph for mothers
      mean_offspring_by_motherHI_by_gen <- offspring_count_per_motherHI_by_gen / mother_count_per_HI_by_gen
      mean_offspring_by_fatherHI_by_gen <- offspring_count_per_fatherHI_by_gen / father_count_per_HI_by_gen
      mean_offspring_by_parentHI_by_gen <- ((offspring_count_per_motherHI_by_gen + offspring_count_per_fatherHI_by_gen)/2) / (mother_count_per_HI_by_gen + father_count_per_HI_by_gen)
      gen_window <- gen_window_set[i]
      times_with_HI_category1 <- which(!is.na(mean_offspring_by_parentHI_by_gen[1,]))
      offspring_at_times_with_HI_category1 <- mean_offspring_by_parentHI_by_gen[1, !is.na(mean_offspring_by_parentHI_by_gen[1,])]
      windowed_times_category1 <- rollapply(times_with_HI_category1, width=gen_window, FUN=mean, by=gen_window)
      windowed_mean_offspring_category1 <- rollapply(offspring_at_times_with_HI_category1, width=gen_window, FUN=mean, by=gen_window) 
      plot(x=windowed_times_category1, xlim = c(1, total_generations), xaxs = "i", xaxt = "n", 
           y=windowed_mean_offspring_category1, ylim=c(0.6, 1.2), yaxs = "i", yaxt = "n",
           type="l", col=colours[1])
      for (j in 2:length(HI_category_values)) {
        times_with_HI_category <- which(!is.na(mean_offspring_by_parentHI_by_gen[j,]))
        offspring_at_times_with_HI_category <- mean_offspring_by_parentHI_by_gen[j, !is.na(mean_offspring_by_parentHI_by_gen[j,])]
        windowed_times_category <- rollapply(times_with_HI_category, width=gen_window, FUN=mean, by=gen_window)
        windowed_mean_offspring_category <- rollapply(offspring_at_times_with_HI_category, width=gen_window, FUN=mean, by=gen_window) 
        points(x=windowed_times_category, 
               y=windowed_mean_offspring_category,
               type="l", col=colours[j])
      } 
      axis(1, at=gens_to_show[[i]], labels=T, tcl = -0.25, cex.axis = axis_label_size)
      axis(2, at=c(0.7,0.8,0.9,1,1.1), labels=T, cex.axis = axis_label_size, las=1, tcl = -0.25, mgp = c(3,0.5,0))
      mtext("Generations", side=1, line=0.5, outer=TRUE, cex=1.5,
            at=0.53)
      mtext("N", side=2, line = 0.25, outer=TRUE, cex=1.5, las=1, 
            at=0.84)
      mtext("HI", side=2, line = 0.25, outer=TRUE, cex=1.5, las=1, 
            at=0.65)
      mtext("fitness", side=2, line = 0.25, outer=TRUE, cex=1.5, las=1, 
            at=0.3)
    }
  }
}



#### Graphs of summary of set of outcomes ----
# For Figs. 3, 4, 6, S1

hybrid_fitness_set <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 1)   # for just one run, just put one number in this and next line
# mate_pref_set <- c(0.999999, 0.333333, 0.1, 0.0333333, 0.01, 0.00333333, 0.001, 0)   # ratio in mating pref between hetero- and homospecific
mate_pref_set <- c(1, 0.333333, 0.1, 0.0333333, 0.01, 0.00333333, 0.001, 0)   # ratio in mating pref between hetero- and homospecific
AM_strengths <- c(1, 3, 10, 30, 100, 300, 1000, 5000)   # this is the inverse of the above

# hybrid_fitness_set <- c(0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.98, 1)
# mate_pref_set <- c(0.999999, 0.333333, 0.1, 0.01, 0.001)
# AM_strengths <- c(1, 3, 10, 100, 1000)   # this is the inverse of the above
log10_AM_strengths <- log10(AM_strengths)
AM_axis_labels <- c(paste0(AM_strengths[1:length(AM_strengths)-1], "x"), "complete")

run_set_numbers <- c(1:5)     # c(1:5)

run_set_IDs <- NULL

# to run the code for the appropriate figure, uncomment just the one line below for the appropriate figure as indicated below, at the right of each line:
for (i in run_set_numbers) {
  # Fig. 3:
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff0_K500_UDLMTL3_gen1000")  # Figure 3A
  run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")   # Figure 3B
  
  # Fig. 4:
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.025_ecolDiff1_K500_UDLMTL3_gen1000")   # Figure 4A
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.2_ecolDiff1_K500_UDLMTL3_gen1000")    # Figure 4B
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate2.6_ecolDiff1_K500_UDLMTL3_gen1000")    # Figure 4C
  
  # Fig. 6:
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff0_K500_UDLMTL3_gen1000")  # Figure 6A
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff0.25_K500_UDLMTL3_gen1000")   # Figure 6B
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff0.5_K500_UDLMTL3_gen1000")   # Figure 6C
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff0.75_K500_UDLMTL3_gen1000")   # Figure 6D
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")   # Figure 6E
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate2.6_ecolDiff0_K500_UDLMTL3_gen1000")    # Figure 6F
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate2.6_ecolDiff0.25_K500_UDLMTL3_gen1000")  # Figure 6G
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate2.6_ecolDiff0.5_K500_UDLMTL3_gen1000")   # Figure 6H
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate2.6_ecolDiff0.75_K500_UDLMTL3_gen1000")   # Figure 6I
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_rep",i,"_growthrate2.6_ecolDiff1_K500_UDLMTL3_gen1000")  # Figure 6J
  
  # Fig. S1: 
  # run_set_IDs[i] <- paste0("simulation_data.HZAM_animation_runA_noSexSel_rep",i,"_growthrate1.05_ecolDiff1_K500_UDLMTL3_gen1000")   # Figure S1
  }

# load data files and determine most common outcome for each condition:
outcome_array <- matrix(nrow = length(mate_pref_set), ncol = length(hybrid_fitness_set))  # pre-allocate a matrix with NA entries
NL_outcome_array <- matrix(nrow = length(mate_pref_set), ncol = length(hybrid_fitness_set))  # pre-allocate a matrix with NA entries
for (i in 1:length(mate_pref_set)) {
  for (j in 1:length(hybrid_fitness_set)) {
    one_outcome <- NULL
    one_NL_outcome <- NULL
    for (which_run in 1:length(run_set_IDs)) {
      filename <- paste0(data_folder,"/",run_set_IDs[which_run],"_hybridfitness",hybrid_fitness_set[j],"_prefratio",mate_pref_set[i],".RData")
      print(paste0("loading ",filename))
      load(filename)  # load the data for a single run (contains final_distribution, bimodality, bimodality_NL)
      if (outcome=="extinction") {
        one_outcome[which_run] <- "extinction"  # integer code 1
        one_NL_outcome[which_run] <- "extinction"  # integer code 1
      } else {
        species0_proportion <- sum(final_distribution$y == 0) / length(final_distribution$y)
        species1_proportion <- sum(final_distribution$y == 1) / length(final_distribution$y)
        if (species0_proportion >= 0.85 | species1_proportion >= 0.85) {
          one_outcome[which_run] <- "one_species"  # integer code 2
        } else if ((species0_proportion + species1_proportion > 0.85) & (species0_proportion >= 0.15) & (species1_proportion >= 0.15)) {
          one_outcome[which_run] <- "two_species"  # integer code 3
        } else {
          one_outcome[which_run] <- "blended"  # integer code 4
        } 
      }
    }
    # Choose the most common outcome.
    most_common_outcomes <- table(one_outcome)[table(one_outcome) == max(table(one_outcome))]
    if (length(most_common_outcomes) == 1) {
      common_outcome <- names(most_common_outcomes)
    }
    # If there is a tie for most common outcome, this chooses randomly from the most common outcome:
    if (length(most_common_outcomes) >= 2) {
      common_outcome <- names(sample(table(one_outcome)[table(one_outcome) == max(table(one_outcome))],1) ) 
    }
    outcome_array[i,j] <- common_outcome
  } 
}

# print out array of outcomes:
outcome_array

# This is the color combo we like best:
outcome_key <- rbind(c("extinction", colorBlindBlackPlasma41[1]),  # integer code 1
                     c("one_species", colorBlindBlackPlasma41[9]),  # integer code 2, etc.
                     c("blended", colorBlindBlackPlasma41[22]),     # 3
                     c("two_species", colorBlindBlackPlasma41[37]))   # 4

#convert to numeric:
outcome_array_numeric <- matrix(nrow=nrow(outcome_array), ncol=ncol(outcome_array))
for (row in 1:nrow(outcome_key)) {
  outcome_array_numeric[outcome_array==outcome_key[row,1]] <- row
}

# graph outcomes:
quartz(title="outcomes in response to mate pref and hybrid fitness width", width=6, height=5, bg="white")
min_category <- min(outcome_array_numeric)  # this and next line needed to make colour key work right
max_category <- max(outcome_array_numeric)
image(x=log10_AM_strengths, y=hybrid_fitness_set, z=outcome_array_numeric, col=outcome_key[min_category:max_category, 2],
      xlab="Strength of conspecific mate preference",
      ylab="Hybrid fitness", cex.lab=1.25,
      main="Outcomes of contact of two species", cex.main=1.5,
      xaxt="n", yaxt="n")
axis(1, at=log10_AM_strengths, labels=AM_axis_labels)
axis(2, at=hybrid_fitness_set, labels=format(hybrid_fitness_set, digits=3), las=1)
# add white lines to separate cases of complete RI from the rest:
x_line <- mean(log10_AM_strengths[(length(log10_AM_strengths)-1):length(log10_AM_strengths)])
lines(x=c(x_line, x_line), y=c(-0.1,1.1), lwd=5, col="white")
y_line <- mean(hybrid_fitness_set[1:2])
lines(x=c(-10, 10), y=c(y_line,y_line), lwd=5, col="white")



#### Part 3: Code for generating other possibly useful figures ----


# Graphs of offspring per HI class over time ----

colours <- viridis(length(HI_category_values))

# graph for mothers
mean_offspring_by_motherHI_by_gen <- offspring_count_per_motherHI_by_gen / mother_count_per_HI_by_gen
quartz(width=5, height=4, title="offspring by mother HI")
times_with_HI_category1 <- which(!is.na(mean_offspring_by_motherHI_by_gen[1,]))
offspring_at_times_with_HI_category1 <- mean_offspring_by_motherHI_by_gen[1, !is.na(mean_offspring_by_motherHI_by_gen[1,])]
plot(x=times_with_HI_category1, xlab="generation time",
     y=offspring_at_times_with_HI_category1, ylab="mean offspring per mother", ylim=c(0, 5),
     type="l", col=colours[1])
for (i in 2:length(HI_category_values)) {
  times_with_HI_category <- which(!is.na(mean_offspring_by_motherHI_by_gen[i,]))
  offspring_at_times_with_HI_category <- mean_offspring_by_motherHI_by_gen[i, !is.na(mean_offspring_by_motherHI_by_gen[i,])]
  points(x=times_with_HI_category, 
         y=offspring_at_times_with_HI_category,
         type="l", col=colours[i])
}

# graph for fathers
mean_offspring_by_fatherHI_by_gen <- offspring_count_per_fatherHI_by_gen / father_count_per_HI_by_gen
quartz(width=5, height=4, title="offspring by father HI")
times_with_HI_category1 <- which(!is.na(mean_offspring_by_fatherHI_by_gen[1,]))
offspring_at_times_with_HI_category1 <- mean_offspring_by_fatherHI_by_gen[1, !is.na(mean_offspring_by_fatherHI_by_gen[1,])]
plot(x=times_with_HI_category1, xlab="generation time",
     y=offspring_at_times_with_HI_category1, ylab="mean offspring per father", ylim=c(0, 5),
     type="l", col=colours[1])
for (i in 2:length(HI_category_values)) {
  times_with_HI_category <- which(!is.na(mean_offspring_by_fatherHI_by_gen[i,]))
  offspring_at_times_with_HI_category <- mean_offspring_by_fatherHI_by_gen[i, !is.na(mean_offspring_by_fatherHI_by_gen[i,])]
  points(x=times_with_HI_category, 
         y=offspring_at_times_with_HI_category,
         type="l", col=colours[i])
} 

# graph for total reproductive fitness (mothers + fathers)
mean_offspring_by_parentHI_by_gen <- (offspring_count_per_motherHI_by_gen + offspring_count_per_fatherHI_by_gen) / (mother_count_per_HI_by_gen + father_count_per_HI_by_gen)
quartz(width=5, height=4, title="offspring by parent HI")
times_with_HI_category1 <- which(!is.na(mean_offspring_by_parentHI_by_gen[1,]))
offspring_at_times_with_HI_category1 <- mean_offspring_by_parentHI_by_gen[1, !is.na(mean_offspring_by_parentHI_by_gen[1,])]
plot(x=times_with_HI_category1, xlab="generation time",
     y=offspring_at_times_with_HI_category1, ylab="mean offspring per parent", ylim=c(0, 5),
     type="l", col=colours[1])
for (i in 2:length(HI_category_values)) {
  times_with_HI_category <- which(!is.na(mean_offspring_by_parentHI_by_gen[i,]))
  offspring_at_times_with_HI_category <- mean_offspring_by_parentHI_by_gen[i, !is.na(mean_offspring_by_parentHI_by_gen[i,])]
  points(x=times_with_HI_category, 
         y=offspring_at_times_with_HI_category,
         type="l", col=colours[i])
} 

# rolling mean graph for total reproductive fitness (mothers + fathers)
gen_window <- 10
quartz(width=5, height=4, title="rolling means, offspring by parent HI")
times_with_HI_category1 <- which(!is.na(mean_offspring_by_parentHI_by_gen[1,]))
offspring_at_times_with_HI_category1 <- mean_offspring_by_parentHI_by_gen[1, !is.na(mean_offspring_by_parentHI_by_gen[1,])]
windowed_times_category1 <- rollapply(times_with_HI_category1, width=gen_window, FUN=mean, by=gen_window)
windowed_mean_offspring_category1 <- rollapply(offspring_at_times_with_HI_category1, width=gen_window, FUN=mean, by=gen_window) 
plot(x=windowed_times_category1, xlab="generation time",
     y=windowed_mean_offspring_category1, ylab=paste0("mean offspring per parent (window=", gen_window,")"), ylim=c(0, 3),
     type="l", col=colours[1])
for (i in 2:length(HI_category_values)) {
  times_with_HI_category <- which(!is.na(mean_offspring_by_parentHI_by_gen[i,]))
  offspring_at_times_with_HI_category <- mean_offspring_by_parentHI_by_gen[i, !is.na(mean_offspring_by_parentHI_by_gen[i,])]
  windowed_times_category <- rollapply(times_with_HI_category, width=gen_window, FUN=mean, by=gen_window)
  windowed_mean_offspring_category <- rollapply(offspring_at_times_with_HI_category, width=gen_window, FUN=mean, by=gen_window) 
  points(x=windowed_times_category, 
         y=windowed_mean_offspring_category,
         type="l", col=colours[i])
} 

# graph for male matings by HI
mean_mating_by_maleHI_by_gen <- matings_count_per_maleHI_by_gen / father_count_per_HI_by_gen
quartz(width=5, height=4, title="matings by male HI")
times_with_HI_category1 <- which(!is.na(mean_mating_by_maleHI_by_gen[1,]))
matings_at_times_with_HI_category1 <- mean_mating_by_maleHI_by_gen[1, !is.na(mean_mating_by_maleHI_by_gen[1,])]
plot(x=times_with_HI_category1, xlab="generation time",
     y=matings_at_times_with_HI_category1, ylab="mean matings per male", ylim=c(0, 2),
     type="l", col=colours[1])
for (i in 2:length(HI_category_values)) {
  times_with_HI_category <- which(!is.na(mean_mating_by_maleHI_by_gen[i,]))
  matings_at_times_with_HI_category <- mean_mating_by_maleHI_by_gen[i, !is.na(mean_mating_by_maleHI_by_gen[i,])]
  points(x=times_with_HI_category, 
         y=matings_at_times_with_HI_category,
         type="l", col=colours[i])
} 

# rolling mean graph for male matings by HI
gen_window <- 10
mean_mating_by_maleHI_by_gen <- matings_count_per_maleHI_by_gen / father_count_per_HI_by_gen
quartz(width=5, height=4, title="matings by male HI")
times_with_HI_category1 <- which(!is.na(mean_mating_by_maleHI_by_gen[1,]))
matings_at_times_with_HI_category1 <- mean_mating_by_maleHI_by_gen[1, !is.na(mean_mating_by_maleHI_by_gen[1,])]
windowed_times_category1 <- rollapply(times_with_HI_category1, width=gen_window, FUN=mean, by=gen_window)
windowed_mean_matings_category1 <- rollapply(matings_at_times_with_HI_category1, width=gen_window, FUN=mean, by=gen_window) 
plot(x=windowed_times_category1, xlab="generation time",
     y=windowed_mean_matings_category1, ylab="mean matings per male", ylim=c(0, 2),
     type="l", col=colours[1])
for (i in 2:length(HI_category_values)) {
  times_with_HI_category <- which(!is.na(mean_mating_by_maleHI_by_gen[i,]))
  matings_at_times_with_HI_category <- mean_mating_by_maleHI_by_gen[i, !is.na(mean_mating_by_maleHI_by_gen[i,])]
  windowed_times_category <- rollapply(times_with_HI_category, width=gen_window, FUN=mean, by=gen_window)
  windowed_mean_matings_category <- rollapply(matings_at_times_with_HI_category, width=gen_window, FUN=mean, by=gen_window) 
  points(x=windowed_times_category, 
         y=windowed_mean_matings_category,
         type="l", col=colours[i])
} 

#### End of file ----

