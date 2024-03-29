################################### 
####### function definitions:
################################### 

# get the neighbourhood of a site:
nhood <- function(x, y, type = "von Neumann") {
  l1 <- list(c(x - 1, y), c(x + 1, y), c(x, y - 1), c(x, y + 1))
  if(type == "von Neumann") l2 <- list()
  else if(type == "Moore") l2 <- list(c(x - 1, y - 1), c(x + 1, y + 1), c(x - 1, y + 1), c(x + 1, y - 1))
  else stop("Neighbourhood type must be either Moore or von Neumann")
  return(c(l1, l2))
}

# count the empty spaces in a neighbourhood:
count_spaces <- function(candidate, nhood_type) {
  neighbours <- nhood(candidate[1], candidate[2], nhood_type)
  return(nhood_size - sum(sign(sapply(neighbours, function(e) sites[e[1], e[2]]))))
}

################################### 
####### parametrisation:
###################################

grid_width <- 100 # width of the grid
nhood_type <- "von Neumann" # neighbourhood type
set.seed(5) # seed for random number generator (for replicable results)
init_diameter <- 1 # population diameter
max_iter <- 1E4 # set a max number of iterations (just in case)
mutation_rate <- 0.1 # mutation rate per cell division
mean_mutation_effect <- 0.1 # mean fitness effect per mutation

################################### 
####### intialisation:
###################################

# initialise the site states:
sites <- matrix(0, nrow = grid_width, ncol = grid_width)
centre <- (init_diameter + 1)/2
occupied <- list()
for(i in 1:init_diameter) for(j in 1:init_diameter) if((i - centre)^2 + (j - centre)^2 < (init_diameter/2)^2) {
  occupied <- c(occupied, list(c(floor(grid_width/2 - init_diameter/2 + i), floor(grid_width/2 - init_diameter/2 + j))))
}
num_occupied <- length(occupied) # number of occupied sites
for(i in 1:num_occupied) sites[occupied[[i]][1], occupied[[i]][2]] <- 1

# neighbourhood size:
nhood_size <- length(nhood(1, 1, nhood_type))

# initialise the array that maps from grid sites to list entries:
index_map <- matrix(NA, nrow = grid_width, ncol = grid_width)

# initial list of occupied sites that are next to at least one empty site,
# and initial list of how many empty neighbours each has:
num_has_space <- 0
has_space <- list()
how_many_spaces <- list()
for(i in 1:num_occupied) {
  candidate <- occupied[[i]]
  spaces <- count_spaces(candidate, nhood_type)
  # if not all neighbours are occupied then add the site to the list:
  if(sites[candidate[1], candidate[2]] > 0 && spaces > 0) {
    num_has_space <- num_has_space + 1 # number of sites that have space
    has_space[[num_has_space]] <- candidate # list of sites that have space
    how_many_spaces[[num_has_space]] <- spaces # list of numbers of spaces for the sites that have space
    index_map[candidate[1], candidate[2]] <- num_has_space # map from from grid sites to list entries
  }
}

# initialise the timer and the output dataframe:
timer <- 0
output_df <- data.frame(Time = timer, Population = num_occupied)

################################### 
####### run the model:
################################### 

for(iter in 1:max_iter) {
  # pick a cell to divide, where probability of picking a cell is proportional to its number of neighbouring empty sites:
  candidate <- sample(1:num_has_space, 1, prob = unlist(how_many_spaces) * sapply(has_space, function(e) sites[e[1], e[2]]))
  # find the cell's neighbours:
  neighbours <- nhood(has_space[[candidate]][1], has_space[[candidate]][2], nhood_type)
  # randomise the order in which neighbours will be tested:
  test_order <- sample(1:length(neighbours))
  
  # test neighbours one by one:
  num_tested <- 0
  for(i in test_order) {
    # keep count of the number of occupied neighbours (for error checking):
    if(sites[neighbours[[i]][1], neighbours[[i]][2]] > 0) num_tested <- num_tested + 1
    # if the neighbouring site is empty then make it occupied:
    else {
      num_occupied <- num_occupied + 1
      new_site <- neighbours[[i]]
      parent_site <- has_space[[candidate]]
      # mutation:
      num_mutations <- rpois(1, mutation_rate) # number of new mutations drawn from a Poisson distribution
      sites[new_site[1], new_site[2]] <- sites[parent_site[1], parent_site[2]]
      if(num_mutations > 0) for(j in 1:num_mutations) {
        mutation_effect <- rexp(1, 1/mean_mutation_effect) # fitness effect per mutation drawn from an exponential distribution
        sites[new_site[1], new_site[2]] <- (1 + mutation_effect) * sites[new_site[1], new_site[2]]
      }
      break
    }
  }
  # error check:
  if(num_tested == length(neighbours)) stop("Site wrongly listed as having space to divide")
  
  # time between events is drawn from an exponential distribution,
  # consistent with the Gillespie algorithm:
  sum_of_rates <- sum(sapply(has_space, function(e) sites[e[1], e[2]]) * unlist(how_many_spaces) / nhood_size)
  timer <- timer + rexp(1, sum_of_rates)
  
  # record the time and the population size:
  output_df <- rbind(output_df, c(timer, num_occupied))
  
  # stop if a cell has reached the edge of the grid:
  if(new_site[1] >= grid_width || new_site[1] <= 1 || new_site[2] >= grid_width || new_site[2] <= 1) break
  
  # update the other lists for the newly occupied site, if it has empty neighbours:
  spaces <- count_spaces(new_site, nhood_type)
  if(spaces > 0) {
    num_has_space <- num_has_space + 1 # number of sites that have space
    has_space[[num_has_space]] <- new_site # list of sites that have space
    how_many_spaces[[num_has_space]] <- spaces # list of numbers of spaces for the sites that have space
    index_map[new_site[1], new_site[2]] <- num_has_space # map from from grid sites to list entries
  }
  
  # update list entries for neighbouring occupied sites:
  neighbours <- nhood(new_site[1], new_site[2], nhood_type)
  for(candidate in neighbours) if(sites[candidate[1], candidate[2]] > 0) {
    index <- index_map[candidate[1], candidate[2]]
    # error check:
    if(is.na(index)) stop("Site wrongly listed in index_map")
    
    # decrement the number of empty spaces neighbouring the site:
    how_many_spaces[[index]] <- how_many_spaces[[index]] - 1
    
    # if there are no neighbouring empty sites then remove the site from the lists:
    if(how_many_spaces[[index]] == 0) {
      has_space[[index]] <- NULL # removes the entry and shifts subsequent entries down one slot
      how_many_spaces[[index]] <- NULL # removes the entry and shifts subsequent entries down one slot
      num_has_space <- num_has_space - 1
      index_map[candidate[1], candidate[2]] <- NA
      index_map[which(index_map > index)] <- index_map[which(index_map > index)] - 1 # because entries have shifted by one slot
    }
  }
}

################################### 
####### plot the results:
################################### 

par(mfrow = c(2, 2))

# plot the occupied sites,
# after rescaling their values (for clearer plotting):
num_types <- length(unique(as.numeric(sites))) # number of unique fitness values
sites_for_plot <- sites
if(num_types > 2) {
  ll <- sites_for_plot[which(sites_for_plot > 0)]
  sites_for_plot[which(sites_for_plot > 0)] <- min(num_types - 2, 4) * (ll - min(ll)) / (max(ll) - min(ll)) + 1
}
image(sites_for_plot, main = "Max relative\nfitness", col = c("black", heat.colors(num_types - 1)))

# plot the numbers of empty neighbours:
space_grid <- matrix(0, nrow = grid_width, ncol = grid_width)
for(x in 1:grid_width) for(y in 1:grid_width) {
  index <- index_map[x, y]
  if(!is.na(index)) space_grid[x, y] <- how_many_spaces[[index]]
}
image(space_grid, main = "Number of\nempty neighbours", col = c("black", heat.colors(nhood_size)))

# plot the growth curves of population and equivalent radius:
plot(Population ~ Time, data = output_df, type = "l", main = "Population growth")
plot(sqrt(Population / pi) ~ Time, data = output_df, type = "l", ylab = "Equivalent radius", main = "Equivalent\nradius growth")

