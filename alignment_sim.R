# Load required packages
library(animation)
library(tidyverse)
library(qgraph)
library(igraph)



n_updates <- 500

n_dim = 2 #lattice dimensions
length_dim = 10 #number of nodes per dimension. 

lattice = make_lattice(
  length = length_dim,
  dim = n_dim) %>% 
  as.matrix() 
config <-sample(seq(1,20), 
                nrow(lattice), 
                replace=TRUE)
potts_q = 20



config <- simulate_alignment_function(adj_matrix = lattice,
                                      beta = beta,
                                      config = config,
                                      n_updates = n_updates,
                                      alignment_function = "Q_ordered",
                                      Q_states = potts_q)

############################################################################
############################################################################
###                                                                      ###
###                                SIM                                   ###
###                                                                      ###
############################################################################
############################################################################

n_updates <- 500

n_dim = 2 #lattice dimensions
length_dim = 10 #number of nodes per dimension. 

lattice = make_lattice(
  length = length_dim,
  dim = n_dim) %>% 
  as.matrix() 
#qgraph(lattice) #visualizing

ini_config <- sample(c(-1,1), 
                     nrow(lattice), 
                     replace=TRUE)
ini_config <- sample(c(0,1), 
                     nrow(lattice), 
                     replace=TRUE)
ini_config <- sample(c(-1,0,1), 
                     nrow(lattice), 
                     replace=TRUE)

ini_config <- sample(seq(1,4), #potts and ordered
                     nrow(lattice), #potts and ordered
                     replace=TRUE)
ini_config <- sample(seq(1,4), #potts and ordered
                     nrow(lattice), #potts and ordered
                     replace=TRUE)
ini_config <- sample(seq(1,20), 
                     nrow(lattice), 
                     replace=TRUE)
ini_config <- sample(seq(1,20), 
                     nrow(lattice), 
                     replace=TRUE)
ini_config <- runif(nrow(lattice), 
                    0, 
                    2*pi)
#c("Ising", "Potts", "Percolation-Ising", "Blume-capel", "XY", "Q-order"),



setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/test_results")
# Create and save GIF animation
saveGIF({config <- sample(encoding, 
                          nrow(lattice), 
                          replace=TRUE)

  animate(config = config, adj_matrix = lattice)
  for (i in 1:50) { #number of pictures
    
    
    config <- simulate_alignment_function(adj_matrix = lattice,
                       beta = 3,
                       config = config,
                       n_updates = n_updates,
                       alignment_function = "Q-ordered",
                       q_ordered = 20,
                       potts_q = NULL)
    
    animate(config = config, adj_matrix = lattice)
  }
}, interval = 0.2)






##################################################################
##                    Simulation all models                     ##
##################################################################

config_list <- list(sample(c(-1,1), 
                           nrow(lattice), 
                           replace=TRUE),
                    sample(c(0,1), 
                           nrow(lattice), 
                           replace=TRUE),
                    sample(c(-1,0,1), 
                           nrow(lattice), 
                           replace=TRUE),
                    sample(seq(1,4), #potts and ordered
                           nrow(lattice), #potts and ordered
                           replace=TRUE),
                    sample(seq(1,4), #potts and ordered
                           nrow(lattice), #potts and ordered
                           replace=TRUE),
                    sample(seq(1,20), 
                           nrow(lattice), 
                           replace=TRUE),
                    sample(seq(1,20), 
                           nrow(lattice), 
                           replace=TRUE),
                    sample(x = seq(0, 2*pi, length.out = 100), 
                                   size = 100,
                           replace = TRUE)) #continuous is approximated by a 100?? Plausible


models <- c("Ising", "Percolation-ising", "Blume-capel","Q-order", "Potts", "Q-order","Potts", "XY")
Potts_q_list <- c(NA, NA, NA, 4,4,20,20,NA)


res_df <- expand.grid(model = 1:8, beta = c(0.1,2)) %>% 
  mutate(res_mean = NA,
         res_potts = NA,
         potts_q = rep(Potts_q_list,2),
         model_name = rep(models, 2),
         n_distinct = rep(c(2, 2,3,4,4,20,20,100),2))
         #min = rep(c(-1,0,-1,1,1,1,1,0), 2),
         #max = rep(c(1,1,1,4,4,20,20,2*pi),2))
         


ani.options(nmax = 9999) #removes pop up
setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/test_results")
for (m in 8:nrow(res_df)){
    
    saveGIF({
    
      color <- scico(n = res_df$n_distinct[m])
      config <- rep(config_list,2)[[res_df$model[m]]]
      animate(config = config, adj_matrix = lattice, color_use = color)
      
      for (i in 1:50) { #number of pictures
      
      
        config <- simulate_alignment_function(adj_matrix = lattice,
                                              beta = res_df$beta[m],
                                              config = config,
                                              n_updates = n_updates,
                                              alignment_function = res_df$model_name[m],
                                              Q_states = res_df$potts_q[m])
      
      animate(config = config, adj_matrix = lattice, 
              color_use = color)
      
    }
    },
    interval = 0.2,
    movie.name = paste0(res_df$model_name[m],m,"_",res_df$beta[m],".gif"),
    ani.width = 600,
    ani.height = 600
    )
    
    res_df$res_mean[m] <- mean(config)
    if (!is.na( res_df$potts_q[m])) res_df$res_potts[m] <- potts_order_param(n_nodes = nrow(lattice), 
                                                                             Q_states =  res_df$potts_q[m],
                                                                             config = config)
    
}




#################################################################
##                           XY ONLY                           ##
#################################################################


n_updates = 10000

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/test_results")
for (m in c(16)){
  
  saveGIF({
    
    color <- scico(n = res_df$n_distinct[m])
    config <- rep(config_list,2)[[res_df$model[m]]]
    animate(config = config, adj_matrix = lattice, color_use = color)
    
    for (i in 1:100) { #number of pictures
      
      
      config <- simulate_alignment_function(adj_matrix = lattice,
                                            beta = res_df$beta[m],
                                            config = config,
                                            n_updates = n_updates,
                                            alignment_function = res_df$model_name[m],
                                            Q_states = res_df$potts_q[m])
      
      animate(config = config, adj_matrix = lattice, 
              color_use = color)
      
    }
  },
  interval = 0.2,
  movie.name = paste0(res_df$model_name[m],m,"_",res_df$beta[m],".gif"),
  ani.width = 600,
  ani.height = 600
  )
  
  res_df$res_mean[m] <- mean(config)
  if (!is.na( res_df$potts_q[m])) res_df$res_potts[m] <- potts_order_param(n_nodes = nrow(lattice), 
                                                                           Q_states =  res_df$potts_q[m],
                                                                           config = config)
  
}





len = 50
res_df <- data.frame(beta = seq(0,1.5,length.out = len),
                     mean_x = rep(NA,len))
for (i in 1:nrow(res_df)){
  
  potts_q = 3
  
  res<- simulate_alignment_function(
    adj_matrix = adj_matrix,
    beta = res_df$beta[i],
    n_updates = 20000, 
    alignment_function = "Ising",
    potts_q = potts_q)
  
  res_df$mean_x[i]<-  mean(res)
  #res_df$mean_x[i]<-  potts_order_param(n_nodes = length(res),potts_q,config = res)
  
}

res_df %>% 
  ggplot(aes(x = beta, y = mean_x)) +
  geom_point()

# critical point around beta = 0.4


#################################################################
##                            Potts                            ##
#################################################################


len = 50
res_df <- data.frame(beta = seq(-2,2,length.out = len),
                     mean_x = rep(NA,len))
for (i in 1:nrow(res_df)){
  
  potts_q = 3
  
  res<- simulate_alignment_function(
    adj_matrix = adj_matrix,
    beta = res_df$beta[i],
    n_updates = 10000, 
    alignment_function = "Potts",
    potts_q = 4)
  
  res_df$mean_x[i]<-  potts_order_param(n_nodes = length(res),potts_q,config = res)
  
}




##################################################################################################
# Set up parameters
nx <- 20      # Number of columns in the lattice
ny <- 20      # Number of rows in the lattice
Temp <- 0.01      # Temperature parameter

n_levels <- 4  # Number of levels or states each lattice site can take

# Generate initial lattice configuration
data <- matrix(sample(0:(n_levels-1), nx*ny, replace=TRUE), nrow=nx, ncol=ny)

# Define initialization function for the plot
init <- function() {
  ggplot() +
    theme_void() +
    geom_tile(aes(x = rep(1:nx, ny), y = rep(1:ny, each = nx), fill = factor(data)), width = 1, height = 1) +
    scale_fill_manual(values = rainbow(n_levels))
}

# Define animation function for each frame
animate <- function(i) {
  print(i)
  data <- nxt_state(data)
  plot <- ggplot() +
    theme_void() +
    geom_tile(aes(x = rep(1:nx, ny), y = rep(1:ny, each = nx), fill = factor(data)), width = 1, height = 1) +
    scale_fill_manual(values = rainbow(n_levels))
  print(plot)
}

# Define function to calculate nxt state of the lattice
nxt_state <- function(data) {
  for (i in 1:(1e5)) {
    n <- sample(1:nx, 1)                    # Select random column index
    m <- sample(1:ny, 1)                    # Select random row index
    nxt <- sample(0:(n_levels-1), 1)        # Select random new state
    
    # Calculate energy based on neighboring sites
    a <- (data[m, (n-1+nx) %% nx + 1] == nxt) * 2 - 1
    b <- (data[m, (n+1) %% nx + 1] == nxt) * 2 - 1
    c <- (data[(m-1+ny) %% ny + 1, n] == nxt) * 2 - 1
    d <- (data[(m+1) %% ny + 1, n] == nxt) * 2 - 1
    d1 <- (data[(m+1) %% ny + 1, (n+1) %% nx + 1] == nxt) * 2 - 1
    d2 <- (data[(m+1) %% ny + 1, (n-1+nx) %% nx + 1] == nxt) * 2 - 1
    d3 <- (data[(m-1+ny) %% ny + 1, (n+1) %% nx + 1] == nxt) * 2 - 1
    d4 <- (data[(m-1+ny) %% ny + 1, (n-1+nx) %% nx + 1] == nxt) * 2 - 1
    E <- a + b + c + d + d1 + d2 + d3 + d4
    
    # Update the state of the lattice site based on energy comparison
    if (E > 0) {
      data[m, n] <- nxt
    } else {
      r <- runif(1)
      if (E <= 0 && r < exp(E) / Temp) {
        data[m, n] <- sample(0:(n_levels-1), 1)
      }
    }
  }
  
  return(data)
}

# Set animation options
ani.options(interval = 1)

# Create and save GIF animation
saveGIF({
  init()
  for (i in 1:10) { #number of pictures
    animate(i)
  }
}, interval = 0.1)
