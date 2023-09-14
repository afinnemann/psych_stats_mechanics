library(tidyverse)
library(igraph)
library(parSim)

#### Main func

simulate_alignment_function <- function(adj_matrix,
                                        beta,
                                        n_updates, 
                                        initial_config=NULL, 
                                        config,
                                        alignment_function = c("Ising", "Potts", "Percolation-Ising", "Blume-capel", "XY", "Q-order"),
                                        Q_states = NULL) {
  
  ##################################################################
  ##        Glauber update rules for the 5 standard models        ##
  ##################################################################
  
  if (alignment_function == "Ising"){
    encoding = c(-1,1)
    update_element = function(node, 
                              config,
                              adj_matrix, 
                              encoding, 
                              Q_states = NULL)
    {
      
      # delta_E is the difference between the new and current alignment
      delta_E <- 2 * (-config[node] * sum(adj_matrix[node, ] * config)) 
      
      # Calculate the transition probability
      transition_prob <- exp(delta_E * beta) / (1 + exp(delta_E * beta))
      
      # If the new alignment is less than the current alignment, or a random number
      # is less than the transition probability, then flip the node's spin
      #if (delta_E < 0 || runif(1) < transition_prob) {
      if (runif(1) < transition_prob) {
        return(-config[node])  # Flip the spin
      }else{
        config[node]
      }
    } 
    
    #Potts Glauber: The random element's state is updated to one of the Q states
    #with a probability proportional to the number neighbors with that state
  }else if(alignment_function == "Percolation-ising"){
    encoding = c(0,1)
    update_element = function(node, 
                              config,
                              adj_matrix, 
                              encoding, 
                              Q_states = NULL){
      
      # delta_E is the difference between the new and current alignment
      delta_E <- 2 * (-config[node] * sum(adj_matrix[node, ] * config)) 
      
      # Calculate the transition probability
      transition_prob <- exp(delta_E * beta) / (1 + exp(delta_E * beta))
      
      # If the new alignment is less than the current alignment, or a random number
      # is less than the transition probability, then flip the node's spin
      #if (delta_E < 0 || runif(1) < transition_prob) {
      if (runif(1) < transition_prob) {
        return((config[node] - 1)^2)  # Flip the spin
      }else{
        config[node]
      }
    } 
    
  } else if(alignment_function == "Potts"){
    encoding <- seq(1,Q_states)
    update_element <- function(Q_states,
                              node, 
                              config, 
                              adj_matrix, 
                              encoding){
      
      
      # Calculate the current alignment for the node
      current_alignment<- sum(adj_matrix[node, ] * (config == config[node])) 
      
      # Calculate the possible new alignment for each state
      new_alignment <- sapply(encoding, function(new_state) {
        sum(adj_matrix[node, ] * (config == new_state))
      })
      
      # delta_E is the difference between the new and current alignment
      delta_E <- new_alignment - current_alignment
      
      # Calculate the transition probabilities for each state
      transition_probs <- exp(delta_E * beta)
      
      # Normalize the transition probabilities so they sum to 1
      transition_probs <- transition_probs / sum(transition_probs)
      
      # Sample a new state for the node
      new_state <- sample(encoding, 1, prob=transition_probs)
      
      return(new_state)
      
    }
    
  } else if(alignment_function == "Blume-capel"){
    
    encoding <- c(-1,0,1)
    update_element <- function(node, 
                              config,
                              adj_matrix, 
                              encoding, 
                              Q_states = NULL){
      
      current_alignment<- sum(adj_matrix[node, ] * (config == config[node])) 
      
      new_alignment <- sapply(encoding, function(new_state) {
        sum(adj_matrix[node, ] * (config == new_state))
      })
      
      # delta_E is the difference between the new and current alignment
      delta_E <- new_alignment - current_alignment
      
      # Calculate the transition probabilities for each state
      transition_probs <- exp(delta_E * beta)
      
      # Normalize the transition probabilities so they sum to 1
      transition_probs <- transition_probs / sum(transition_probs)
      
      # Sample a new state for the node
      new_state <- sample(encoding, 1, prob=transition_probs)
      
      return(new_state)
      
    }
    
  } else if(alignment_function =="Q-order"){

    encoding <- seq(1,Q_states)
    
    update_element = function(node, 
                              config,
                              adj_matrix, 
                              encoding, 
                              Q_states = NULL){
      
      current_alignment<- sum(adj_matrix[node, ] * (config == config[node])) 
      
      new_alignment <- sapply(encoding, function(new_state) {
        sum(adj_matrix[node, ] * (config == new_state))
      })
      
      # delta_E is the difference between the new and current alignment
      delta_E <- new_alignment - current_alignment
      
      # Calculate the transition probabilities for each state
      transition_probs <- exp(delta_E * beta)
      
      # Normalize the transition probabilities so they sum to 1
      transition_probs <- transition_probs / sum(transition_probs)
      
      # Sample a new state for the node
      new_state <- sample(encoding, 1, prob=transition_probs)
      
      return(new_state)
    }
      
    }else if(alignment_function == "XY"){
      
      encoding = sample(x = seq(0, 2*pi, length.out = 100), size = 1)
      
      update_element <- function(Q_states = NULL,node, config, adj_matrix, encoding){
        
        # Propose a new angle for this node
        #new_angle <- runif(1, min=0, max=2*pi)
        new_angle <- sample(encoding, size = 1)
        
        # Calculate the change in alignment if we update this node's angle
        delta_E <- -sum(adj_matrix[node, ] * (cos(new_angle - config) - cos(config[node] - config)))
        
        # Calculate the transition probability
        transition_prob <- min(1, exp(-delta_E * beta))
        
        # Update the node's angle with this probability
        if (runif(1) < transition_prob) {
          return(new_angle)
        }else{
          return(config[node])
        }
      }
    }else{stop("No proper alignment function chosen")}  
    
  n_nodes <- nrow(adj_matrix)
  
    for (i in 1:n_updates){
      
      node <- sample(1:n_nodes, 1)
      
      # Update the node's state
      config[node] <- update_element(node = node, 
                                     config = config, 
                                     adj_matrix = adj_matrix, 
                                     encoding = encoding, 
                                     Q_states = Q_states)
    }
    return(config)
  }
  
  #library(qgraph)
  #qgraph(adj_matrix, color = config)
  
  
  
  
  ### Misc fun
  
  
  animate <- function(adj_matrix, config, color_use) {
    #qgraph color argument requires all numbers to be positive
    if (min(config) < 1) config = config - min(config) + 2
    
    # Generate colors using scico
    #color_use <- scico(n = length(encoding))
    
    # Map config values to colors
    config_colors <- color_use[match(config, unique(config))]
    
    qgraph(input = adj_matrix, color = config_colors, layout = "spring", theme = "colorblind")
  }
  
  potts_order_param <- function(n_nodes, Q_states, config){
    #from https://physics.stackexchange.com/questions/59663/how-to-define-the-order-parameter-of-the-q-state-potts-model
    #Computes the proportion of nodes in the dominant q state
    state_count <- sapply(1:Q_states, function(x) sum(config == x))
    
    ((Q_states/n_nodes) * max(state_count))/ (Q_states - 1)
  }
  
  
  
  # Define the Blume-Capel alignment function
  blume_capel_alignment <- function(current_state, adjacency_matrix, pars) {
    J <- pars[1]  # Coupling constant
    D <- pars[2]  # Crystal field constant
    
    alignment <- -sum(J * current_state * current_state[adjacency_matrix == 1]) - D * sum(current_state^2)
    return(alignment)
  }
  
  
  H_extend_fun <- function(node, param){
    d <- param[1]
    if (node == 0) d else 0
  }
 
  
  
  