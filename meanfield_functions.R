#Ising Mean field approximation solver
#This function returns a data frame of equlibrium points and stability given beta and average number of neighbors
library(tidyverse)
library(rootSolve)
library(plotly)
library(IsingSampler)
library(parSim)
library(igraph)


#################################################################
##                          Functions                          ##
#################################################################





MFA_sim <- function(root_MFA,
                    alpha_sim = 0,
                    beta_sim = 1,
                    average_density = 4,
                    MFA_to_solve = MFA_to_solve){
  
  
  root_MFA_as.list = list(root_MFA)
  if(MFA_to_solve != "custom"){
    
    result =  parSim(
      ### SIMULATION CONDITIONS
      alpha = alpha_sim,
      beta = beta_sim, 
      average_density = average_density,
      root_MFA_as.list = root_MFA_as.list,
      reps = 1, # repetitions per condition
      write = FALSE, # Writing to a file
      nCores = 1, # Number of cores to use
      expression = {
        
        #print(alpha)
        #print(beta)
        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, average_density = average_density,alpha = alpha, beta = beta)
        
        # Find all roots of the root function in the interval [-1,1]
        all_roots <- uniroot.all(solve_MFA,interval = c(-1, 1))
        
        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "unstable", "stable")),
                          value = all_roots)
      }
    )
    result$value = as.numeric(as.character(result$value))
    return(result)
    
    
  }else if(MFA_to_solve == "custom"){
    result =  parSim(
      ### SIMULATION CONDITIONS
      alpha = alpha_sim,
      beta = beta_sim, 
      average_density = average_density,
      root_MFA_as.list = root_MFA_as.list,
      reps = 1, # repetitions per condition
      write = FALSE, # Writing to a file
      nCores = 1, # Number of cores to use
      expression = {
        
        #print(alpha)
        #print(beta)
        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, param = list(alpha,beta,average_density))
        
        # Find all roots of the root function in the interval [-1,1]
        all_roots <- uniroot.all(solve_MFA,interval = c(-1, 1))
        
        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "unstable", "stable")),
                          value = all_roots)
      }
    )
    result$value = as.numeric(as.character(result$value))
    return(result)
  }
}



plot_MFA_sim <- function(alpha_sim = NA,
                         beta_sim = NA,
                         average_density = NA,
                         MFA_to_solve = c("minus", "zero", "han", 
                                          "zero_corrected", "indepedent_a",
                                          "blume_capel"),
                         view = c("front", "side", "top", "custom"),
                         custom_MFA = NA,
                         custom_pars = NA){
  
  
  #Errors
  
  if (!(MFA_to_solve %in% c("minus", "zero", "han", 
                            "zero_corrected", "indepedent_a",
                            "blume_capel"))) stop("Incorrect mean-field selected")
  
  
  #Define the root function based on the type of mean field to solve for.
  #The root function is the mean field function minus the mean field.
  #The root function is zero when the mean field is equal to the mean field.
  if(MFA_to_solve != "custom")
  { #custom MFA?
    if ((!is.na(alpha_sim[1]) & !is.na(beta_sim[1])) & !is.na(average_density[1])){ #custom parameters?
      
      if (MFA_to_solve == "minus"){
        root_MFA = function (x, average_density, alpha, beta) ((exp(2 * beta * (average_density * x + alpha)) - 1) / (1 + exp(2 * beta * (average_density * x + alpha)))) - x
      } else if (MFA_to_solve == "zero"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + alpha)) / (1 + exp(beta * (average_density * x + alpha)))) - x
      } else if (MFA_to_solve == "han"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + (alpha - 1) / 2)) / (1 + exp(beta * (average_density * x + (alpha - 1) / 2)))) - x
      } else if (MFA_to_solve == "zero_corrected"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + alpha - 0.5 * average_density)) / (1 + exp(beta * (average_density * x + alpha - 0.5 * average_density)))) - x
      } else if (MFA_to_solve == "indepedent_a"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * average_density * x + alpha - 0.5 * average_density) / (1 + exp(beta * average_density * x + alpha - 0.5 * average_density))) - x
      } else if (MFA_to_solve == "blume_capel"){
        #For blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * ((average_density * x) - alpha)) - exp(beta * (( - average_density * x) - alpha))) / (exp(beta * ((average_density * x - alpha))) + exp( beta * ((- average_density * x - alpha))) + 1) - x
      }
      
      #default MFA with custom parameters
      #result <- MFA_sim(alpha_sim = seq(-3,3,length.out = 200),
      #                  beta_sim = seq(0,4,length.out = 200),
      #                  average_density = 2,
      #                  root_MFA = root_MFA)
      result <- MFA_sim(alpha_sim = alpha_sim ,
                        beta_sim = beta_sim,
                        average_density = average_density,
                        root_MFA = root_MFA,
                        MFA_to_solve = MFA_to_solve)
      
    }else{
      #default
      
      if (MFA_to_solve == "Blume_capel"){
        result <- read.csv("Blume_capel_default.csv")
      }
    }
  }else{
    #Custom H with custom parameters  
    
    result <- MFA_sim(alpha_sim = ifelse(is.null(custom_pars[[1]]),NA,custom_pars[[1]]),
                      beta_sim = ifelse(is.null(custom_pars[[2]]),NA,custom_pars[[2]]),
                      average_den = ifelse(is.null(custom_pars[[3]]),NA,custom_pars[[3]]),
                      root_MFA = custom_MFA,
                      MFA_to_solve = MFA_to_solve)
    
  }
  
  
  ###Plotting settings
  
  
  if (view[1] == "front") camera_angle = list(eye = list(x = 2, y = 0, z = 1)) #default
  if (view[1] == "side") camera_angle = list(eye = list(x = 0, y = 2, z = 0))
  if (view[1] == "top") camera_angle = list(eye = list(x = 0, y = 0, z = 2))
  
  
  
  #Plotting
  
  plot_ly(z = ~result$value,
          x = ~result$beta,
          y = ~result$alpha,
          mode="markers",
          type = "scatter3d",
          marker=list(color=~result$value,
                      colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                      showscale=TRUE)
          #line=list(width=2,color='DarkSlateGrey'))
  ) %>% 
    #layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
    layout(scene = list(xaxis = list(title = "beta"),
                        yaxis = list(title = "alpha"),
                        zaxis = list(title = "Mean field"),
                        aspectratio = list(x = 1, y = 1, z = 0.6),
                        camera = camera_angle)) %>% 
    hide_colorbar()
  
}

#rewrite the following in a concise way:
#Universities can play a vital role in effective climate action - from researching and informing solutions and educating students to contribute to a livable future to reducing their own direct environmental impacts and divesting from fossil fuels.



##################################################################
##                           MFA sim                           ##
##################################################################



t <- plot_MFA_sim(MFA_to_solve = "minus",
                  alpha_sim = seq(-3,3,length.out = 30),
                  beta_sim = seq(0,2,length.out = 30),
                  average_density = 2,
                  view = "side")



top <-plot_MFA_sim(MFA_to_solve = "minus",
                   alpha_sim = seq(-3,3,length.out = 300),
                   beta_sim = seq(0,2,length.out = 300),
                   average_density = 2,
                   view = "top")

front <- plot_MFA_sim(MFA_to_solve = "minus",
                      alpha_sim = seq(-3,3,length.out = 300),
                      beta_sim = seq(0,2,length.out = 300),
                      average_density = 2,
                      view = "front")
side <- plot_MFA_sim(MFA_to_solve = "minus",
                     alpha_sim = seq(-3,3,length.out = 300),
                     beta_sim = seq(0,2,length.out = 300),
                     average_density = 2,
                     view = "side")

library(processx)

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/misc")

#install.packages('reticulate')
library(reticulate)

#reticulate::conda_install('r-reticulate', 'python-kaleido')
#reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')
# Call the save_image function
reticulate::py_run_string("import sys")
save_image(p = front, file = "ising_mfa_front.pdf", width = 600, height = 600)
save_image(p = side, file = "ising_mfa_side.pdf", width = 600, height = 600)
save_image(p = top, file = "ising_mfa_top.pdf", width = 600, height = 600)



###### Example custom function

root_MFA = function (mean_field, param){
  
  param1 <- param[1] #alpha
  param2 <- param[2] #beta
  param3 <- param[3] #average density
  
  ((exp(2 * param2 * (param3 * x + param1)) - 1) / (1 + exp(2 * param2 * (param3 * mean_field + param1)))) - mean_field
}


pars = list(seq(-3,3,length.out = 20),
            seq(0,2, length.out = 20),
            3)





#Bernoulli

library(ggplot2)

d <- seq(-5, 5, by = 0.01)
p <- 2*exp(d) / (1 + 2*exp(d))

df <- data.frame(d, p)

ggplot(df, aes(x = d, y = p)) +
  geom_line() +
  labs(x = "d", y = "Percolation probability")







result =  parSim(
  ### SIMULATION CONDITIONS
  alpha_sim = seq(-2,2, by = 0.05), #prior inclusion probability we loop over
  beta_sim = seq(0,2, by = 0.05), #precisions looped over
  reps = 1, # repetitions per condition
  write = FALSE, # Writing to a file
  nCores = 1, # Number of cores to use
  expression = {
    
    out = ising_MFA_sim(MFA_to_solve = "minus",
                        beta = beta_sim, 
                        average_density = 2, 
                        alpha = alpha_sim)
    
  }
)


result$value = as.numeric(as.character(result$value))
#setwd("~/Documents/PhD/growing_networks/Ising_transformation/draft_plots")
#pdf("minus_cusp.pdf",width= 9, height= 5)

plot_ly(z = ~result$value,
        x = ~result$beta,
        y = ~result$alpha,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
  layout(scene = list(xaxis = list(title = "Beta"),
                      yaxis = list(title = "Alpha"),
                      zaxis = list(title = "Mean field"))) %>%
  layout(scene = list(xaxis = list(tickvals = c(0, 1,2)),
                      yaxis = list(tickvals = c(-2, -1, 0, 1,2)),
                      zaxis = list(tickvals = c(-1, 0,1))))

##Rotated cusp
to_zero_value <- function(x) (x + 1) / 2
to_zero_beta <- function(b) 4 * b
to_zero_alpha <- function(a,b) 2  * b * a - 2 * b

#add line to plot_ly of a = - 2beta
plot_ly(z = ~ to_zero_value(result$value),
        x = ~ to_zero_beta(result$beta_sim),
        y = ~ to_zero_alpha(result$alpha_sim, result$beta_sim),
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
  layout(scene = list(xaxis = list(title = "Beta"),
                      yaxis = list(title = "Alpha"),
                      zaxis = list(title = "Mean field"))) %>%
  layout(scene = list(xaxis = list(tickvals = c(2, 4, 8)),
                      yaxis = list(tickvals = c(-5 -3, -1, 0, 1,2)),
                      zaxis = list(tickvals = c( 0,1)))) 


##################################################################
##                          (0,1) case                          ##
##################################################################



result_zero =  parSim(
  ### SIMULATION CONDITIONS
  alpha_sim = seq(-2,2, by = 0.05), #prior inclusion probability we loop over
  beta_sim = seq(0,10, by = 0.05), #precisions looped over
  reps = 1, # repetitions per condition
  write = FALSE, # Writing to a file
  nCores = 1, # Number of cores to use
  expression = {
    
    out = ising_MFA_equlibrium_zero(beta = beta_sim, average_density = 2, alpha = alpha_sim)
    
    #out2 = matrix(c(out[,2] %>% unlist(), 
    #                rep(alpha[i], nrow(out))),
    #              nrow = nrow(out),
    #              ncol = 2)
    
    #res = rbind(res, out2)
    
    
  }
)

plot_ly(z = ~result_zero$value,
        x = ~result_zero$beta_sim,
        y = ~result_zero$alpha_sim,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result_zero$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
  layout(scene = list(xaxis = list(title = "Beta"),
                      yaxis = list(title = "Alpha"),
                      zaxis = list(title = "Mean field"))) %>%
  layout(scene = list(xaxis = list(tickvals = c(2, 4)),
                      yaxis = list(tickvals = c(-2, -1, 0, 1,2)),
                      zaxis = list(tickvals = c(0, 0.5,1))))

#setwd("~/Documents/PhD/growing_networks/Ising_transformation/draft_plots")
#pdf("minus_cusp.pdf",width= 9, height= 5)

#central line







###########################################################################
###########################################################################
###                                                                     ###
###                           ISING MODEL SIM                           ###
###                                                                     ###
###########################################################################
###########################################################################
#' Plot mean spins as a function of beta
#'
#'The function generates a small world network using the igraph package, 
#'#then runs a parallel simulation using the parSim function with the given parameters. 
#'#The simulation uses the Ising model with the IsingSampler function from the IsingSampler package. 
#'#The results of the simulation are plotted using either plotly or ggplot2, 
#'depending on the values of alpha and beta. If alpha is a sequence and beta is a single value, 
#'a 2D plot is generated with alpha on the x-axis and mean spin on the y-axis. 
#'If alpha and beta are both sequences, a 3D plot is generated with alpha, beta, 
#'and mean spin on the x, y, and z axes, respectively. If alpha is a single value and beta is a sequence,
#' a 2D plot is generated with beta on the x-axis and mean spin on the y-axis.
#'
#' @param dim Dimension of small world network
#' @param size Size of small world network
#' @param nei Number of neighbors connected to
#' @param p Probability of rewiring in small world network
#' @param alpha Thresholds for Ising model
#' @param beta_sequence Sequence of betas to loop over
#' @param reps Number of repetitions per beta
#' @param nCores Number of cores to use
#' @param nIter Number of iterations in IsingSampler
#' @param thresholds Thresholds for IsingSampler
#' @param responses Responses for IsingSampler
#' @return ggplot object
#' @export
plot_ising_sim <- function(dim = 1, 
                           size = 40, 
                           nei = 10, 
                           p = 0, 
                           alpha = 0, 
                           beta = 1, 
                           reps = 1, 
                           nCores = 1, 
                           nIter = 100, 
                           responses = c(-1,1),
                           normalise = F
){
  
  # Generate small world network
  sw_network = igraph::sample_smallworld(dim = dim, 
                                         size = size, 
                                         nei = nei, 
                                         p = p) %>%  
    as_adjacency_matrix() %>% 
    as.matrix() #%>% 
  #apply(c(1,2),function(x) x / (size))
  if (normalise == TRUE){
    sw_network <- sw_network %>% apply(c(1,2),function(x) x / size) #normalising factor
  }
  
  sw_network[lower.tri(sw_network)] = t(sw_network)[lower.tri(sw_network)]
  #print(sw_network)
  #Run simulation with given parameters
  
  #assinging sw_network to global environment to use it inside parSim.
  assign('sw_network', sw_network, envir = .GlobalEnv)
  assign('responses', responses, envir = .GlobalEnv)
  
  beta_sim = parSim(
    beta = beta, 
    alpha = alpha,
    #size = size,
    reps = reps, 
    write = FALSE, 
    nCores = nCores, 
    nIter = nIter,
    expression = {
      
      samp = IsingSampler(n = 50, 
                          graph = sw_network,
                          beta = beta, 
                          nIter = nIter, 
                          thresholds = rep(alpha, nrow(sw_network)),
                          responses = responses)
      
      mean_spin = apply(samp, 1, function(x) mean(x)) 
      
      Results <- data.frame(
        "mean_spin" = mean_spin
      )
    }
  )
  print(responses)
  
  #Plots the result. Type depends on wether alpha and beta are constant or vary
  plot <- plot_ly(z = ~ beta_sim$mean_spin,
                  x = ~ beta_sim$beta,
                  y = ~ beta_sim$alpha,
                  mode="markers",
                  type = "scatter3d",
                  marker=list(color=~beta_sim$mean_spin,
                              colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                              showscale=TRUE)
  ) %>% 
    layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))
  
  return(list(plot, beta_sim))
}


size = 100

res200 <- plot_ising_sim(alpha = seq(-1, 1, length.out = 10),
                         beta = seq(0, 3, length.out = 10),
                         responses = c(0, 1L),
                         size = size,
                         nei = size,
                         normalise = TRUE)

size = 30

res30 <- plot_ising_sim(alpha = seq(-1, 1, length.out = 20),
                        beta = seq(0, 2, length.out = 30),
                        responses = c(-1L, 1L),
                        size = size,
                        nei = size,
                        normalise = TRUE)











#Extend the plot_ly plot with the points specified by bifurcation_set_y and x in grid
grid <- expand.grid(x = seq(-3,3, length.out = 100), #create grid of numbers
                    y = seq(0, 4, length.out = 40)) %>% 
  
  #computing bifurcation set: isolating beta in 4bta^3 = 27alpha^2
  mutate(bifurcation_set_y = ((27 * (x^2) )/4)^(1/3),
         #transform using: beta^new = 4 * beta_old
         transformed_y <- 4 * bifurcation_set_y,
         transformed_x = 2 * x - 2 * bifurcation_set_y,
  )

#how to change this plot so it only shows alpha_sim > 0



alpha_pos <- result_zero %>% 
  filter(alpha_sim > 0)
#pdf("alpha_pos.pdf", width = 12, height = 8)
plot_ly(z = ~alpha_pos$value,
        x = ~alpha_pos$beta_sim,
        y = ~alpha_pos$alpha_sim,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~alpha_pos$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) 
#dev.off()


result_zero %>% 
  mutate(value = value %>% as.character() %>% as.numeric()) %>% 
  filter(value > 0.1 & (value < 0.9)) %>% 
  filter(beta_sim > 2) %>% 
  group_by(beta_sim) %>% 
  summarise(min_alpha = min(alpha_sim),
            max_alpha = max(alpha_sim)) -> t2

#view of bifurcation set
t2 %>% 
  ggplot() +
  geom_line(aes(x = beta_sim, y = max_alpha)) +
  geom_line(aes(x = beta_sim, y = min_alpha))

#The plot_ly plot is extended with the points specified by bifurcation_set_y and x in grid. The points are added as a trace with mode set to "lines" and line width and color set to 2 and "DarkSlateGrey" respectively.

#solve for a (27a^2 /4)^(1/3)

#solve for a (27a^2 /4)^(1/3)


#a = (4/27)^(1/2) * (27a^2 /4)^(1/3)



result_zero = result_zero %>% 
  mutate(value = value %>% as.character() %>% as.numeric(),
         beta_sim = beta_sim %>% as.character() %>% as.numeric(),
         alpha_sim = alpha_sim %>% as.character() %>% as.numeric(),
         value_trans = (value *2)-1,
         beta_trans = beta_sim / 4,
         alpha_trans = (alpha_sim + 2 * beta_sim) / 2)


plot_ly(z = ~result_zero$value_trans,
        x = ~result_zero$beta_trans,
        y = ~result_zero$alpha_trans,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result_zero$value_trans,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))




###########################################################################
###########################################################################
###                                                                     ###
###                           ISING MODEL SIM                           ###
###                                                                     ###
###########################################################################
###########################################################################
#' Plot mean spins as a function of beta
#'
#'The function generates a small world network using the igraph package, 
#'#then runs a parallel simulation using the parSim function with the given parameters. 
#'#The simulation uses the Ising model with the IsingSampler function from the IsingSampler package. 
#'#The results of the simulation are plotted using either plotly or ggplot2, 
#'depending on the values of alpha and beta. If alpha is a sequence and beta is a single value, 
#'a 2D plot is generated with alpha on the x-axis and mean spin on the y-axis. 
#'If alpha and beta are both sequences, a 3D plot is generated with alpha, beta, 
#'and mean spin on the x, y, and z axes, respectively. If alpha is a single value and beta is a sequence,
#' a 2D plot is generated with beta on the x-axis and mean spin on the y-axis.
#'
#' @param dim Dimension of small world network
#' @param size Size of small world network
#' @param nei Number of neighbors connected to
#' @param p Probability of rewiring in small world network
#' @param alpha Thresholds for Ising model
#' @param beta_sequence Sequence of betas to loop over
#' @param reps Number of repetitions per beta
#' @param nCores Number of cores to use
#' @param nIter Number of iterations in IsingSampler
#' @param thresholds Thresholds for IsingSampler
#' @param responses Responses for IsingSampler
#' @return ggplot object
#' @export
plot_ising_sim <- function(dim = 1, 
                           size = 40, 
                           nei = 10, 
                           p = 0, 
                           alpha = 0, 
                           beta = 1, 
                           reps = 1, 
                           nCores = 1, 
                           nIter = 100, 
                           responses = c(-1,1),
                           normalise = F
){
  
  # Generate small world network
  sw_network = igraph::sample_smallworld(dim = dim, 
                                         size = size, 
                                         nei = nei, 
                                         p = p) %>%  
    as_adjacency_matrix() %>% 
    as.matrix() #%>% 
  #apply(c(1,2),function(x) x / (size))
  if (normalise == TRUE){
    sw_network <- sw_network %>% apply(c(1,2),function(x) x / size) #normalising factor
  }
  
  sw_network[lower.tri(sw_network)] = t(sw_network)[lower.tri(sw_network)]
  #print(sw_network)
  #Run simulation with given parameters
  
  #assinging sw_network to global environment to use it inside parSim.
  assign('sw_network', sw_network, envir = .GlobalEnv)
  assign('responses', responses, envir = .GlobalEnv)
  
  beta_sim = parSim(
    beta = beta, 
    alpha = alpha,
    #size = size,
    reps = reps, 
    write = FALSE, 
    nCores = nCores, 
    nIter = nIter,
    expression = {
      
      samp = IsingSampler(n = 50, 
                          graph = sw_network,
                          beta = beta, 
                          nIter = nIter, 
                          thresholds = rep(alpha, nrow(sw_network)),
                          responses = responses)
      
      mean_spin = apply(samp, 1, function(x) mean(x)) 
      
      Results <- data.frame(
        "mean_spin" = mean_spin
      )
    }
  )
  print(responses)
  
  #Plots the result. Type depends on wether alpha and beta are constant or vary
  plot <- plot_ly(z = ~ beta_sim$mean_spin,
                  x = ~ beta_sim$beta,
                  y = ~ beta_sim$alpha,
                  mode="markers",
                  type = "scatter3d",
                  marker=list(color=~beta_sim$mean_spin,
                              colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                              showscale=TRUE)
  ) %>% 
    layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))
  
  return(list(plot, beta_sim))
}


size = 200

res200 <- plot_ising_sim(alpha = seq(-0.3, 0.3, length.out = 20),
                         beta = seq(0.5, 1.5, length.out = 30),
                         responses = c(-1L, 1L),
                         size = size,
                         nei = size,
                         normalise = TRUE)

size = 30

res30 <- plot_ising_sim(alpha = seq(-0.3, 0.3, length.out = 20),
                        beta = seq(0, 2, length.out = 30),
                        responses = c(-1L, 1L),
                        size = size,
                        nei = size,
                        normalise = TRUE)

res30_2 <- plot_ising_sim(alpha = seq(-0.3, 0.3, length.out = 20),
                          beta = seq(0, 2, length.out = 30),
                          responses = c(-1L, 1L),
                          size = size,
                          nei = size,
                          normalise = TRUE)


df_res <- cbind(res[[2]], res_zero[[2]])



plot <- plot_ly(z = ~beta_sim$mean_spin,
                x = ~beta_sim$beta,
                y = ~beta_sim$alpha,
                mode="markers",
                type = "scatter3d",
                marker=list(color=~beta_sim$mean_spin,
                            colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                            showscale=TRUE)
                #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))




###########################################################################
###########################################################################
###                                                                     ###
###                        VECTOR FIELD ANALYSIS                        ###
###                                                                     ###
###########################################################################
###########################################################################

#write the following code in a simpler way in R
grid <- expand.grid(x = seq(-3,3, length.out = 1000)) %>%  #create grid of numbers
  #computing bifurcation set: isolating beta in 4bta^3 = 27alpha^2
  mutate(bifurcation_set_y = ((27 * (x^2) )/4)^(1/3),
         #transform using: beta^new = 4 * beta_old
         transformed_y = 4 * bifurcation_set_y,
         transformed_x = 2 * x - 2 * bifurcation_set_y,
  )

grid %>% ggplot(aes(x = x, y = bifurcation_set_y)) + #plot results
  geom_line(color = "black") +
  geom_point(aes(x = transformed_x, y = transformed_y),color = "red", size = 2)  +
  xlim(-3,3) + ylim(0,4) + xlab(expression(~ alpha)) + ylab(expression(~ beta)) + theme_classic()



a_fun <- function(b = seq(0, 14, by = 0.1)) {
  a_plus <- b * sqrt(b) /(6 * sqrt(3)) - b/2
  a_minus <- b * (- sqrt(b)) /(6 * sqrt(3)) - b/2
  return(data.frame("beta" = b, "alpha_plus" = a_plus, "alpha_minus" = a_minus))
}
res = a_fun() #run function

#Transform (0,Y) using x_new = 2x - 2y and y_new = 4 y ->
# transformed (-2y, 4Y)

pdf("bifurcation_sets.pdf", width = 12, height = 8)
#incrase axis label size on the followin plot
grid %>% ggplot(aes(x = x, y = bifurcation_set_y)) + #plot results
  geom_line(size = 2, color = "black") +
  geom_point(aes(x = transformed_x, y = transformed_y),color = "orange", size = 2)  +
  #geom_point(data = res, aes(x = alpha_plus, y = beta), color = "orange") +
  #geom_point(data = res, aes(x = alpha_minus, y = beta), color = "orange") +
  xlim(-3,3) + ylim(0,4) + xlab(expression(~ alpha)) + ylab(expression(~ beta)) + theme_classic()  +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
#geom_abline(intercept = 0, slope = -2, color = "blue")
dev.off()






theme(axis.title = element_text(size = 20))
grid %>% ggplot(aes(x = transformed_x, y = transformed_y)) +
  geom_point()


#Compute b for each alpha value
#Write a function that computes b for the range of alpha -10 to 10 in R
#b = (108 a^2 + (9072 a^3)/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + (122472 a^2)/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + (472392 a)/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + 81 (128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + 531441/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + 2916 a + 6561)^(1/3)

#change the color of the line to blue and make it smaller. Also detect and fix bugs
plot_ly(z = ~result_zero$value,
        y = ~result_zero$beta_sim,
        x = ~result_zero$alpha_sim,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result_zero$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  #change the trace to a line. Also detect and fix bugs in the following code
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
  add_trace(x = ~res$alpha_plus,
            y = ~res$beta,
            z = ~rep(0, nrow(res)),
            mode = "lines",
            line = list(width = 0.2, color = "blue")) %>% 
  add_trace(x = ~res$alpha_minus,
            y = ~res$beta,
            z = ~rep(1, nrow(res)),
            mode = "lines",
            line = list(width = 0.2, color = "blue")) 

#Write a function that computes b for the range of alpha -10 to 10 in R
#b = (108 a^2 + (9072 a^3)/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + (122472 a^2)/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + (472392 a)/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + 81 (128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + 531441/(128 a^5 + 6264 a^4 + 71928 a^3 + 341172 a^2 + 8 sqrt(256 a^10 + 3104 a^9 + 11745 a^8 + 13122 a^7) + 708588 a + 531441)^(1/3) + 2916 a + 6561)^(1/3)


b <- function(a = seq(-10, 10, by = 0.1)) {
  computed_beta <- (108*a^2 + (9072*a^3)/(128*a^5 + 6264*a^4 + 71928*a^3 + 341172*a^2 + 8*sqrt(256*a^10 + 3104*a^9 + 11745*a^8 + 13122*a^7) + 708588*a + 531441)^(1/3) + 
                      (122472*a^2)/(128*a^5 + 6264*a^4 + 71928*a^3 + 341172*a^2 +8*sqrt(256*a^10 + 3104*a^9 + 11745*a^8 + 13122*a^7) + 708588*a + 531441)^(1/3) + 
                      (472392*a)/(128*a^5 + 6264*a^4 + 71928*a^3 + 341172*a^2 + 8*sqrt(256*a^10 + 3104*a^9 + 11745*a^8 + 13122*a^7) + 708588*a + 531441)^(1/3) + 
                      81*(128*a^5 + 6264*a^4 + 71928*a^3 + 341172*a^2 + 8*sqrt(256*a^10 + 3104*a^9 + 11745*a^8 + 13122*a^7) + 708588*a + 531441)^(1/3) + 531441/(128*a^5 + 6264*a^4 + 71928*a^3 + 341172*a^2 + 8*sqrt(256*a^10 + 3104*a^9 + 11745*a^8 + 13122*a^7) + 708588*a + 531441)^(1/3) + 
                      2916*a + 6561)^(1/3)
  
  return(data.frame("beta" = computed_beta, "alpha" = a))
}

res <- b()

res %>% 
  ggplot(aes(x = alpha, y = beta)) +
  geom_point()



###
size = 20
sw_network = igraph::sample_smallworld(dim = 1, 
                                       size = size, 
                                       nei = 20, 
                                       p = 0.4) %>%  
  as_adjacency_matrix() %>% 
  as.matrix() #%>% 

sw_network[lower.tri(sw_network)] = t(sw_network)[lower.tri(sw_network)]
#Run simulation with given parameters

#The size of network is sum(1:n). The contribution to the hamiliation for a fully connected positve node
sum(sw_network) / 2 == (size - 1) * ((size - 1) + 1) / 2

##validing predictions of the model

betas = seq(0.1, 1, by = 0.001)
res = matrix(NA,nrow = 0, ncol = 2)
average_density = 2


for (i in 1:length(betas)) {
  out = ising_MFA_equlibrium_minus(beta = betas[i], average_density = 2)
  
  out2 = matrix(c(out[,2] %>% unlist(), 
                  rep(betas[i], nrow(out))),
                nrow = nrow(out),
                ncol = 2)
  
  res = rbind(res, out2)
  colnames(res) = c("value", "beta")
}

#setwd("~/Documents/PhD/growing_networks/Ising_transformation/draft_plots")
#pdf("minus_beta.pdf",width= 9, height= 5)
res %>%
  as.data.frame() %>%
  mutate(value = as.numeric(value),
         beta = as.numeric(beta)) %>%
  ggplot(aes(beta, value)) +
  geom_point() +
  theme_classic()
#dev.off()

#for 2-D case: average density of 4
#We find critical temperature (1/beta) is 4. Accords with http://tuvalu.santafe.edu/~simon/practical.pdf




##Alpha simulation
beta = 1
alpha = seq(-4, 4, by = 0.001)
res = matrix(NA,nrow = 0, ncol = 2)

for (i in 1:length(alpha)) {
  out = ising_MFA_equlibrium_minus(beta = beta, average_density = average_density, alpha = alpha[i])
  
  out2 = matrix(c(out[,2] %>% unlist(), 
                  rep(alpha[i], nrow(out))),
                nrow = nrow(out),
                ncol = 2)
  
  res = rbind(res, out2)
  colnames(res) = c("value", "alpha")
}

#setwd("~/Documents/PhD/growing_networks/Ising_transformation/draft_plots")
#pdf("minus_alpha.pdf",width= 9, height= 5)
res %>%
  as.data.frame() %>%
  mutate(value = as.numeric(value),
         alpha = as.numeric(alpha)) %>%
  ggplot(aes(alpha, value)) +
  geom_point() +
  theme_classic()
#dev.off()
#for 2-D case: average density of 4
#We find critical temperature (1/beta) is 4. Accords with http://tuvalu.santafe.edu/~simon/practical.pdf



### Dynamics as a function of beta

betas = seq(0.1, 10, by = 0.01)

res = matrix(NA,nrow = 0, ncol = 2)
average_density = 2
alpha_use = -2

for (i in 1:length(betas)) {
  out = ising_MFA_equlibrium_zero(beta = betas[i], average_density = average_density, alpha = alpha_use)
  
  out2 = matrix(c(out[,2] %>% unlist(), 
                  rep(betas[i], nrow(out))),
                nrow = nrow(out),
                ncol = 2)
  
  res = rbind(res, out2)
  colnames(res) = c("value", "beta")
}

res %>%
  as.data.frame() %>%
  mutate(value = as.numeric(value),
         beta = as.numeric(beta)) %>%
  ggplot(aes(beta, value)) +
  geom_point() +
  theme_classic()


##Alpha simulation
beta = 3
alpha = seq(-5, 5, by = 0.01)
res = matrix(NA,nrow = 0, ncol = 2)

for (i in 1:length(alpha)) {
  out = ising_MFA_equlibrium_zero(beta = beta, average_density = average_density, alpha = alpha[i])
  
  out2 = matrix(c(out[,2] %>% unlist(), 
                  rep(alpha[i], nrow(out))),
                nrow = nrow(out),
                ncol = 2)
  
  res = rbind(res, out2)
  colnames(res) = c("value", "alpha")
}

res %>%
  as.data.frame() %>%
  mutate(value = as.numeric(value),
         alpha = as.numeric(alpha)) %>%
  ggplot(aes(alpha, value)) +
  geom_point()
#for 2-D case: average density of 4
#We find critical temperature (1/beta) is 4. Accords with http://tuvalu.santafe.edu/~simon/practical.pdf







##################################################################
##                  Simulation of Cramer model                  ##
##################################################################


ising_MFA_equlibrium_zero_cramer <- function(threshold, alpha, beta, average_density) {
  #Ising mean field approximation
  #rotated by - x to find roots
  
  #cramer_fun <- function(threshold, alpha, beta, average_density,x) threshold + alpha - (average_density * x * beta)
  
  root_MFA <- function (x) ( 1 /(1 +exp(alpha -  x * beta)))  - x
  #root_MFA <- function (x) ( 1 /(1 +exp(cramer_fun(threshold, alpha, beta, average_density,x))) ) - x
  
  #MFA solved for all roots relevant areas
  all_roots <- uniroot.all(root_MFA, c(-2, 2))
  
  #determining stable equilibriums
  #testing if function is negative prior to fix point.
  #if negative value means curve is going upwards -> unstable
  #if positive value then curve is downwards -> unstable
  lapply(all_roots, function(x) root_MFA(x - 0.001)) %>% 
    lapply( function(x) ifelse(x > 0, "stable", "unstable")) %>% 
    unlist() %>%
    cbind(all_roots) %>% 
    as.data.frame()-> res
  
  colnames(res) = c("stability", "value")
  
  return(res)
}

result_cramer_3 =  parSim(
  ### SIMULATION CONDITIONS
  threshold = -0,
  alpha_sim = seq(-5,20, by = 0.1), #prior inclusion probability we loop over
  beta_sim = seq(0,20, by = 0.05), #precisions looped over
  reps = 1, # repetitions per condition
  write = FALSE, # Writing to a file
  nCores = 1, # Number of cores to use
  expression = {
    
    out = ising_MFA_equlibrium_zero_cramer(beta = beta_sim, 
                                           average_density = 1, 
                                           alpha = alpha_sim,
                                           threshold = threshold)
    
  }
)

plot_ly(z = ~result_cramer_3$value,
        x = ~result_cramer_3$beta,
        y = ~result_cramer_3$alpha,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result_cramer_3$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))
#layout(scene = list(xaxis = list(title = "Beta"),
#                    yaxis = list(title = "Alpha"),
#                    zaxis = list(title = "Mean field"))) %>%
#layout(scene = list(xaxis = list(tickvals = c(0, 1,2)),
#                    yaxis = list(tickvals = c(-2, -1, 0, 1,2)),
#                    zaxis = list(tickvals = c(-1, 0,1))))




######## Cramer derivative solver


ising_MFA_equlibrium_zero_cramer_deriv <- function(threshold, alpha, beta, average_density) {
  #Ising mean field approximation
  #rotated by - x to find roots
  
  #print(c(alpha, beta))
  
  #cramer_fun <- function(threshold, alpha, beta, average_density,x)  alpha -  x * beta
  
  root_MFA <- function (x) (( beta * exp(alpha - beta * x)) / (1 + exp(alpha - beta * x))^2) - 1
  #1 /(1 +exp(cramer_fun(threshold, alpha, beta, average_density,x))) ) - x
  
  #MFA solved for all roots relevant areas
  all_roots <- uniroot.all(root_MFA, c(0, 1))
  
  #print(length(all_roots) )
  #determining stable equilibriums
  #testing if function is negative prior to fix point.
  #if negative value means curve is going upwards -> unstable
  #if positive value then curve is downwards -> unstable
  if (length(all_roots) > 0){
    lapply(all_roots, function(x) root_MFA(x - 0.001)) %>% 
      lapply( function(x) ifelse(x > 0, "stable", "unstable")) %>% 
      unlist() %>%
      cbind(all_roots) %>% 
      as.data.frame()-> res
    
    colnames(res) = c("stability", "value")
    
    return(res)
  }else{return(data_frame("stability" = NA, "value" = NA))}
}


result_cramer_deriv =  parSim(
  ### SIMULATION CONDITIONS
  threshold = -0,
  alpha_sim = seq(-5,20, by = 0.5), 
  beta_sim = seq(0,20, by = 0.5), 
  reps = 1, # repetitions per condition
  write = FALSE, # Writing to a file
  nCores = 1, # Number of cores to use
  expression = {
    
    out = ising_MFA_equlibrium_zero_cramer_deriv(beta = beta_sim, 
                                                 average_density = 1, 
                                                 alpha = alpha_sim,
                                                 threshold = threshold)
    
  }
)

result_cramer_3$value = as.numeric(as.character(result_cramer_3$value))

result_cramer_deriv_clean = result_cramer_deriv %>% 
  filter(complete.cases(result_cramer_deriv))

range(result_cramer_deriv_clean$beta_sim)
range(result_cramer_deriv_clean$alpha_sim)

plot_ly(z = ~result_cramer_3$value,
        x = ~result_cramer_3$beta,
        y = ~result_cramer_3$alpha,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result_cramer_3$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  add_trace(z = ~ result_cramer_deriv_clean$value,#rep(0, length(result_cramer_deriv_clean$value)),
            x = ~result_cramer_deriv_clean$beta_sim,
            y = ~result_cramer_deriv_clean$alpha_sim,
            mode="markers",
            type = "scatter3d",
            marker=list(color=~ "black",#result_cramer_deriv_clean$value,#"black",
                        colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                        showscale=TRUE)
            #line=list(width=2,color='DarkSlateGrey'))
  ) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))

#Count roots

#Extend the R code so it only keeps rows under one condition. 
#The condition is that there are multiple rows with the "id" column
cramer_bif_set <- result_cramer_deriv_clean %>% 
  group_by(id) %>% 
  filter(row_number() == 2) #keep the second row of rows with multiple od


plot_ly(z = ~result_cramer_3$value,
        x = ~result_cramer_3$beta_sim,
        y = ~result_cramer_3$alpha_sim,
        mode="markers",
        type = "scatter3d",
        marker=list(color=~result_cramer_3$value,
                    colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                    showscale=TRUE)
        #line=list(width=2,color='DarkSlateGrey'))
) %>% 
  add_trace(z = ~ rep(0,length(cramer_bif_set$value)),#rep(0, length(result_cramer_deriv_clean$value)),
            x = ~cramer_bif_set$beta_sim,
            y = ~cramer_bif_set$alpha_sim,
            mode="markers",
            type = "scatter3d",
            marker=list(color=~ "black",#result_cramer_deriv_clean$value,#"black",
                        colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                        showscale=TRUE)
            #line=list(width=2,color='DarkSlateGrey'))
  ) %>% 
  layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6)))


