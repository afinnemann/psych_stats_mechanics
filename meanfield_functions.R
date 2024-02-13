#Mean field approximation solver
#This function returns a data frame of equlibrium points and stability given beta and average number of neighbors
library(tidyverse)
library(rootSolve)
library(plotly)
library(IsingSampler)
library(parSim)
library(igraph)
library(scales)



#################################################################
##                          Functions                          ##
#################################################################

### Note Continuous 01 Ising model can give erros. 
#I.e. alpha = (-3,3, length.out =100) does not work - root solver doesnt find solution. but alpha = (-3,3, length.out =102) works


MFA_sim <- function(alpha_sim = NA,
                    beta_sim = NA,
                    average_density = NA,
                    MFA_to_solve = c("Ising", "Percolation", "han", 
                                     "zero_corrected", "independent_a",
                                     "Blume_capel", "Potts",  "Continuous_01"),
                    view = c("front", "side", "top", "custom"),
                    custom_MFA = NA,
                    custom_pars = NA){
  
  
  #Errors
  if (!(MFA_to_solve %in% c("Ising", "Percolation", "han", 
                            "zero_corrected", "independent_a",
                            "Blume_capel", "Potts",  "Continuous_01"))) stop("Incorrect mean-field selected")
  
  
  #Define the root function based on the type of mean field to solve for.
  #The root function is the mean field function minus the mean field.
  #The root function is zero when the mean field is equal to the mean field.
  if(MFA_to_solve != "custom"){ #custom MFA?
    if ((!is.na(alpha_sim[1]) & !is.na(beta_sim[1])) & !is.na(average_density[1])){ #custom parameters? 
      
      
      if (MFA_to_solve == "Ising"){
        root_MFA = function (x, average_density, alpha, beta) ((exp(2 * beta * (average_density * x + alpha)) - 1) / (1 + exp(2 * beta * (average_density * x + alpha)))) - x
      } else if (MFA_to_solve == "Percolation"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + alpha)) / (1 + exp(beta * (average_density * x + alpha)))) - x
      } else if (MFA_to_solve == "han"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + (alpha - 1) / 2)) / (1 + exp(beta * (average_density * x + (alpha - 1) / 2)))) - x
      } else if (MFA_to_solve == "zero_corrected"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * (average_density * x + alpha - 0.5 * average_density)) / (1 + exp(beta * (average_density * x + alpha - 0.5 * average_density)))) - x
      } else if (MFA_to_solve == "independent_a"){
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * average_density * x - alpha - 0.5 * average_density) / (1 + exp(beta * average_density * x - alpha - 0.5 * average_density))) - x
      } else if (MFA_to_solve == "Blume_capel"){
        #For Blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta) (exp(beta * ((average_density * x) - alpha)) - exp(beta * (( - average_density * x) - alpha))) / (exp(beta * ((average_density * x - alpha))) + exp( beta * ((- average_density * x - alpha))) + 1) - x
      } else if (MFA_to_solve == "Potts"){
        #For Blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta) (1 /(1 + (alpha - 1.0001) * exp((beta * (1 - alpha * x)/(alpha - 1.0001))))+0.0001) - x
      } else if (MFA_to_solve == "Continuous_01"){
        #For Blume capel, beta is the coupling constant, and alpha is the costs to presence. Temp is assumed to be 1
        root_MFA = function (x, average_density, alpha, beta){
          # Define theta
          theta = beta * (alpha + average_density * x) 
          
          # continuous Ising in form: 0 = equation - x
          result = ((1 / theta) * (1 - (theta * exp(-theta)) / (1 - exp(-theta)))) - x
          
          return(result)
        }
      }
      
      #default MFA with custom parameters
      #result <- MFA_sim(alpha_sim = seq(-3,3,length.out = 200),
      #                  beta_sim = seq(0,4,length.out = 200),
      #                  average_density = 2,
      #                  root_MFA = root_MFA)
      result <- MFA_solver(alpha_sim = alpha_sim ,
                           beta_sim = beta_sim,
                           average_density = average_density,
                           root_MFA = root_MFA,
                           MFA_to_solve = MFA_to_solve)
      
    }else{
      #load one of the default data sets to avoid simulating
      
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
}


MFA_solver <- function(root_MFA,
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
        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, average_density, alpha, beta)
        #print(beta * (alpha + average_density * x)+ 0.0001)
        # Find all roots of the root function in the interval [-1,1]
        all_roots <- uniroot.all(solve_MFA,interval = c(-1, 1))
        
        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "stable", "unstable")),
                          value = all_roots)
        #res <- data.frame(value = all_roots)
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
      #root_MFA = root_MFA(),
      reps = 1, # repetitions per condition
      write = FALSE, # Writing to a file
      nCores = 1, # Number of cores to use
      expression = {
        
        #solve_MFA <- function(x)  root_MFA_as.list[[1]](x, param = list(alpha,beta,average_density))
        solve_MFA <- function(x)  root_MFA_as.list[[1]](x, average_density = average_density,alpha = alpha, beta = sim_df)
        
        # Find all roots of the root function in the interval [-1,1]
        all_roots <- uniroot.all(solve_MFA,interval = c(-10, 10))
        
        # Determine the stability of each root by testing if the function is positive or negative prior to the root
        res <- data.frame(stability = sapply(all_roots, function(x) ifelse(solve_MFA(x - 0.001) > 0, "stable", "unstable")),
                          value = all_roots)
      }
    )
    result$value = as.numeric(as.character(result$value))
    return(result)
  }
}


plot_MFA <- function(result,
                     camera_x = 2,
                     camera_y = 0,
                     camera_z = 0,
                     y_label = "Alpha",
                     x_label = "Beta",
                     z_label = "Mean field"){
  
  
  ###Plotting settings
  
  
  #if (view[1] == "front") camera_angle = list(eye = list(x = 2, y = 0, z = 0)) #default
  #if (view[1] == "side") camera_angle = list(eye = list(x = 0, y = 2, z = 0))
  #if (view[1] == "top") camera_angle = list(eye = list(x = 0, y = 0, z = 2))
  #if (view[1] == "other-side") camera_angle = list(eye = list(x = 0, y = -2, z = 0))
  
  
  #Coloring based on stability and magnitude
  reds <- colorRampPalette(c("lightpink", "darkred"))(100) #color pallete for unstable
  purp <- colorRampPalette(c("lavender", "darkorchid4"))(100) #color palette for stable
  
  # Creating a color mapping function
  color_reds <- col_numeric(palette = reds, domain = c(min(result$value), max(result$value)))
  color_purp <- col_numeric(palette = purp, domain = c(min(result$value), max(result$value)))
  
  # Apply the color function based on stability
  result$color <- ifelse(result$stability == "stable", color_purp(result$value), color_reds(result$value))
  
  plot_ly(z = ~result$value,
          x = ~result$beta,
          y = ~result$alpha,
          mode="markers",
          type = "scatter3d",
          marker=list(color=result$color, #~result$stability,
                      size = 6,
                      #colorscale=c("rgb(244, 244, 244)","rgb(65, 65, 65)"),
                      showscale=TRUE)
          #line=list(width=2,color='DarkSlateGrey'))
  ) %>% 
    #layout(scene = list(aspectratio = list(x = 1, y = 1, z = 0.6))) %>%
    layout(scene = list(xaxis = list(title = x_label),
                        yaxis = list(title = y_label),
                        zaxis = list(title = z_label),
                        aspectratio = list(x = 1, y = 1, z = 0.6),
                        camera = list(eye = list(x = camera_x, y = camera_y, z = camera_z)))) %>% 
    hide_colorbar()
  
}

