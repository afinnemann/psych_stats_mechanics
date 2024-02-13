library(processx)
library(reticulate)
library(tidyverse)

set.seed(183)

##################################################################
##                         Illustration                        ##
##################################################################


#For potts alpha has to be integers

MFA_output <- MFA_sim(MFA_to_solve = "Blume_capel",
                     alpha_sim = seq(-3,3,length.out = 500),#seq(-3,3, length.out = 50), 
                     beta_sim = seq(0,4,length.out =500),
                     average_density = 2)
#Runs fast for length.ouput up till around 200-300. A couple of minutes 300 to 600. 
# 500 * 500 works, but does not work for 600 * 600. 


#
#write.csv(MFA_output %>% select(- c(root_MFA_as.list,errorMessage)),
#          "Blume_capel_default_500x500.csv")

bc <- read.csv("Blume_capel_default_500x500.csv")

plot_MFA(result = MFA_output,
         camera_x = 2,
         camera_y = 0,
         camera_z = 1)

##################################################################
##                            Sims                              ##
##################################################################
setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/mfa_sim_default_df")

Ising_MFA <- MFA_sim(MFA_to_solve = "Ising",
                     alpha_sim = seq(-3,3,length.out = 500),
                     beta_sim = seq(0,1,length.out =500),
                     average_density = 2)
write.csv(Ising_MFA %>% select(- c(root_MFA_as.list,errorMessage)),
          "Ising_default_500x500.csv")

Percolation_MFA <- MFA_sim(MFA_to_solve = "Percolation",
                     alpha_sim = c(seq(-3,-1.5,length.out = 100), seq(-1.5,-0.5, length.out = 400), seq(-0.5,3,length.out = 200)),
                     beta_sim = seq(0,4,length.out = 300),
                     average_density = 2)
write.csv(Percolation_MFA %>% select(- c(root_MFA_as.list,errorMessage)),
          "Percolation_default_500x500.csv")

BC_output <- MFA_sim(MFA_to_solve = "Blume_capel",
                      alpha_sim = c(seq(-3,0.8,length.out = 200), seq(0.8,1, length.out = 200), seq(1,3,length.out = 100)),#seq(-3,3, length.out = 50), 
                      beta_sim = c(seq(0,0.5,length.out = 50), seq(0.5,0.8, length.out = 200), seq(0.8,4,length.out = 200)),
                      average_density = 2)
write.csv(BC_output %>% select(- c(root_MFA_as.list,errorMessage)),
          "Blume_capel_default_500x500.csv")


continuous01_output <- MFA_sim(MFA_to_solve = "Percolation",
                     alpha_sim = seq(-2, 2, length.out = 200), 
                     beta_sim = seq(0.0001,4, length.out = 200),
                     average_density = 2)

plot_MFA(result = continuous01_output,
         camera_x = 1.8,
         camera_y = 0,
         camera_z = 0.2)


write.csv(continuous01_output %>% select(- c(root_MFA_as.list,errorMessage)),
          "continuous01_default_300x300.csv")


t <- continuous01_output %>% filter(round(alpha,1) == 0)


##################################################################
##                        Visualize                             ##
##################################################################


Ising_MFA <- read.csv("Ising_default_500x500.csv")
front <- plot_MFA(result = Ising_MFA,
                  camera_x = 1.8,
                  camera_y = 0,
                  camera_z = 0.2)

side <- plot_MFA(result = Ising_MFA,
                 camera_x = 0,
                 camera_y = -1.8,
                 camera_z = 0.2)

top <- plot_MFA(result = Ising_MFA,
                camera_x = 0,
                camera_y = 0,
                camera_z = 2)

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/misc")

reticulate::use_miniconda('r-reticulate')
reticulate::py_run_string("import sys")
save_image(p = front, file = "ising_mfa_front.pdf", width = 600, height = 600)
save_image(p = side, file = "ising_mfa_side.pdf", width = 600, height = 600)
save_image(p = top, file = "ising_mfa_top.pdf", width = 600, height = 600)

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/mfa_sim_default_df")
Blume_capel_MFA <- read.csv("Blume_capel_default_500x500.csv")

BC_front <- plot_MFA(result = Blume_capel_MFA,
                  camera_x = 1.8,
                  camera_y = 0,
                  camera_z = 0.2,
                  y_label = "Delta")

BC_side <- plot_MFA(result = Blume_capel_MFA,
                 camera_x = 0,
                 camera_y = -1.8,
                 camera_z = 0.2,
                 y_label = "Delta")

BC_top <- plot_MFA(result = Blume_capel_MFA,
                camera_x = 0.001,
                camera_y = 0,
                camera_z = 2,
                y_label = "Delta")

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/misc")
reticulate::use_miniconda('r-reticulate')
reticulate::py_run_string("import sys")
save_image(p = BC_front, file = "BC_mfa_front.pdf", width = 600, height = 600)
save_image(p = BC_side, file = "BC_mfa_side.pdf", width = 600, height = 600)
save_image(p = BC_top, file = "BC_mfa_top.pdf", width = 600, height = 600)


Percolation_MFA <- read.csv("Percolation_default_500x500.csv")
percolation_front <- plot_MFA(result = Percolation_MFA,
                     camera_x = 1.8,
                     camera_y = 0,
                     camera_z = 0.2)

percolation_side <- plot_MFA(result = Percolation_MFA,
                    camera_x = 0,
                    camera_y = -1.8,
                    camera_z = 0.2)

percolation_top <- plot_MFA(result = Percolation_MFA,
                   camera_x = 0,
                   camera_y = 0,
                   camera_z = 2)

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/misc")

reticulate::use_miniconda('r-reticulate')
reticulate::py_run_string("import sys")
save_image(p = percolation_front, file = "percolation_mfa_front.pdf", width = 600, height = 600)
save_image(p = percolation_side, file = "percolation_mfa_side.pdf", width = 600, height = 600)
save_image(p = percolation_top, file = "percolation_mfa_top.pdf", width = 600, height = 600)



continuous01_MFA <- read.csv("continuous01_default_300x300.csv")

front <- plot_MFA(result = continuous01_MFA,
                  camera_x = 1.8,
                  camera_y = 0,
                  camera_z = 0.2)

side <- plot_MFA(result = continuous01_MFA,
                 camera_x = 0,
                 camera_y = -1.8,
                 camera_z = 0.2)

top <- plot_MFA(result = continuous01_MFA,
                camera_x = 0,
                camera_y = 0,
                camera_z = 2)

setwd("/Users/adam/Documents/PhD/stats_mechanic_psych/misc")
reticulate::use_miniconda('r-reticulate')
reticulate::py_run_string("import sys")
save_image(p = front, file = "continuous01_mfa_front.pdf", width = 600, height = 600)
save_image(p = side, file = "continuous01_mfa_side.pdf", width = 600, height = 600)
#save_image(p = continuous01_top, file = "continuous01_mfa_top.pdf", width = 600, height = 600)





#### Is continuous symmetric around alpha = 0
continuous01_output <- MFA_sim(MFA_to_solve = "Continuous_01",
                               alpha_sim = seq(-3.001,1, length.out = 300), 
                               beta_sim = seq(0.0001,100, length.out = 300),
                               average_density = 2)

plot_MFA(result = continuous01_output,
         camera_x = 0,
         camera_y = 0,
         camera_z = 2)


t <- continuous01_output %>% filter(round(alpha,1) == 0)



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



