source("code/sim2_script.R")
library(reshape2)
library(data.table)
library(ggplot2)

get_gamma_transition_props <- function(eps_no_scale_data, tying_grp_vect, n_pheno){
  gamma_and_parent <- matrix(c(eps_no_scale_data[get_parent_indices(1:(n_pheno-1))],eps_no_scale_data)
                             , nrow = 2, ncol = (n_pheno-1)
                             , byrow = T)
  gamma_and_parent_dt <- as.data.table(t(gamma_and_parent))
  setnames(gamma_and_parent_dt,names(gamma_and_parent_dt),c("Parent","Child"))
  gamma_and_parent_dt[,"Transition" := paste0(Child,Parent)]
  gamma_and_parent_dt[,"TreeLvl" := findInterval(.I + 1,tying_grp_vect[-1])]

  # Exclude the root
  gamma_and_parent_dt <- copy(gamma_and_parent_dt[-1])
  gamma_and_parent_stats <- dcast.data.table(
    gamma_and_parent_dt[,.N,by = .(Transition,TreeLvl)]
    ,formula = TreeLvl ~ Transition
    ,value.var = "N"
    ,fill = 0
  )

  transitions_not_there <- setdiff(c("11","01","10","00"),names(gamma_and_parent_stats))
  if(length(transitions_not_there) > 0){
    for(col in transitions_not_there){
      gamma_and_parent_stats[,(col) := 0]
    }
  }

  gamma_and_parent_stats[,"Total" := apply(.SD,1,sum),.SDcols = 2:5]
  return(gamma_and_parent_stats[])
}

for(i in 2:9){
  dat <- readRDS(paste0("data/results_attempt",i,".RDS"))
  assign(x = paste0("results_",i)
         ,dat)
}

results_2$tying_grp <- c(1,2,3,5,9,17,33,65)
results_3$tying_grp <- c(1,2,3,17,33,65)
results_4$tying_grp <- c(1,2,3,17,33,65)
results_5$tying_grp <- c(1,2,3,17,33,65)
results_6$tying_grp <- c(1,2,3,17,33,65)
results_7$tying_grp <- c(1,2,3,33,65)
results_8$tying_grp <- c(1,2,3,33,65)
results_9$tying_grp <- c(1,2,3,33,65)

# Scatterplot of sampled transitions vs algorithm epsilons
sample_vs_algo <- function(results_set){
  list_results <- list()
  for(j in 1:length(results_set$results_gamma_seq)){
    a <- get_gamma_transition_props(eps_no_scale_data = results_set$results_gamma_seq[[j]][-1]
                                    ,tying_grp_vect = results_set$tying_grp
                                    ,n_pheno = 128)
    a[,c("sim_eps_11","sim_eps_10","sim_eps_01","sim_eps_00") :=
        .(`11`/(`11`+`01`)
          ,`10`/(`10`+`00`)
          ,`01`/(`11`+`01`)
          ,`00`/(`10`+`00`))]
    a[,c("algo_eps_11","algo_eps_10") :=
        .(round(c(0,0,results_set$results_eps_11[[j]])[results_set$tying_grp][-(1:2)],2)
          ,round(c(0,0,results_set$results_eps_10[[j]])[results_set$tying_grp][-(1:2)],2))]
    a1 <- copy(a)
    a1 <- a1[,.(TreeLvl
                ,`00`,`01`,`10`,`11`
                ,sim_eps_11,algo_eps_11
                ,sim_eps_10,algo_eps_10
                ,sim_eps_01,algo_eps_01 = 1-algo_eps_11
                ,sim_eps_00,algo_eps_00 = 1-algo_eps_10)]
    a1[,"simNo" := j]
    list_results[[j]] <- a1
  }
  return(rbindlist(list_results))
}

r2 <- sample_vs_algo(results_2)
r3 <- sample_vs_algo(results_3)
r4 <- sample_vs_algo(results_4)
r5 <- sample_vs_algo(results_5)
r6 <- sample_vs_algo(results_6)
r7 <- sample_vs_algo(results_7)
r8 <- sample_vs_algo(results_8)
r9 <- sample_vs_algo(results_9)

create_sample_vs_algo_plots <- function(dataset){
  sim_vs_algo_plot <- list()
  sim_distn_plot <- list()
  algo_distn_plot <- list()
  n = 1
  for(k in unique(dataset$TreeLvl)){
    melted <- melt.data.table(dataset[TreeLvl == k][,-(2:5)]
                              , id.vars = c("TreeLvl","simNo")
                              , measure.vars = patterns("^sim_","^algo_"))
    setnames(melted,c("value1","value2")
             ,c("sim","algo"))
    melted[,"variable" := factor(variable
                                 ,levels = 1:4
                                 ,labels = c("11","10","01","00"))]

    sim_vs_algo_plot[[n]] <- ggplot(melted[variable %in% c("11","10")]) +
      geom_point(aes(x = sim, y = algo, colour = factor(variable)), alpha = 0.5) +
      labs(colour = "transition") +
      ggtitle(paste0("sim vs algo - treeLvl ",k)) +
      scale_x_continuous(breaks = seq(0,1,0.1)) +
      scale_y_continuous(breaks = seq(0,1,0.1))
    sim_distn_plot[[n]] <- ggplot(melted) +
      geom_histogram(aes(x = sim, fill = factor(variable))
                     ,binwidth = 0.05) +
      facet_grid(variable ~ ., scales = "free") +
      guides(fill = FALSE) +
      ggtitle(paste0("simulated transition distribution - treeLvl ",k)) +
      scale_x_continuous(breaks = seq(0,1,0.1))
    algo_distn_plot[[n]] <- ggplot(melted) +
      geom_histogram(aes(x = algo, fill = factor(variable))
                     ,binwidth = 0.05) +
      facet_grid(variable ~ ., scales = "free") +
      guides(fill = FALSE) +
      ggtitle(paste0("algorithm transition distribution - treeLvl ",k)) +
      scale_x_continuous(breaks = seq(0,1,0.1))

    n = n + 1
  }
  return(list(sim_vs_algo_plot = sim_vs_algo_plot
              ,sim_distn_plot = sim_distn_plot
              ,algo_distn_plot = algo_distn_plot))
}

r2_plots <- create_sample_vs_algo_plots(r2)
r3_plots <- create_sample_vs_algo_plots(r3)
r4_plots <- create_sample_vs_algo_plots(r4)
r5_plots <- create_sample_vs_algo_plots(r5)
r6_plots <- create_sample_vs_algo_plots(r6)
r7_plots <- create_sample_vs_algo_plots(r7)
r8_plots <- create_sample_vs_algo_plots(r8)
r9_plots <- create_sample_vs_algo_plots(r9)
