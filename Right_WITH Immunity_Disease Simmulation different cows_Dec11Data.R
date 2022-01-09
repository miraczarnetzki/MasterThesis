library(igraph)

#loading our data of Network
load("/Users/miraczarnetzki/OneDrive - Högskolan Dalarna/Thesis/Disease Simulation/Data 11 Dec/network_1.25m_30min_right.Rda")

summary(igraph_df_30)

#making adjacency matrix
from_cow <- igraph_df_30[, c('from_ID')]
to_cow <- igraph_df_30[, c('to_ID')]
dat <- data.frame(from_cow, to_cow)
adj_matrix <- get.adjacency(graph.edgelist(as.matrix(dat), directed=FALSE))
adj_matrix
#column_names <- data.frame(colnames(adj_matrix)) - not needed anymore
#column_names 


#colnames(adj_matrix) <- c(1:dim(adj_matrix)[1]) - not needed anymore
#rownames(adj_matrix) <- c(1:dim(adj_matrix)[1])
adj_matrix
class(adj_matrix)
dim(adj_matrix)
adj_matrix[1,1]

#simulation function with immunity
sim.network.transmission.immunity <- function(adj_matrix, first.infected = 1, immune = NULL, prob.infected = 0.03, prob.cured = 0.1, n.weeks = 50) {
  all.sims <- list()
  df <- data.frame(first.infected,1) #Infection starts with individual 1
  names(df) = c("infected","weeks")
  t=0
  while(nrow(df) >0 & t<n.weeks) {
    t=t+1
    remove=NULL
    #First check who are cured and then who are infected
    ## Which are cured
    for(i in 1:nrow(df)) {
      cured.sim <- rbinom(1,1,prob.cured)
      if (cured.sim == 1) remove <- c(remove, i)
    }
    if(!is.null(remove)) df <- df[-remove,]
    if(nrow(df) > 0) {
      ## Which are infected
      for(i in 1:nrow(df)) {
        n.edges <- sum(adj_matrix[df[i,1],])
        if(n.edges != 0 ) {
          for(j in 1:n.edges) {
            indx <- which(adj_matrix[df[i,1],]==1)
            ind <- indx[j]
            not.infected <- !( ind  %in%  df[,1] ) 
            check.immune = FALSE
            if (!is.null(immune)) check.immune <- indx[j]==immune
            if(not.infected & !check.immune) {
              infected.sim <- rbinom(1,1,prob.infected)
              if(infected.sim == 1) df <- rbind(df, c(ind, 1))
            }
          }
        }
        df[i,2] <- df[i,2]+1
      }
    }
    cat("Week", t, "\n")
    rownames(df)=NULL
    ifelse(nrow(df) >0, print(df), print("No infected"))
    all.sims[[t]] <- df
  }
  return(all.sims)
}

####### Simulating disease transmission with immunity on different cows #######

##### Simulations on Cow 136 #########
#COW 136 - Central to community 1

### 1. scenario cow 139 immune - one with highest degree centrality and eigen vector centrality (most amount of connections)####
#Start infection Cow 136 -> Simulating 500 times disease with immunity for cow 139
simu136.i139 <- list()
for (i in 1:500){
  simu136.i139[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=136, immune = 139) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week136.i139 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu136.i139[[i]])
  for (j in 1: N_weeks){
    cows.per.week136.i139[i,j] <- length(simu136.i139[[i]][[j]][[1]])
    cows.per.week136.i139[is.na(cows.per.week136.i139)] <- 0
    print(cows.per.week136.i139[i,j])
  }
}

View(cows.per.week136.i139)

#at a prob.infect = 0.03 in 291 cases it dies out 

#calculating the average
mean.infected.per.day136.i139 <- as.data.frame(colMeans(cows.per.week136.i139))

View(mean.infected.per.day136.i139)

### Next scenario cow 112 immune - highest betweeness centrality #####
#Start infection Cow 27 -> Simulating 500 times disease with immunity for cow 9
simu136.i112 <- list()
for (i in 1:500){
  simu136.i112[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=136, immune = 112) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week136.i112 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu136.i112[[i]])
  for (j in 1: N_weeks){
    cows.per.week136.i112[i,j] <- length(simu136.i112[[i]][[j]][[1]])
    cows.per.week136.i112[is.na(cows.per.week136.i112)] <- 0
    print(cows.per.week136.i112[i,j])
  }
}

View(cows.per.week136.i112)

#at a prob.infect = 0.03 in 289 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day136.i112 <- as.data.frame(colMeans(cows.per.week136.i112))

View(mean.infected.per.day136.i112)

#### Next Scenario Cow 92 out of smallest community immune #####
#Start infection Cow 136 -> Simulating 500 times disease with immunity for cow 92
simu136.i92 <- list()
for (i in 1:500){
  simu136.i92[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=136, immune = 92) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week136.i92 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu136.i92[[i]])
  for (j in 1: N_weeks){
    cows.per.week136.i92[i,j] <- length(simu136.i92[[i]][[j]][[1]])
    cows.per.week136.i92[is.na(cows.per.week136.i92)] <- 0
    print(cows.per.week136.i92[i,j])
  }
}

View(cows.per.week136.i92)

#at a prob.infect = 0.03 in 282 cases it dies out 

#calculating the average
mean.infected.per.day136.i92 <- as.data.frame(colMeans(cows.per.week136.i92))

View(mean.infected.per.day136.i92)


#### Next Scenario Cow 158 out of second largest community immune #####
#Start infection Cow 136 -> Simulating 500 times disease with immunity for cow 158
simu136.i158 <- list()
for (i in 1:500){
  simu136.i158[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=136, immune = 158) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week136.i158 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu136.i158[[i]])
  for (j in 1: N_weeks){
    cows.per.week136.i158[i,j] <- length(simu136.i158[[i]][[j]][[1]])
    cows.per.week136.i158[is.na(cows.per.week136.i158)] <- 0
    print(cows.per.week136.i158[i,j])
  }
}

View(cows.per.week136.i158)

#at a prob.infect = 0.03 in 304 cases it dies out compared 

#calculating the average
mean.infected.per.day136.i158 <- as.data.frame(colMeans(cows.per.week136.i158))

View(mean.infected.per.day136.i158)

#### End of Simulation on Cow 136 #######


##### Simulations on Cow 158 #########
#COW 158 - Central to community 2

### 1. scenario cow 139 immune - one with highest higher degree centrality and eigen vector centrality (most amount of connections)####
#Start infection Cow 158 -> Simulating 500 times disease with immunity for cow 139
simu158.i139 <- list()
for (i in 1:500){
  simu158.i139[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=158, immune = 139) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week158.i139 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu158.i139[[i]])
  for (j in 1: N_weeks){
    cows.per.week158.i139[i,j] <- length(simu158.i139[[i]][[j]][[1]])
    cows.per.week158.i139[is.na(cows.per.week158.i139)] <- 0
    print(cows.per.week158.i139[i,j])
  }
}

View(cows.per.week158.i139)

#at a prob.infect = 0.03 in 291 cases it dies out 

#calculating the average
mean.infected.per.day158.i139 <- as.data.frame(colMeans(cows.per.week158.i139))

View(mean.infected.per.day158.i139)


### Next scenario cow 112 immune - highest betweeness centrality #####
#Start infection Cow 158 -> Simulating 500 times disease with immunity for cow 9
simu158.i112 <- list()
for (i in 1:500){
  simu158.i112[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=158, immune = 112) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week158.i112 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu158.i112[[i]])
  for (j in 1: N_weeks){
    cows.per.week158.i112[i,j] <- length(simu158.i112[[i]][[j]][[1]])
    cows.per.week158.i112[is.na(cows.per.week158.i112)] <- 0
    print(cows.per.week158.i112[i,j])
  }
}

View(cows.per.week158.i112)

#at a prob.infect = 0.03 in 289 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day158.i112 <- as.data.frame(colMeans(cows.per.week158.i112))

View(mean.infected.per.day158.i112)

#### Next Scenario Cow 92 out of smallest community immune #####
#Start infection Cow 158 -> Simulating 500 times disease with immunity for cow 92
simu158.i92 <- list()
for (i in 1:500){
  simu158.i92[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=158, immune = 92) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week158.i92 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu158.i92[[i]])
  for (j in 1: N_weeks){
    cows.per.week158.i92[i,j] <- length(simu158.i92[[i]][[j]][[1]])
    cows.per.week158.i92[is.na(cows.per.week158.i92)] <- 0
    print(cows.per.week158.i92[i,j])
  }
}

View(cows.per.week158.i92)

#at a prob.infect = 0.03 in 282 cases it dies out 

#calculating the average
mean.infected.per.day158.i92 <- as.data.frame(colMeans(cows.per.week158.i92))

View(mean.infected.per.day158.i92)


#### Next Scenario Cow 136 out of second largest community immune #####
#Start infection Cow 158 -> Simulating 500 times disease with immunity for cow 136
simu158.i136 <- list()
for (i in 1:500){
  simu158.i136[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=158, immune = 136) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week158.i136 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu158.i136[[i]])
  for (j in 1: N_weeks){
    cows.per.week158.i136[i,j] <- length(simu158.i136[[i]][[j]][[1]])
    cows.per.week158.i136[is.na(cows.per.week158.i136)] <- 0
    print(cows.per.week158.i136[i,j])
  }
}

View(cows.per.week158.i136)

#at a prob.infect = 0.03 in 304 cases it dies out compared 

#calculating the average
mean.infected.per.day158.i136 <- as.data.frame(colMeans(cows.per.week158.i136))

View(mean.infected.per.day158.i136)

#### End of Simulation on Cow 158 #######


##### Simulations on Cow 92 #########
#COW 92 - Central to community 3

### 1. scenario cow 139 immune - one with highest higher degree centrality and eigen vector centrality (most amount of connections)####
#Start infection Cow 92 -> Simulating 500 times disease with immunity for cow 139
simu92.i139 <- list()
for (i in 1:500){
  simu92.i139[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=92, immune = 139) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week92.i139 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu92.i139[[i]])
  for (j in 1: N_weeks){
    cows.per.week92.i139[i,j] <- length(simu92.i139[[i]][[j]][[1]])
    cows.per.week92.i139[is.na(cows.per.week92.i139)] <- 0
    print(cows.per.week92.i139[i,j])
  }
}

View(cows.per.week92.i139)

#at a prob.infect = 0.03 in 291 cases it dies out 

#calculating the average
mean.infected.per.day92.i139 <- as.data.frame(colMeans(cows.per.week92.i139))

View(mean.infected.per.day92.i139)


### Next scenario cow 112 immune - highest betweeness centrality #####
#Start infection Cow 92 -> Simulating 500 times disease with immunity for cow 112
simu92.i112 <- list()
for (i in 1:500){
  simu92.i112[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=92, immune = 112) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week92.i112 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu92.i112[[i]])
  for (j in 1: N_weeks){
    cows.per.week92.i112[i,j] <- length(simu92.i112[[i]][[j]][[1]])
    cows.per.week92.i112[is.na(cows.per.week92.i112)] <- 0
    print(cows.per.week92.i112[i,j])
  }
}

View(cows.per.week92.i112)

#at a prob.infect = 0.03 in 289 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day92.i112 <- as.data.frame(colMeans(cows.per.week92.i112))

View(mean.infected.per.day92.i112)

#### Next Scenario Cow 116 out of second smallest community immune #####
#Start infection Cow 92 -> Simulating 500 times disease with immunity for cow 116
simu92.i116 <- list()
for (i in 1:500){
  simu92.i116[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=92, immune = 116) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week92.i116 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu92.i116[[i]])
  for (j in 1: N_weeks){
    cows.per.week92.i116[i,j] <- length(simu92.i116[[i]][[j]][[1]])
    cows.per.week92.i116[is.na(cows.per.week92.i116)] <- 0
    print(cows.per.week92.i116[i,j])
  }
}

View(cows.per.week92.i116)

#at a prob.infect = 0.03 in 282 cases it dies out 

#calculating the average
mean.infected.per.day92.i116 <- as.data.frame(colMeans(cows.per.week92.i116))

View(mean.infected.per.day92.i116)


#### Next Scenario Cow 136 out of largest community immune #####
#Start infection Cow 92 -> Simulating 500 times disease with immunity for cow 136
simu92.i136 <- list()
for (i in 1:500){
  simu92.i136[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=92, immune = 136) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week92.i136 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu92.i136[[i]])
  for (j in 1: N_weeks){
    cows.per.week92.i136[i,j] <- length(simu92.i136[[i]][[j]][[1]])
    cows.per.week92.i136[is.na(cows.per.week92.i136)] <- 0
    print(cows.per.week92.i136[i,j])
  }
}

View(cows.per.week92.i136)

#at a prob.infect = 0.03 in 304 cases it dies out compared 

#calculating the average
mean.infected.per.day92.i136 <- as.data.frame(colMeans(cows.per.week92.i136))

View(mean.infected.per.day92.i136)

#### End of Simulation on Cow 92 #######


##### Simulations on Cow 116 #########
#COW 116 - Central to community 4

### 1. scenario cow 139 immune - one with highest higher degree centrality and eigen vector centrality (most amount of connections)####
#Start infection Cow 116 -> Simulating 500 times disease with immunity for cow 139
simu116.i139 <- list()
for (i in 1:500){
  simu116.i139[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=116, immune = 139) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week116.i139 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu116.i139[[i]])
  for (j in 1: N_weeks){
    cows.per.week116.i139[i,j] <- length(simu116.i139[[i]][[j]][[1]])
    cows.per.week116.i139[is.na(cows.per.week116.i139)] <- 0
    print(cows.per.week116.i139[i,j])
  }
}

View(cows.per.week116.i139)

#at a prob.infect = 0.03 in 291 cases it dies out 

#calculating the average
mean.infected.per.day116.i139 <- as.data.frame(colMeans(cows.per.week116.i139))

View(mean.infected.per.day116.i139)


### Next scenario cow 112 immune - highest betweeness centrality #####
#Start infection Cow 116 -> Simulating 500 times disease with immunity for cow 112
simu116.i112 <- list()
for (i in 1:500){
  simu116.i112[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=116, immune = 112) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week116.i112 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu116.i112[[i]])
  for (j in 1: N_weeks){
    cows.per.week116.i112[i,j] <- length(simu116.i112[[i]][[j]][[1]])
    cows.per.week116.i112[is.na(cows.per.week116.i112)] <- 0
    print(cows.per.week116.i112[i,j])
  }
}

View(cows.per.week116.i112)

#at a prob.infect = 0.03 in 289 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day116.i112 <- as.data.frame(colMeans(cows.per.week116.i112))

View(mean.infected.per.day116.i112)

#### Next Scenario Cow 92 out of smallest community immune #####
#Start infection Cow 116 -> Simulating 500 times disease with immunity for cow 92
simu116.i92 <- list()
for (i in 1:500){
  simu116.i92[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=116, immune = 92) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week116.i92 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu116.i92[[i]])
  for (j in 1: N_weeks){
    cows.per.week116.i92[i,j] <- length(simu116.i92[[i]][[j]][[1]])
    cows.per.week116.i92[is.na(cows.per.week116.i92)] <- 0
    print(cows.per.week116.i92[i,j])
  }
}

View(cows.per.week116.i92)

#at a prob.infect = 0.03 in 282 cases it dies out 

#calculating the average
mean.infected.per.day116.i92 <- as.data.frame(colMeans(cows.per.week116.i92))

View(mean.infected.per.day116.i92)


#### Next Scenario Cow 136 out of second largest community immune #####
#Start infection Cow 116 -> Simulating 500 times disease with immunity for cow 158
simu116.i136 <- list()
for (i in 1:500){
  simu116.i136[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=116, immune = 136) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week116.i136 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu116.i136[[i]])
  for (j in 1: N_weeks){
    cows.per.week116.i136[i,j] <- length(simu116.i136[[i]][[j]][[1]])
    cows.per.week116.i136[is.na(cows.per.week116.i136)] <- 0
    print(cows.per.week116.i136[i,j])
  }
}

View(cows.per.week116.i136)

#at a prob.infect = 0.03 in 304 cases it dies out compared 

#calculating the average
mean.infected.per.day116.i136 <- as.data.frame(colMeans(cows.per.week116.i136))

View(mean.infected.per.day116.i136)

#### End of Simulation on Cow 116 #######


##### Simulations on Cow 84 #########
#COW 84 - Central to community 5

### 1. scenario cow 139 immune - one with highest higher degree centrality and eigen vector centrality (most amount of connections)####
#Start infection Cow 84 -> Simulating 500 times disease with immunity for cow 139
simu84.i139 <- list()
for (i in 1:500){
  simu84.i139[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=84, immune = 139) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week84.i139 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu84.i139[[i]])
  for (j in 1: N_weeks){
    cows.per.week84.i139[i,j] <- length(simu84.i139[[i]][[j]][[1]])
    cows.per.week84.i139[is.na(cows.per.week84.i139)] <- 0
    print(cows.per.week84.i139[i,j])
  }
}

View(cows.per.week84.i139)

#at a prob.infect = 0.03 in 291 cases it dies out 

#calculating the average
mean.infected.per.day84.i139 <- as.data.frame(colMeans(cows.per.week84.i139))

View(mean.infected.per.day84.i139)


### Next scenario cow 112 immune - highest betweeness centrality #####
#Start infection Cow 84 -> Simulating 500 times disease with immunity for cow 112
simu84.i112 <- list()
for (i in 1:500){
  simu84.i112[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=84, immune = 112) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week84.i112 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu84.i112[[i]])
  for (j in 1: N_weeks){
    cows.per.week84.i112[i,j] <- length(simu84.i112[[i]][[j]][[1]])
    cows.per.week84.i112[is.na(cows.per.week84.i112)] <- 0
    print(cows.per.week84.i112[i,j])
  }
}

View(cows.per.week84.i112)

#at a prob.infect = 0.03 in 289 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day84.i112 <- as.data.frame(colMeans(cows.per.week84.i112))

View(mean.infected.per.day84.i112)

#### Next Scenario Cow 92 out of smallest community immune #####
#Start infection Cow 84 -> Simulating 500 times disease with immunity for cow 92
simu84.i92 <- list()
for (i in 1:500){
  simu84.i92[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=84, immune = 92) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week84.i92 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu84.i92[[i]])
  for (j in 1: N_weeks){
    cows.per.week84.i92[i,j] <- length(simu84.i92[[i]][[j]][[1]])
    cows.per.week84.i92[is.na(cows.per.week84.i92)] <- 0
    print(cows.per.week84.i92[i,j])
  }
}

View(cows.per.week84.i92)

#at a prob.infect = 0.03 in 282 cases it dies out 

#calculating the average
mean.infected.per.day84.i92 <- as.data.frame(colMeans(cows.per.week84.i92))

View(mean.infected.per.day84.i92)


#### Next Scenario Cow 136 out of largest community immune #####
#Start infection Cow 84 -> Simulating 500 times disease with immunity for cow 136
simu84.i136 <- list()
for (i in 1:500){
  simu84.i136[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=84, immune = 136) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week84.i136 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu84.i136[[i]])
  for (j in 1: N_weeks){
    cows.per.week84.i136[i,j] <- length(simu84.i136[[i]][[j]][[1]])
    cows.per.week84.i136[is.na(cows.per.week84.i136)] <- 0
    print(cows.per.week84.i136[i,j])
  }
}

View(cows.per.week84.i136)

#at a prob.infect = 0.03 in 304 cases it dies out compared 

#calculating the average
mean.infected.per.day84.i136 <- as.data.frame(colMeans(cows.per.week84.i136))

View(mean.infected.per.day84.i136)

#### End of Simulation on Cow 84 #######


### END OF SIMULATIONS #####

##### Plotting #######

#Combining all the cows mean data together
days <- as.data.frame(1:50)
mean.infected.per.day.immunity.total.right <- cbind(days, mean.infected.per.day136.i139, mean.infected.per.day136.i112, mean.infected.per.day136.i158, mean.infected.per.day136.i92,
                                                   mean.infected.per.day158.i139, mean.infected.per.day158.i112, mean.infected.per.day158.i136, mean.infected.per.day158.i92, 
                                                   mean.infected.per.day92.i139, mean.infected.per.day92.i112, mean.infected.per.day92.i136, mean.infected.per.day92.i116,
                                                   mean.infected.per.day116.i139, mean.infected.per.day116.i112, mean.infected.per.day116.i136, mean.infected.per.day116.i92,
                                                   mean.infected.per.day84.i139, mean.infected.per.day84.i112, mean.infected.per.day84.i136, mean.infected.per.day84.i92)

View(mean.infected.per.day.immunity.total.left)

#Plotting & Analyzing the Means of all the Simulations on different Cows with Immunity



#plotting different immunity scenarios vs. no immunity on individual cows

#Cow 136
plot(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week136.i139)`, type = "l", col = 2, lty = 1, ylim = c(0, 22), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 136", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week136.i112)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week136.i158)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week136.i92)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week136)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 136 First Infected with Immunity Cow 139", "Cow 136 First Infected with Immunity Cow 112", "Cow 136 First Infected with Immunity Cow 158", 
                "Cow 136 First Infected with Immunity Cow 92", "Cow 136 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 158
plot(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week158.i139)`, type = "l", col = 2, lty = 1, ylim = c(0, 25), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 158", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week158.i112)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week158.i136)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week158.i92)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week158)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 158 First Infected with Immunity Cow 139", "Cow 158 First Infected with Immunity Cow 112", "Cow 158 First Infected with Immunity Cow 136", 
                "Cow 158 First Infected with Immunity Cow 92", "Cow 158 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 92
plot(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week92.i139)`, type = "l", col = 2, lty = 1, ylim = c(0, 25), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 92", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week92.i112)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week92.i136)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week92.i116)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week92)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 92 First Infected with Immunity Cow 139", "Cow 92 First Infected with Immunity Cow 112", "Cow 92 First Infected with Immunity Cow 136", 
                "Cow 92 First Infected with Immunity Cow 116", "Cow 92 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 116
plot(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week116.i139)`, type = "l", col = 2, lty = 1, ylim = c(0, 22), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 116", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week116.i112)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week116.i136)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week116.i92)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week116)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 116 First Infected with Immunity Cow 139", "Cow 116 First Infected with Immunity Cow 112", "Cow 116 First Infected with Immunity Cow 136", 
                "Cow 116 First Infected with Immunity Cow 92", "Cow 116 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 84
plot(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week84.i139)`, type = "l", col = 2, lty = 1, ylim = c(0, 25), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 84", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week84.i112)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week84.i136)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.immunity.total.right$`colMeans(cows.per.week84.i92)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week84)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 84 First Infected with Immunity Cow 139", "Cow 84 First Infected with Immunity Cow 112", "Cow 84 First Infected with Immunity Cow 136", 
                "Cow 84 First Infected with Immunity Cow 92", "Cow 84 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)


