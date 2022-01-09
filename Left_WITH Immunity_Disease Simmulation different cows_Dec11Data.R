library(igraph)

#loading our data of Network
load("/Users/miraczarnetzki/OneDrive - Högskolan Dalarna/Thesis/Disease Simulation/Data 11 Dec/network_1.25m_30min_left.Rda")

summary(igraph_df_30)

#making adjacency matrix
from_cow <- igraph_df_30[, c('from')]
to_cow <- igraph_df_30[, c('to')]
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

##### Simulations on Cow 27 #########
#COW 27 - Central to community 1

### 1. scenario cow 31 immune ####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 31
simu27.i31 <- list()
for (i in 1:500){
  simu27.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=27, immune = 31) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week27.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu27.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week27.i31[i,j] <- length(simu27.i31[[i]][[j]][[1]])
    cows.per.week27.i31[is.na(cows.per.week27.i31)] <- 0
    print(cows.per.week27.i31[i,j])
  }
}

View(cows.per.week27.i31)

#at a prob.infect = 0.03 in 291 cases it dies out 

#calculating the average
mean.infected.per.day27.i31 <- as.data.frame(colMeans(cows.per.week27.i31))

View(mean.infected.per.day27.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 27 -> Simulating 500 times disease with immunity for cow 9
simu27.i9 <- list()
for (i in 1:500){
  simu27.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=27, immune = 9) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week27.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu27.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week27.i9[i,j] <- length(simu27.i9[[i]][[j]][[1]])
    cows.per.week27.i9[is.na(cows.per.week27.i9)] <- 0
    print(cows.per.week27.i9[i,j])
  }
}

View(cows.per.week27.i9)

#at a prob.infect = 0.03 in 289 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day27.i9 <- as.data.frame(colMeans(cows.per.week27.i9))

View(mean.infected.per.day27.i9)

#### Next Scenario Cow 63 out of biggest community immune #####
#Start infection Cow 27 -> Simulating 500 times disease with immunity for cow 63
simu27.i63 <- list()
for (i in 1:500){
  simu27.i63[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=27, immune = 63) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week27.i63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu27.i63[[i]])
  for (j in 1: N_weeks){
    cows.per.week27.i63[i,j] <- length(simu27.i63[[i]][[j]][[1]])
    cows.per.week27.i63[is.na(cows.per.week27.i63)] <- 0
    print(cows.per.week27.i63[i,j])
  }
}

View(cows.per.week27.i63)

#at a prob.infect = 0.03 in 282 cases it dies out 

#calculating the average
mean.infected.per.day27.i63 <- as.data.frame(colMeans(cows.per.week27.i63))

View(mean.infected.per.day27.i63)


#### Next Scenario Cow 40 out of smallest community immune #####
#Start infection Cow 27 -> Simulating 500 times disease with immunity for cow 40
simu27.i40 <- list()
for (i in 1:500){
  simu27.i40[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=27, immune = 40) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week27.i40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu27.i40[[i]])
  for (j in 1: N_weeks){
    cows.per.week27.i40[i,j] <- length(simu27.i40[[i]][[j]][[1]])
    cows.per.week27.i40[is.na(cows.per.week27.i40)] <- 0
    print(cows.per.week27.i40[i,j])
  }
}

View(cows.per.week27.i40)

#at a prob.infect = 0.03 in 304 cases it dies out compared 

#calculating the average
mean.infected.per.day27.i40 <- as.data.frame(colMeans(cows.per.week27.i40))

View(mean.infected.per.day27.i40)

#### End of Simulation on Cow 27 #######


##### Simulations on Cow 40 #########
#COW 40 - Central to community 2

### 1 scenario cow 31 immune ####
#Start infection Cow 40 -> Simulating 500 times disease with immunity for cow 31
simu40.i31 <- list()
for (i in 1:500){
  simu40.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=40, immune = 31) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week40.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu40.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week40.i31[i,j] <- length(simu40.i31[[i]][[j]][[1]])
    cows.per.week40.i31[is.na(cows.per.week40.i31)] <- 0
    print(cows.per.week40.i31[i,j])
  }
}

View(cows.per.week40.i31)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day40.i31 <- as.data.frame(colMeans(cows.per.week40.i31))

View(mean.infected.per.day40.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 40 -> Simulating 500 times disease with immunity for cow 9
simu40.i9 <- list()
for (i in 1:500){
  simu40.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=40, immune = 9) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week40.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu40.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week40.i9[i,j] <- length(simu40.i9[[i]][[j]][[1]])
    cows.per.week40.i9[is.na(cows.per.week40.i9)] <- 0
    print(cows.per.week40.i9[i,j])
  }
}

View(cows.per.week40.i9)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day40.i9 <- as.data.frame(colMeans(cows.per.week40.i9))

View(mean.infected.per.day40.i9)

#### Next Scenario Cow 63 out of biggest community immune #####
#Start infection Cow 40 -> Simulating 500 times disease with immunity for cow 63
simu40.i63 <- list()
for (i in 1:500){
  simu40.i63[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=40, immune = 63) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week40.i63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu40.i63[[i]])
  for (j in 1: N_weeks){
    cows.per.week40.i63[i,j] <- length(simu40.i63[[i]][[j]][[1]])
    cows.per.week40.i63[is.na(cows.per.week40.i63)] <- 0
    print(cows.per.week40.i63[i,j])
  }
}

View(cows.per.week40.i63)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day40.i63 <- as.data.frame(colMeans(cows.per.week40.i63))

View(mean.infected.per.day40.i63)


#### Next Scenario Cow 1 out of second smallest community & important in cliques immune #####
#Start infection Cow 40 -> Simulating 500 times disease with immunity for cow 1
simu40.i1 <- list()
for (i in 1:500){
  simu40.i1[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=40, immune = 1) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week40.i1 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu40.i1[[i]])
  for (j in 1: N_weeks){
    cows.per.week40.i1[i,j] <- length(simu40.i1[[i]][[j]][[1]])
    cows.per.week40.i1[is.na(cows.per.week40.i1)] <- 0
    print(cows.per.week40.i1[i,j])
  }
}

View(cows.per.week40.i1)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day40.i1 <- as.data.frame(colMeans(cows.per.week40.i1))

View(mean.infected.per.day40.i1)

#### End of Simulation on Cow 40 #######

##### Simulations on Cow 63 #########
#COW 63 - Central to community 3

### 1 scenario cow 31 immune ####
#Start infection Cow 63 -> Simulating 500 times disease with immunity for cow 31
simu63.i31 <- list()
for (i in 1:500){
  simu63.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=63, immune = 31) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week63.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu63.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week63.i31[i,j] <- length(simu63.i31[[i]][[j]][[1]])
    cows.per.week63.i31[is.na(cows.per.week63.i31)] <- 0
    print(cows.per.week63.i31[i,j])
  }
}

View(cows.per.week63.i31)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day63.i31 <- as.data.frame(colMeans(cows.per.week63.i31))

View(mean.infected.per.day63.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 63 -> Simulating 500 times disease with immunity for cow 9
simu63.i9 <- list()
for (i in 1:500){
  simu63.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=63, immune = 9) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week63.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu63.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week63.i9[i,j] <- length(simu63.i9[[i]][[j]][[1]])
    cows.per.week63.i9[is.na(cows.per.week63.i9)] <- 0
    print(cows.per.week63.i9[i,j])
  }
}

View(cows.per.week63.i9)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day63.i9 <- as.data.frame(colMeans(cows.per.week63.i9))

View(mean.infected.per.day63.i9)

#### Next Scenario Cow 15 out of second biggest community immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 63
simu63.i15 <- list()
for (i in 1:500){
  simu63.i15[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=63, immune = 15) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week63.i15 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu63.i15[[i]])
  for (j in 1: N_weeks){
    cows.per.week63.i15[i,j] <- length(simu63.i15[[i]][[j]][[1]])
    cows.per.week63.i15[is.na(cows.per.week63.i15)] <- 0
    print(cows.per.week63.i15[i,j])
  }
}

View(cows.per.week63.i15)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day63.i15 <- as.data.frame(colMeans(cows.per.week63.i15))

View(mean.infected.per.day63.i15)


#### Next Scenario Cow 40 out of smallest community immune #####
#Start infection Cow 63 -> Simulating 500 times disease with immunity for cow 40
simu63.i40 <- list()
for (i in 1:500){
  simu63.i40[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=63, immune = 40) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week63.i40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu63.i40[[i]])
  for (j in 1: N_weeks){
    cows.per.week63.i40[i,j] <- length(simu63.i40[[i]][[j]][[1]])
    cows.per.week63.i40[is.na(cows.per.week63.i40)] <- 0
    print(cows.per.week63.i40[i,j])
  }
}

View(cows.per.week63.i40)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day63.i40 <- as.data.frame(colMeans(cows.per.week63.i40))

View(mean.infected.per.day63.i40)

#### End of Simulation on Cow 63 #######


##### Simulations on Cow 15 #########
#COW 15 - Central to community 4

### 1 scenario cow 31 immune ####
#Start infection Cow 15 -> Simulating 500 times disease with immunity for cow 31
simu15.i31 <- list()
for (i in 1:500){
  simu15.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=15, immune = 31) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week15.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu15.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week15.i31[i,j] <- length(simu15.i31[[i]][[j]][[1]])
    cows.per.week15.i31[is.na(cows.per.week15.i31)] <- 0
    print(cows.per.week15.i31[i,j])
  }
}

View(cows.per.week15.i31)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day15.i31 <- as.data.frame(colMeans(cows.per.week15.i31))

View(mean.infected.per.day15.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 15 -> Simulating 500 times disease with immunity for cow 9
simu15.i9 <- list()
for (i in 1:500){
  simu15.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=15, immune = 9) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week15.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu15.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week15.i9[i,j] <- length(simu15.i9[[i]][[j]][[1]])
    cows.per.week15.i9[is.na(cows.per.week15.i9)] <- 0
    print(cows.per.week15.i9[i,j])
  }
}

View(cows.per.week15.i9)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day15.i9 <- as.data.frame(colMeans(cows.per.week15.i9))

View(mean.infected.per.day15.i9)

### Next scenario cow 63 out of biggest community immune #####
#Start infection Cow 15 -> Simulating 500 times disease with immunity for cow 63
simu15.i63 <- list()
for (i in 1:500){
  simu15.i63[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=15, immune = 63) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week15.i63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu15.i63[[i]])
  for (j in 1: N_weeks){
    cows.per.week15.i63[i,j] <- length(simu15.i63[[i]][[j]][[1]])
    cows.per.week15.i63[is.na(cows.per.week15.i63)] <- 0
    print(cows.per.week15.i63[i,j])
  }
}

View(cows.per.week15.i63)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day15.i63 <- as.data.frame(colMeans(cows.per.week15.i63))

View(mean.infected.per.day15.i63)


#### Next Scenario Cow 40 out of smallest community immune #####
#Start infection Cow 15 -> Simulating 500 times disease with immunity for cow 40
simu15.i40 <- list()
for (i in 1:500){
  simu15.i40[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=15, immune = 40) #need to specify which cow gets immune by immune = 
}

#Adding NA ifstatement
cows.per.week15.i40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu15.i40[[i]])
  for (j in 1: N_weeks){
    cows.per.week15.i40[i,j] <- length(simu15.i40[[i]][[j]][[1]])
    cows.per.week15.i40[is.na(cows.per.week15.i40)] <- 0
    print(cows.per.week15.i40[i,j])
  }
}

View(cows.per.week15.i40)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day15.i40 <- as.data.frame(colMeans(cows.per.week15.i40))

View(mean.infected.per.day15.i40)

#### End of Simulation on Cow 15 #######


##### Simulations on Cow 1 #########
#COW 1 - Central to community 5

### 1 scenario cow 31 immune ####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 31
simu1.i31 <- list()
for (i in 1:500){
  simu1.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 31) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i31[i,j] <- length(simu1.i31[[i]][[j]][[1]])
    cows.per.week1.i31[is.na(cows.per.week1.i31)] <- 0
    print(cows.per.week1.i31[i,j])
  }
}

View(cows.per.week1.i31)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i31 <- as.data.frame(colMeans(cows.per.week1.i31))

View(mean.infected.per.day1.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 9
simu1.i9 <- list()
for (i in 1:500){
  simu1.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 9) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i9[i,j] <- length(simu1.i9[[i]][[j]][[1]])
    cows.per.week1.i9[is.na(cows.per.week1.i9)] <- 0
    print(cows.per.week1.i9[i,j])
  }
}

View(cows.per.week1.i9)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i9 <- as.data.frame(colMeans(cows.per.week1.i9))

View(mean.infected.per.day1.i9)

#### Next Scenario Cow 63 out of biggest community immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 63
simu1.i63 <- list()
for (i in 1:500){
  simu1.i63[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 63) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i63[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i63[i,j] <- length(simu1.i63[[i]][[j]][[1]])
    cows.per.week1.i63[is.na(cows.per.week1.i63)] <- 0
    print(cows.per.week1.i63[i,j])
  }
}

View(cows.per.week1.i63)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i63 <- as.data.frame(colMeans(cows.per.week1.i63))

View(mean.infected.per.day1.i63)


#### Next Scenario Cow 40 out of smallest community immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 40
simu1.i40 <- list()
for (i in 1:500){
  simu1.i40[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 40) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i40[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i40[i,j] <- length(simu1.i40[[i]][[j]][[1]])
    cows.per.week1.i40[is.na(cows.per.week1.i40)] <- 0
    print(cows.per.week1.i40[i,j])
  }
}

View(cows.per.week1.i40)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i40 <- as.data.frame(colMeans(cows.per.week1.i40))

View(mean.infected.per.day1.i40)

#### End of Simulation on Cow 1 #######


##### Simulations on Cow 1 #########
#COW 1 - Central to community 5

### 1 scenario cow 31 immune ####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 31
simu1.i31 <- list()
for (i in 1:500){
  simu1.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 31) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i31[i,j] <- length(simu1.i31[[i]][[j]][[1]])
    cows.per.week1.i31[is.na(cows.per.week1.i31)] <- 0
    print(cows.per.week1.i31[i,j])
  }
}

View(cows.per.week1.i31)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i31 <- as.data.frame(colMeans(cows.per.week1.i31))

View(mean.infected.per.day1.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 9
simu1.i9 <- list()
for (i in 1:500){
  simu1.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 9) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i9[i,j] <- length(simu1.i9[[i]][[j]][[1]])
    cows.per.week1.i9[is.na(cows.per.week1.i9)] <- 0
    print(cows.per.week1.i9[i,j])
  }
}

View(cows.per.week1.i9)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i9 <- as.data.frame(colMeans(cows.per.week1.i9))

View(mean.infected.per.day1.i9)

#### Next Scenario Cow 63 out of biggest community immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 63
simu1.i63 <- list()
for (i in 1:500){
  simu1.i63[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 63) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i63[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i63[i,j] <- length(simu1.i63[[i]][[j]][[1]])
    cows.per.week1.i63[is.na(cows.per.week1.i63)] <- 0
    print(cows.per.week1.i63[i,j])
  }
}

View(cows.per.week1.i63)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i63 <- as.data.frame(colMeans(cows.per.week1.i63))

View(mean.infected.per.day1.i63)


#### Next Scenario Cow 40 out of smallest community immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 40
simu1.i40 <- list()
for (i in 1:500){
  simu1.i40[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=1, immune = 40) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week1.i40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1.i40[[i]])
  for (j in 1: N_weeks){
    cows.per.week1.i40[i,j] <- length(simu1.i40[[i]][[j]][[1]])
    cows.per.week1.i40[is.na(cows.per.week1.i40)] <- 0
    print(cows.per.week1.i40[i,j])
  }
}

View(cows.per.week1.i40)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day1.i40 <- as.data.frame(colMeans(cows.per.week1.i40))

View(mean.infected.per.day1.i40)

#### End of Simulation on Cow 1 #######


#### Start Simulation Cow 38 #####

#COW 38 - Central to community 6

### 1 scenario cow 31 immune ####
#Start infection Cow 38 -> Simulating 500 times disease with immunity for cow 31
simu38.i31 <- list()
for (i in 1:500){
  simu38.i31[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=38, immune = 31) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week38.i31 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu38.i31[[i]])
  for (j in 1: N_weeks){
    cows.per.week38.i31[i,j] <- length(simu38.i31[[i]][[j]][[1]])
    cows.per.week38.i31[is.na(cows.per.week38.i31)] <- 0
    print(cows.per.week38.i31[i,j])
  }
}

View(cows.per.week38.i31)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day38.i31 <- as.data.frame(colMeans(cows.per.week38.i31))

View(mean.infected.per.day38.i31)

### Next scenario cow 9 immune #####
#Start infection Cow 38 -> Simulating 500 times disease with immunity for cow 9
simu38.i9 <- list()
for (i in 1:500){
  simu38.i9[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=38, immune = 9) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week38.i9 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu38.i9[[i]])
  for (j in 1: N_weeks){
    cows.per.week38.i9[i,j] <- length(simu38.i9[[i]][[j]][[1]])
    cows.per.week38.i9[is.na(cows.per.week38.i9)] <- 0
    print(cows.per.week38.i9[i,j])
  }
}

View(cows.per.week38.i9)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day38.i9 <- as.data.frame(colMeans(cows.per.week38.i9))

View(mean.infected.per.day38.i9)

#### Next Scenario Cow 63 out of biggest community immune #####
#Start infection Cow 38 -> Simulating 500 times disease with immunity for cow 63
simu38.i63 <- list()
for (i in 1:500){
  simu38.i63[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=38, immune = 63) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week38.i63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu38.i63[[i]])
  for (j in 1: N_weeks){
    cows.per.week38.i63[i,j] <- length(simu38.i63[[i]][[j]][[1]])
    cows.per.week38.i63[is.na(cows.per.week38.i63)] <- 0
    print(cows.per.week38.i63[i,j])
  }
}

View(cows.per.week38.i63)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day38.i63 <- as.data.frame(colMeans(cows.per.week38.i63))

View(mean.infected.per.day38.i63)


#### Next Scenario Cow 40 out of smallest community immune #####
#Start infection Cow 1 -> Simulating 500 times disease with immunity for cow 40
simu38.i40 <- list()
for (i in 1:500){
  simu38.i40[[i]] <- sim.network.transmission.immunity(adj_matrix, first.infected=38, immune = 40) #need to specify which cow gets immune by immune = 1
}

#Adding NA ifstatement
cows.per.week38.i40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu38.i40[[i]])
  for (j in 1: N_weeks){
    cows.per.week38.i40[i,j] <- length(simu38.i40[[i]][[j]][[1]])
    cows.per.week38.i40[is.na(cows.per.week38.i40)] <- 0
    print(cows.per.week38.i40[i,j])
  }
}

View(cows.per.week38.i40)

#at a prob.infect = 0.03 in 216 cases it dies out compared to 52 previously without immunity

#calculating the average
mean.infected.per.day38.i40 <- as.data.frame(colMeans(cows.per.week38.i40))

View(mean.infected.per.day38.i40)

#### End of Simulation on Cow 38 #######

### END OF SIMULATIONS #####

##### Plotting #######

#Combining all the cows mean data together
days <- as.data.frame(1:50)
mean.infected.per.day.immunity.total.left <- cbind(days, mean.infected.per.day27.i31, mean.infected.per.day27.i9, mean.infected.per.day27.i63, mean.infected.per.day27.i40,
                                                   mean.infected.per.day40.i31, mean.infected.per.day40.i9, mean.infected.per.day40.i63, mean.infected.per.day40.i1, 
                                                   mean.infected.per.day63.i31, mean.infected.per.day63.i9, mean.infected.per.day63.i15, mean.infected.per.day63.i40,
                                                   mean.infected.per.day15.i31, mean.infected.per.day15.i9, mean.infected.per.day15.i63, mean.infected.per.day15.i40,
                                                   mean.infected.per.day1.i31, mean.infected.per.day1.i9, mean.infected.per.day1.i63, mean.infected.per.day1.i40,
                                                   mean.infected.per.day38.i31, mean.infected.per.day38.i9, mean.infected.per.day38.i63, mean.infected.per.day38.i40)

View(mean.infected.per.day.immunity.total.left)

#Plotting & Analyzing the Means of all the Simulations on different Cows with Immunity



#plotting different immunity scenarios vs. no immunity on individual cows

#Cow 27
plot(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week27.i31)`, type = "l", col = 2, lty = 1, ylim = c(0, 10), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 27", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week27.i9)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week27.i63)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week27.i40)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week27)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 27 First Infected with Immunity Cow 31", "Cow 27 First Infected with Immunity Cow 9", "Cow 27 First Infected with Immunity Cow 63", 
                "Cow 27 First Infected with Immunity Cow 40", "Cow 27 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 40
plot(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week40.i31)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 40", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week40.i9)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week40.i63)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week40.i1)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week40)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 40 First Infected with Immunity Cow 31", "Cow 40 First Infected with Immunity Cow 9", "Cow 40 First Infected with Immunity Cow 63", 
                "Cow 40 First Infected with Immunity Cow 1", "Cow 40 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 63
plot(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week63.i31)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 63", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week63.i9)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week63.i15)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week63.i40)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week63)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 63 First Infected with Immunity Cow 31", "Cow 63 First Infected with Immunity Cow 9", "Cow 63 First Infected with Immunity Cow 15", 
                "Cow 63 First Infected with Immunity Cow 40", "Cow 63 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 15
plot(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week15.i31)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day
     with different immunity scenarios vs. without Immunity - Cow 15", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week15.i9)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week15.i63)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week15.i40)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week15)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 15 First Infected with Immunity Cow 31", "Cow 15 First Infected with Immunity Cow 9", "Cow 15 First Infected with Immunity Cow 63", 
                "Cow 15 First Infected with Immunity Cow 40", "Cow 15 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 1
plot(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week1.i31)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 1", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week1.i9)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week1.i63)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week1.i40)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week1)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 1 First Infected with Immunity Cow 31", "Cow 1 First Infected with Immunity Cow 9", "Cow 1 First Infected with Immunity Cow 63", 
                "Cow 1 First Infected with Immunity Cow 40", "Cow 1 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)

#Cow 38
plot(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week38.i31)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day
     with different Immunity scenarios vs. without Immunity - Cow 38", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week38.i9)`, type = "l", col = 3)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week38.i63)`, type = "l", col = 4)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week38.i40)`, type = "l", col = 5)
lines(mean.infected.per.day.immunity.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week38)`, type = "l", col = 6)
legend("topleft", 
       legend=c("Cow 38 First Infected with Immunity Cow 31", "Cow 38 First Infected with Immunity Cow 9", "Cow 38 First Infected with Immunity Cow 63", 
                "Cow 38 First Infected with Immunity Cow 40", "Cow 38 First Infected without any Immunity"), 
       col = (2:6),
       lwd = 2)


#OLD PLOT TYPE#
#Cow 79
plot(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total$`colMeans(cows.per.week79.i)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day
     with vs. without Immunity Cow 79", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.total$`colMeans(cows.per.week79)`, type = "l", col = 3)
legend("topleft", 
       legend=c("Cow 79 First Infected with Immunity", "Cow 79 First Infected without Immunity"), 
       col = (2:3),
       lwd = 2)


#plotting ALL infected averages - perhaps not necessary
#plot(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week1.i)`, type = "l", col = 2, lty = 1, ylim = c(0, 17), main = "Average Number of Infected Cows Per Day
#with Immunity", xlab = "Days", ylab = "Average Number of Infected Cows")
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week7.i)`, type = "l", col = 3)
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total.left$`colMeans(cows.per.week17.i)`, type = "l", col = 4)
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total$`colMeans(cows.per.week38.i)`, type = "l", col = 5)
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total$`colMeans(cows.per.week40.i)`, type = "l", col = 6)
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total$`colMeans(cows.per.week43.i)`, type = "l", col = 7)
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total$`colMeans(cows.per.week61.i)`, type = "l", col = 8)
#lines(mean.infected.per.day.immunity.total$`1:50`, mean.infected.per.day.immunity.total$`colMeans(cows.per.week79.i)`, type = "l", col = 9)

# Add a legend to the plot
#legend("topleft", 
#legend=c("Cow 1 First Infected", "Cow 7 First Infected", 
#"Cow 17 First Infected","Cow 38 First Infected",
#"Cow 40 First Infected","Cow 43 First Infected",
#"Cow 61 First Infected","Cow 79 First Infected"), 
# col = (2:9),
#lwd = 2)

