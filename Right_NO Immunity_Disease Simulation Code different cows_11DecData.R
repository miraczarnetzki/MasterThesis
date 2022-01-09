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
#column_names <- data.frame(colnames(adj_matrix)) - not needed for new network code


#colnames(adj_matrix) <- c(1:dim(adj_matrix)[1])
#rownames(adj_matrix) <- c(1:dim(adj_matrix)[1])
#adj_matrix
class(adj_matrix)
dim(adj_matrix)
  #79 cows
adj_matrix[1,1]

#simulation function without immunity
#at a probability treshhold of infection probability = 0.03 50% of the time disease dies out and 50% it spreads until 50 days
sim.network.transmission <- function(adj_matrix, first.infected = 1, prob.infected = 0.03, prob.cured = 0.1, n.weeks = 50) {
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
            if(not.infected) {
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

######Simulating on different cows ######

#COW 136 - Central to community 1

#Simulating 500 times disease without immunity for cow ID 136
simu136 <- list()
for (i in 1:500){
  simu136[[i]] <- sim.network.transmission(adj_matrix, first.infected=136)
}

#Adding NA ifstatement 

cows.per.week136 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu136[[i]])
  for (j in 1: N_weeks){
    cows.per.week136[i,j] <- length(simu136[[i]][[j]][[1]])
    cows.per.week136[is.na(cows.per.week136)] <- 0
    print(cows.per.week136[i,j])
  }
}

View(cows.per.week136)

#at a prob.infect = 0.03 in 197 cases it dies out

#calculating the average
mean.infected.per.day136 <- as.data.frame(colMeans(cows.per.week136))


#COW 158 - Central to community 2

#Simulating 500 times disease without immunity for cow 158
simu158 <- list()
for (i in 1:500){
  simu158[[i]] <- sim.network.transmission(adj_matrix, first.infected=158)
}

#Adding NA ifstatement 
cows.per.week158 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu158[[i]])
  for (j in 1: N_weeks){
    cows.per.week158[i,j] <- length(simu158[[i]][[j]][[1]])
    cows.per.week158[is.na(cows.per.week158)] <- 0
    print(cows.per.week158[i,j])
  }
}

View(cows.per.week158)

#at a prob.infect = 0.03 in 171 cases it dies out

#calculating the average
mean.infected.per.day158 <- as.data.frame(colMeans(cows.per.week158))


#COW 92 - Central to community 3

#Simulating 500 times disease without immunity for cow 92
simu92 <- list()
for (i in 1:500){
  simu92[[i]] <- sim.network.transmission(adj_matrix, first.infected=92)
}

#Adding NA ifstatement 
cows.per.week92 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu92[[i]])
  for (j in 1: N_weeks){
    cows.per.week92[i,j] <- length(simu92[[i]][[j]][[1]])
    cows.per.week92[is.na(cows.per.week92)] <- 0
    print(cows.per.week92[i,j])
  }
}
View(cows.per.week92)

#at a prob.infect = 0.03 in 215 cases it dies out

#calculating the average
mean.infected.per.day92 <- as.data.frame(colMeans(cows.per.week92))


#COW 116 - Central to community 4

#Simulating 500 times disease without immunity for cow 95
simu116 <- list()
for (i in 1:500){
  simu116[[i]] <- sim.network.transmission(adj_matrix, first.infected=116)
}

#Adding NA ifstatement 
cows.per.week116 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu116[[i]])
  for (j in 1: N_weeks){
    cows.per.week116[i,j] <- length(simu116[[i]][[j]][[1]])
    cows.per.week116[is.na(cows.per.week116)] <- 0
    print(cows.per.week116[i,j])
  }
}

View(cows.per.week116)

#at a prob.infect = 0.03 in 188 cases it dies out

#calculating the average
mean.infected.per.day116 <- as.data.frame(colMeans(cows.per.week116))


#COW 84 - Central to community 5

#Simulating 500 times disease without immunity for cow 110
simu84 <- list()
for (i in 1:500){
  simu84[[i]] <- sim.network.transmission(adj_matrix, first.infected=84)
}

#Adding NA ifstatement 
cows.per.week84 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu84[[i]])
  for (j in 1: N_weeks){
    cows.per.week84[i,j] <- length(simu84[[i]][[j]][[1]])
    cows.per.week84[is.na(cows.per.week84)] <- 0
    print(cows.per.week84[i,j])
  }
}

View(cows.per.week84)

#at a prob.infect = 0.03 in 276 cases it dies out

#calculating the average
mean.infected.per.day84 <- as.data.frame(colMeans(cows.per.week84))

View(mean.infected.per.day84)

#####End of simulations######

####Plotting####
#Combining all the cows mean data together
days <- as.data.frame(1:50)
mean.infected.per.day.total.right <- cbind(days, mean.infected.per.day136, mean.infected.per.day158, mean.infected.per.day92, mean.infected.per.day116, mean.infected.per.day84)

View(mean.infected.per.day.total.right)

#Plotting & Analyzing the Means of all the Simulations on different Cows

#plotting all infected averages
plot(mean.infected.per.day.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week136)`, type = "l", col = 2, lty = 1, ylim = c(0, 30), main = "Average Number of Infected Cows Per Day", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week158)`, type = "l", col = 3)
lines(mean.infected.per.day.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week92)`, type = "l", col = 4)
lines(mean.infected.per.day.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week116)`, type = "l", col = 5)
lines(mean.infected.per.day.total.right$`1:50`, mean.infected.per.day.total.right$`colMeans(cows.per.week84)`, type = "l", col = 6)

# Add a legend to the plot
legend("topleft", 
       legend=c("Cow 136 First Infected", "Cow 158 First Infected", 
                "Cow 92 First Infected","Cow 116 First Infected",
                "Cow 84 First Infected"),
       col = (2:6),
       lwd = 2)


