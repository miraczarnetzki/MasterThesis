library(igraph)

#loading our data of Network
load("/Users/miraczarnetzki/OneDrive - Högskolan Dalarna/Thesis/Disease Simulation/Data 11 Dec/network_1.25m_30min_left.Rda")

summary(igraph_df_30)

#making adjacency matrix
#from_cow <- igraph_df_30[, c('from_name')] - old network code
from_cow <- igraph_df_30[, c('from')]
#to_cow <- igraph_df_30[, c('to_name')]- old network code
to_cow <- igraph_df_30[, c('to')]
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


##########Running the disease simulation on different central cows in each of the communities #############

#COW 27 - Central to community 1

#Simulating 500 times disease without immunity for cow 27
simu27 <- list()
for (i in 1:500){
  simu27[[i]] <- sim.network.transmission(adj_matrix, first.infected=27)
}

#Adding NA ifstatement 
cows.per.week27 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu27[[i]])
  for (j in 1: N_weeks){
    cows.per.week27[i,j] <- length(simu27[[i]][[j]][[1]])
    cows.per.week27[is.na(cows.per.week27)] <- 0
    print(cows.per.week27[i,j])
  }
}

View(cows.per.week27)

#at a prob.infect = 0.03 in 303 cases it dies out

#calculating the average
mean.infected.per.day27 <- as.data.frame(colMeans(cows.per.week27))


#COW 40 - Central to community 2

#Simulating 500 times disease without immunity for cow 40
simu40 <- list()
for (i in 1:500){
  simu40[[i]] <- sim.network.transmission(adj_matrix, first.infected=40)
}

#Adding NA ifstatement 
cows.per.week40 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu40[[i]])
  for (j in 1: N_weeks){
    cows.per.week40[i,j] <- length(simu40[[i]][[j]][[1]])
    cows.per.week40[is.na(cows.per.week40)] <- 0
    print(cows.per.week40[i,j])
  }
}
View(cows.per.week40)

#at a prob.infect = 0.03 in 267 cases it dies out

#calculating the average
mean.infected.per.day40 <- as.data.frame(colMeans(cows.per.week40))



#COW 63 - Central to community 3

#Simulating 500 times disease without immunity for cow 63
simu63 <- list()
for (i in 1:500){
  simu63[[i]] <- sim.network.transmission(adj_matrix, first.infected=63)
}

#Adding NA ifstatement 
cows.per.week63 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu63[[i]])
  for (j in 1: N_weeks){
    cows.per.week63[i,j] <- length(simu63[[i]][[j]][[1]])
    cows.per.week63[is.na(cows.per.week63)] <- 0
    print(cows.per.week63[i,j])
  }
}

View(cows.per.week63)

#at a prob.infect = 0.03 in 284 cases it dies out

#calculating the average
mean.infected.per.day63 <- as.data.frame(colMeans(cows.per.week63))

View(mean.infected.per.day63)


#COW 15 - Central to community 4

#Simulating 500 times disease without immunity for cow 15
simu15 <- list()
for (i in 1:500){
  simu15[[i]] <- sim.network.transmission(adj_matrix, first.infected=15)
}

#Adding NA ifstatement 
cows.per.week15 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu15[[i]])
  for (j in 1: N_weeks){
    cows.per.week15[i,j] <- length(simu15[[i]][[j]][[1]])
    cows.per.week15[is.na(cows.per.week15)] <- 0
    print(cows.per.week15[i,j])
  }
}

View(cows.per.week15)

#at a prob.infect = 0.03 in 269 cases it dies out

#calculating the average
mean.infected.per.day15 <- as.data.frame(colMeans(cows.per.week15))

View(mean.infected.per.day15)



#COW 1 - Central to community 5

#Simulating 500 times disease without immunity for cow ID 1
simu1 <- list()
for (i in 1:500){
  simu1[[i]] <- sim.network.transmission(adj_matrix, first.infected=1)
}

#Adding NA ifstatement 

cows.per.week1 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu1[[i]])
  for (j in 1: N_weeks){
    cows.per.week1[i,j] <- length(simu1[[i]][[j]][[1]])
    cows.per.week1[is.na(cows.per.week1)] <- 0
    print(cows.per.week1[i,j])
  }
}

View(cows.per.week1)

#at a prob.infect = 0.03 in 231 cases it dies out

#calculating the average
mean.infected.per.day1 <- as.data.frame(colMeans(cows.per.week1))


#COW 38 - Central to community 6

#Simulating 500 times disease without immunity for cow 38
simu38 <- list()
for (i in 1:500){
  simu38[[i]] <- sim.network.transmission(adj_matrix, first.infected=38)
}

#Adding NA ifstatement 
cows.per.week38 <- matrix(,500,50)
for (i in 1:500){
  N_weeks <- length(simu38[[i]])
  for (j in 1: N_weeks){
    cows.per.week38[i,j] <- length(simu38[[i]][[j]][[1]])
    cows.per.week38[is.na(cows.per.week38)] <- 0
    print(cows.per.week38[i,j])
  }
}

View(cows.per.week38)

#at a prob.infect = 0.03 in 243 cases it dies out

#calculating the average
mean.infected.per.day38 <- as.data.frame(colMeans(cows.per.week38))



####### MAKING GRAPHS ###############

#Combining all the cows mean data together
days <- as.data.frame(1:50)
mean.infected.per.day.total.left <- cbind(days, mean.infected.per.day27, mean.infected.per.day40, mean.infected.per.day63, mean.infected.per.day15, mean.infected.per.day1, mean.infected.per.day38)

View(mean.infected.per.day.total.left)

#Plotting & Analyzing the Means of all the Simulations on different Cows

#plotting all infected averages
plot(mean.infected.per.day.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week27)`, type = "l", col = 2, lty = 1, ylim = c(0, 15), main = "Average Number of Infected Cows Per Day", xlab = "Days", ylab = "Average Number of Infected Cows")
lines(mean.infected.per.day.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week40)`, type = "l", col = 3)
lines(mean.infected.per.day.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week63)`, type = "l", col = 4)
lines(mean.infected.per.day.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week15)`, type = "l", col = 5)
lines(mean.infected.per.day.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week1)`, type = "l", col = 6)
lines(mean.infected.per.day.total.left$`1:50`, mean.infected.per.day.total.left$`colMeans(cows.per.week38)`, type = "l", col = 7)


# Add a legend to the plot
legend("topleft", 
       legend=c("Cow 27 First Infected", "Cow 40 First Infected", 
                "Cow 63 First Infected","Cow 15 First Infected",
                "Cow 1 First Infected","Cow 38 First Infected"),
       col = (2:9),
       lwd = 2)


