library(tidyverse)
library(ggplot2)

wide = read.csv("/Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/NetworkAnalysis/output/matrix_communities.csv")
## Create Similarity matrix of Euclidean distance between points
S <- as.matrix(dist(wide[,-1]))

## Create Degree matrix
D <- matrix(0, nrow=nrow(wide[,-1]), ncol = nrow(wide[,-1])) # empty nxn matrix

for (i in 1:nrow(D)) {
  # Find top 10 nearest neighbors using Euclidean distance
  index <- order(S[i,])[2:11]
  
  # Assign value to neighbors
  D[i,][index] <- 1 
}

# find mutual neighbors
D = D + t(D) 
D[ D == 2 ] = 1

# find degrees of vertices
degrees = colSums(D) 
n = nrow(D)


## Compute Laplacian matrix
# Since k > 2 clusters, we normalize the Laplacian matrix:
laplacian = ( diag(n) - diag(degrees^(-1/2)) %*% D %*% diag(degrees^(-1/2)) )


## Compute eigenvectors
eigenvectors = eigen(laplacian, symmetric = TRUE)
n = nrow(laplacian)
eigenvectors = eigenvectors$vectors[,(n - 2):(n - 1)]

set.seed(1748)
## Run Kmeans on eigenvectors
k=5
sc = kmeans(eigenvectors, k)


## Pull clustering results
sc_results = cbind(wide[,1], cluster = as.factor(sc$cluster))
#View((sc_results))
write.csv(sc_results, file="/Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/NetworkAnalysis/output/sc_results.csv", 
          quote=F, row.names = F)
for (i in c(1:k)){
  print(sum(sc$cluster==i))
}
