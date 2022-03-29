



setwd('F:/Corbi/pagerank_paper+引文在ppt的备注/code')
## ------------------------------------------------------------------------
rm(list=ls(all=TRUE))
library(Corbi)
library(Matrix)
options(scipen=0)
source("markrank.R")

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

net <- read.csv("hESCs_to_neurons_matrix_tmp_network.csv", header=F)

dataset <- read.csv("hESCs_to_neurons_matrix.csv", header=T)

dataset<-t(dataset)

label_num<-c(40,504,278,595,502,765)



count<-5
res<-c(1:5)
select_num<-0.05*ncol(dataset)
ln<-1


all<-ncol(dataset)
result_matrix<-matrix(1:all,nrow=count,ncol=ncol(dataset))


for (j in 1:count){
  num<-sum(label_num[j]+label_num[j+1])
  data<-dataset[(1:num),]
  data[1:num,]<-dataset[ln:sum(label_num[1:(j+1)]),]
  ln<-1+sum(label_num[1:j])
  
  sample_num <- num
  normal_num <- label_num[j]
  label <- c(rep(0,  normal_num), rep(1, sample_num- normal_num))
  
  for (n in 1:ncol(data)) {
    x<-data[,n]
    if (sd(x)==0) {
      data[,n]<-runif(nrow(data), min = 0, max = 1)
    }
    
  }
  data<-scale(data)
  colnames(data)<-c(1:ncol(dataset))
  
  time1 <- system.time(
    result <- markrank_revise(data, net, label, normal_num, alpha=0.85)
  )
  
  s1 <- sort(result$score, decreasing=TRUE)
  print(s1)
  res[j]<-mean(s1[1:select_num])
  result_matrix[j,1:ncol(dataset)]<-result$score
}

write.csv(as.matrix(result_matrix),file = "hESCs_to_neurons_DNRS_matrix.csv")
x<-1:5
plot(x,res) 
lines(x,res,type = 'l',col = 'blue',lwd = 2)







