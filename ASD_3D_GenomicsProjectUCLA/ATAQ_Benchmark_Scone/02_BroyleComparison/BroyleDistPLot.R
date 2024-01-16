library(ggplot2)

#load the distribution from the data folder to plot.
#each distribution was yielded comparing all the samples on the Broyle data set to one of the peak sets on the data folder.


DistEncode=readRDS("DistEncdoePipevsBroyle.RDS")
DistJennyNG=readRDS("DistJennyPipevsBroyleNG.RDS")
df=data.frame(cbind(DistEncodeNG,DistJennyNG))
df=reshape2::melt(df)


ggplot(df,aes(x=value,fill=variable)) + geom_histogram(position = "identity",alpha=0.5,bins=100) +  geom_vline(aes(xintercept=mean(DistEncode)), color="red",linetype="dashed") + geom_vline(aes(xintercept=mean(DistJennyNG)),color="blue",linetype="dashed") +labs(x="Overlap percentage with Broyle samples", y="Frequency", fill="Pipeline")



