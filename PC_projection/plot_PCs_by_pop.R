# Define input arguments
args <- commandArgs(TRUE)

pcfile <- args[1] # name of PC file
numpcs <- args[2] # number of PCs

# load ggplot
library(ggplot2)

# read in PCs
dfpc <- read.table(gzfile(pcfile),header=T)

# loop over PC pairs
pcs=seq(1,as.numeric(numpcs)-1)
for(pc in pcs) {
   dfpcp <- data.frame("POPN"=as.factor(dfpc$POP),x=dfpc[[paste0("PC",pc,"_AVG")]],y=dfpc[[paste0("PC",pc+1,"_AVG")]])
   plot_ <- ggplot(dfpcp,aes(x=x,y=y,colour=POPN)) + geom_point() + labs(title=paste0("1000G PCs - PC",pc+1," vs. PC",pc),x=paste0("PC",pc),y=paste0("PC",pc+1)) + theme_bw()
   ggsave(plot=plot_,device="png",filename=paste0("pc_plot_PC",pc,"_and_PC",pc+1,".png"))
}
