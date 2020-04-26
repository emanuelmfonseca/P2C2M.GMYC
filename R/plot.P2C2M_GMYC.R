plot.P2C2M_GMYC <- function(P2C2M_GMYC.results){

x_max <-max(P2C2M_GMYC.plot.info[,1]) + max(P2C2M_GMYC.plot.info[,1])*0.1
max <- length(P2C2M_GMYC.plot.info[,1])

ggplot(data.frame(P2C2M_GMYC.plot.info[2:max,]), aes(x=P2C2M_GMYC.plot.info[2:max,])) + 
  geom_histogram(binwidth=.5, color = "black", fill = "gold3")+
  coord_cartesian(xlim = c(0, x_max)) +
  geom_vline(xintercept=P2C2M_GMYC.plot.info[1,1], color ="gray30", size=2) +
  geom_hline(yintercept= 0, color = "black", size=.5) +
  labs(x = "\nStandard deviation of the simulated number of species",
       y = "Frequency\n") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

}
