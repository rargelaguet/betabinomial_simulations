library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

io <- list()
io$basedir <- "/Users/ricard/betabinomial_simulations/scalability/out"

###############################
## Load pre-computed results ##
###############################

dt.N <- fread(paste0(io$basedir, "/N.txt.gz")) %>% .[,variable:="N"] %>% setnames("N","value")
dt.D <- fread(paste0(io$basedir, "/D.txt.gz")) %>% .[,variable:="D"] %>% setnames("D","value")

to.plot <- rbind(dt.N,dt.D) %>%
  .[,time:=time/60]

##########
## Plot ##
##########

ggline(to.plot, x="value", y="time", color="inference", add = c("mean_se"), 
       palette = "jco", facet="variable") +
  # facet_wrap(~variable, scales = "fixed") +
  labs(x="", y="Time (min)") +
  theme(
    legend.title = element_blank()
  )
