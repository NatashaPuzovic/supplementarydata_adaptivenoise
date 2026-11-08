num_nodes=20

# names
node_names = paste("G", 1:num_nodes, sep = "")
meanGenotColnames <- paste("meanG_", node_names, sep = "")
meanPhenotColnames <- paste("meanP_", node_names, sep = "")
varPhenotColnames <- paste("varP_", node_names, sep = "")


enoise <- read.table("expressionnoise/cpp/tests/simEnoise.summary")
eno1 <- read.table("expressionnoise/cpp/tests/net_1/GPFsel_rep1.summary", sep = " ", header = FALSE)
colnames(eno1) <- c("generation", meanGenotColnames, meanPhenotColnames, varPhenotColnames, "meanF")

netlings <- read.table("expressionnoise/cpp/tests/simNetlings.hist")
colnames(netlings) <- c("generation", meanGenotColnames, meanPhenotColnames, varPhenotColnames, "meanF")

enoise1.4 <- read.table("expressionnoise/cpp/tests/simEnoise1.4.summary", sep = " ", header = FALSE)
colnames(enoise1.4) <- c("generation", meanGenotColnames, meanPhenotColnames, varPhenotColnames, "meanF")

enoise1.4_uneq <- read.table("expressionnoise/cpp/tests/simEnoise1.4_uneq.summary", sep = " ", header = FALSE)
colnames(enoise1.4_uneq) <- c("generation", meanGenotColnames, meanPhenotColnames, varPhenotColnames, "meanF")


firstgenNetlings <- netlings[netlings$generation == 1,]
firstgenEnoise <- eno1[eno1$generation == 1,]


firstgenNetlings[meanPhenotColnames] - firstgenEnoise[meanPhenotColnames]


firstgenNetlings[varPhenotColnames] - firstgenEnoise[varPhenotColnames]

eno1[1, meanPhenotColnames]
enoise1.4[1, meanPhenotColnames]
enoise1.4_uneq[1, meanPhenotColnames]

eno1[5, meanPhenotColnames]
enoise1.4[5, meanPhenotColnames]
enoise1.4_uneq[5, meanPhenotColnames]

eno1[10, meanPhenotColnames]
enoise1.4[10, meanPhenotColnames]
enoise1.4_uneq[10, meanPhenotColnames]

###########

noises <- read.table("supplementarydata_expressionnoise/results/20221123_ER_40node_0.05dens/net_1/GPFsel_rep1.summary")
colnames(noises) <- c("generation", meanGenotColnames)

ggplot(data = melt(noises[,c(1:11)], id.vars = "generation"), 
       aes(x=generation, y = value, group = variable)) +
  geom_line(aes(colour = variable)) +
  ylim(-1, 120)

ggplot(data = melt(noises[,c(1, 11:21)], id.vars = "generation"), 
       aes(x=generation, y = value, group = variable)) +
  geom_line(aes(colour = variable)) +
  ylim(-1, 120)















library(reshape2)
sbasalsColor = "#050075"
################# regular run

# noise genotypes
noises <- read.table("expressionnoise/cpp/netlings2.0/sim.noises")
colnames(noises) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(noises[,c(1:11)], id.vars = "generation"), 
       aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/dynamics10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[,c(1, 12:21)], id.vars = "generation"), 
       aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/dynamics10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[102:502,c(1:11)], id.vars = "generation"), 
       aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 2k-10k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/equilibrium_10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(noises[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 2k-10k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/equilibrium_10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

hist <- read.table("expressionnoise/cpp/netlings2.0/sim.hist")
colnames(hist) <- c("generation", meanPhenotColnames, varPhenotColnames, "meanF")
a <- ggplot(data = hist, 
            aes(x = generation, y = meanF)) +
  geom_line(color = fitnessColor, size = 0.5) +
  ylab("Mean fitness") +
  ylim(-0.1, 1.1)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/dynamics10kgen1kpopsize_fitness.png"),
       plot = a,
       device = "png", scale = 1.5, width = 5, height = 5, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
mean(hist$meanF[101:501])

# sbasal genotypes
sbasals <- read.table("expressionnoise/cpp/netlings2.0/sim.sbasals")
colnames(sbasals) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(sbasals[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/dynamics10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/dynamics10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 2k-10k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/equilibrium_10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(sbasals[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 2k-10k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/equilibrium_10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

################# long run

# noise genotypes
noises <- read.table("expressionnoise/cpp/netlings2.0/simLongRun.noises")
colnames(noises) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(noises[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_dynamics10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_dynamics10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 20k-100k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_equilibrium_10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(noises[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 20k-100k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_equilibrium_10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

hist <- read.table("expressionnoise/cpp/netlings2.0/simLongRun.hist")
colnames(hist) <- c("generation", meanPhenotColnames, varPhenotColnames, "meanF")
a <- ggplot(data = hist, 
            aes(x = generation, y = meanF)) +
  geom_line(color = fitnessColor, size = 0.5) +
  ylab("Mean fitness") +
  ylim(-0.1, 1.1)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_dynamics10kgen1kpopsize_fitness.png"),
       plot = a,
       device = "png", scale = 1.5, width = 5, height = 5, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
mean(hist$meanF[101:501])

# sbasal genotypes
sbasals <- read.table("expressionnoise/cpp/netlings2.0/simLongRun.sbasals")
colnames(sbasals) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(sbasals[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_dynamics10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_dynamics10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 20k-100k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_equilibrium_10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(sbasals[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 20k-100k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/longRun_equilibrium_10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

################# big pop

# noise genotypes
noises <- read.table("expressionnoise/cpp/netlings2.0/simBigPop.noises")
colnames(noises) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(noises[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_dynamics10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_dynamics10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 2k-10k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_equilibrium_10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(noises[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 2k-10k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_equilibrium_10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

hist <- read.table("expressionnoise/cpp/netlings2.0/simBigPop.hist")
colnames(hist) <- c("generation", meanPhenotColnames, varPhenotColnames, "meanF")
a <- ggplot(data = hist, 
            aes(x = generation, y = meanF)) +
  geom_line(color = fitnessColor, size = 0.5) +
  ylab("Mean fitness") +
  ylim(-0.1, 1.1)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_dynamics10kgen1kpopsize_fitness.png"),
       plot = a,
       device = "png", scale = 1.5, width = 5, height = 5, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
mean(hist$meanF[101:501])

# sbasal genotypes
sbasals <- read.table("expressionnoise/cpp/netlings2.0/simBigPop.sbasals")
colnames(sbasals) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(sbasals[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_dynamics10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_dynamics10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 2k-10k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_equilibrium_10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(sbasals[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 2k-10k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/bigPop_equilibrium_10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

################# away from peak

# noise genotypes
noises <- read.table("expressionnoise/cpp/netlings2.0/simPeak.noises")
colnames(noises) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(noises[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_dynamics10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean noise genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_dynamics10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(noises[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 2k-10k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_equilibrium_10kgen1kpopsize_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(noises[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = genotypeColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = genotypeColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = genotypeColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean noise genotype") + xlab("Mean noise genotype gen. 2k-10k") +
  ylim(-1, 10)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_equilibrium_10kgen1kpopsize_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

hist <- read.table("expressionnoise/cpp/netlings2.0/simPeak.hist")
colnames(hist) <- c("generation", meanPhenotColnames, varPhenotColnames, "meanF")
a <- ggplot(data = hist, 
            aes(x = generation, y = meanF)) +
  geom_line(color = fitnessColor, size = 0.5) +
  ylab("Mean fitness") +
  ylim(-0.1, 1.1)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_dynamics10kgen1kpopsize_fitness.png"),
       plot = a,
       device = "png", scale = 1.5, width = 5, height = 5, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
mean(hist$meanF[101:501])

# sbasal genotypes
sbasals <- read.table("expressionnoise/cpp/netlings2.0/simPeak.sbasals")
colnames(sbasals) <- c("generation", meanGenotColnames)

a <- ggplot(data = melt(sbasals[,c(1:11)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_dynamics10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[,c(1, 12:21)], id.vars = "generation"), 
            aes(x = generation, y = value, group = variable)) +
  geom_line(aes(colour = variable), size = 0.5) +
  ylab("Mean sbasal genotype") +
  ylim(-1, 120)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_dynamics10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 10, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")

a <- ggplot(data = melt(sbasals[102:502,c(1:11)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 2k-10k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_equilibrium_10kgen1kpopsize_sbasals_1-10genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")
a <- ggplot(data = melt(sbasals[102:502,c(1, 12:21)], id.vars = "generation"), 
            aes(x = variable, y = value)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_violin(colour = sbasalsColor) +
  stat_summary(fun.data = mean_sdl, geom = "pointrange", color = sbasalsColor, alpha = 0.5) + 
  stat_summary(fun = mean, colour = sbasalsColor, geom = "text", vjust = -0.7, aes( label=round(..y.., digits=1))) +
  ylab("Mean sbasal genotype") + xlab("Mean sbasal genotype gen. 2k-10k") +
  ylim(-1, 100)
ggsave(filename = sprintf("expressionnoise/expressionnoise_master/plots/misc/awayFromPeak_equilibrium_10kgen1kpopsize_sbasals_11-20genes.png"),
       plot = a,
       device = "png", scale = 1.5, width = 14, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE,
       bg = "white")























