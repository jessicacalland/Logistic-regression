

# Modelling post-vaccination strain replacement 

## LIBRARIES
library(dplyr)
library(tidyr)
library(ggplot2)


## DATA

dist_metadata <- read.csv("dist_metadata.csv", sep=";")

dist_metadata$pop <- as.factor(dist_metadata$pop)
dist_metadata$CC <- as.factor(dist_metadata$CC)

# vaccine and non-vaccine CCs
vac.cc <- c("CC206", "CC257", "CC443", "ST-4430")
non.vac.cc <- c("CC21", "CC353", "CC828", "ST-7735")


## LOGISTIC REGRESSION MODEL SCENARIOS

# We design 'survival data' according to the isolates in the vaccine trial such there is a 
# pair of observations y = survival, x = min distance to vaccine composition strains. 
# For ST-206 there would be 44 observations with y = 0, x = 'the min Hamming distance of isolate i in ST-206 to the four vaccine strains'., 
# for ST-21 25 observations with y = 1, x_i = 'the min Hamming distance of isolate i in ST-21 to the four vaccine strains', etc.

# Pre- and post vac frequencies

# Data includes four comparisons
# keep each seq only once

seq.dat <- dist_metadata %>% 
  filter(vac_seq == "6620")

freqs <- as.data.frame(table(seq.dat$CC, seq.dat$pop))
names(freqs) <- c("CC", "pop", "freq")

surv.dat <- dist_metadata %>%
  filter(vac_seq == "6620") %>%
  select(seq, CC, pop) %>%
  mutate(vac.strain = if_else(CC %in% vac.cc, 1, 0))

surv.dat2 <- left_join(surv.dat, freqs, by = c("CC", "pop"))

surv.dat3 <- surv.dat2 %>%
  mutate(surv = ifelse(pop == "post-vac", 1, 0)) %>%
  mutate(vac.dist.all = NA)

for(i in 1:nrow(surv.dat3)){
  surv.dat3$vac.dist.all[i] <- min(dist_metadata$diff[which(dist_metadata$seq == paste(surv.dat3$seq[i]))])
}


# Logistic regression curves for the vaccine trial data, x = minimum distance of each sequence to the four vaccine strains (in 1000 bases). 

mod.dat <- surv.dat3 %>%
            mutate(vac.dist.all.per1k = vac.dist.all/1000)

mod.fit <- glm(surv ~ vac.dist.all.per1k, data = mod.dat, family = "binomial")
coef(mod.fit)

theta.mod <- 0.426379
intercept.mod <- -1.481275

# Now there is too much emphasis on the large distances of the non-surviving isolates,
# let's try different scenarios of the whole population effect prediction by plugging 
# in a range of meaningful values of \theta 

# Baseline survival probability (constant term in logistic regression model) = 
# average fraction of the relative normalised frequencies of the vaccine strains present post vaccine

sum.pre <- sum(freqs$freq[which(freqs$pop == "pre-vac")])
sum.post <- sum(freqs$freq[which(freqs$pop == "post-vac")])

surv.probs <- freqs %>%
  mutate(norm.freq = case_when(pop == "pre-vac" ~ freq/sum.pre,
                          pop == "post-vac" ~ freq/sum.post))

rel.norm.freq.cc206 <- surv.probs$norm.freq[which(surv.probs$CC == "CC206" & surv.probs$pop == "post-vac")] / 
  (surv.probs$norm.freq[which(surv.probs$CC == "CC206" & surv.probs$pop == "post-vac")] + surv.probs$norm.freq[which(surv.probs$CC == "CC206" & surv.probs$pop == "pre-vac")])

rel.norm.freq.cc257 <- surv.probs$norm.freq[which(surv.probs$CC == "CC257" & surv.probs$pop == "post-vac")] / 
  (surv.probs$norm.freq[which(surv.probs$CC == "CC257" & surv.probs$pop == "post-vac")] + surv.probs$norm.freq[which(surv.probs$CC == "CC257" & surv.probs$pop == "pre-vac")])

rel.norm.freq.cc443 <- surv.probs$norm.freq[which(surv.probs$CC == "CC443" & surv.probs$pop == "post-vac")] / 
  (surv.probs$norm.freq[which(surv.probs$CC == "CC443" & surv.probs$pop == "post-vac")] + surv.probs$norm.freq[which(surv.probs$CC == "CC443" & surv.probs$pop == "pre-vac")])

rel.norm.freq.st4430 <- surv.probs$norm.freq[which(surv.probs$CC == "ST-4430" & surv.probs$pop == "post-vac")] / 
  (surv.probs$norm.freq[which(surv.probs$CC == "ST-4430" & surv.probs$pop == "post-vac")] + surv.probs$norm.freq[which(surv.probs$CC == "ST-4430" & surv.probs$pop == "pre-vac")])

bl.surv.prob <- mean(c(rel.norm.freq.cc206, rel.norm.freq.cc257, rel.norm.freq.cc443, rel.norm.freq.st4430))

intercept <- log(bl.surv.prob)
# intercept <- -2.442347
theta2 <- 0.6931472
theta2.5 <- 0.9162907
theta3 <- 1.098612

mod.dat$pred.surv.prob.mod <- 1/(1+exp(-(intercept.mod + theta.mod*(mod.dat$vac.dist.all/1000))))
mod.dat$pred.surv.prob2 <- 1/(1+exp(-(intercept + theta2*(mod.dat$vac.dist.all/1000))))
mod.dat$pred.surv.prob2.5 <- 1/(1+exp(-(intercept + theta2.5*(mod.dat$vac.dist.all/1000))))
mod.dat$pred.surv.prob3 <- 1/(1+exp(-(intercept + theta3*(mod.dat$vac.dist.all/1000))))

mod.dat$vac.strain <- factor(mod.dat$vac.strain, labels = c("non-vaccine strain", "vaccine strain"))

p <- ggplot(mod.dat) + 
  geom_jitter(aes(x=vac.dist.all/1000, y=surv, colour = vac.strain), height = 0.02, alpha = 0.8) +
  geom_line(aes(x=vac.dist.all/1000, y=pred.surv.prob.mod), colour = "grey") + 
  geom_line(aes(x=vac.dist.all/1000, y=pred.surv.prob2)) + 
  geom_line(aes(x=vac.dist.all/1000, y=pred.surv.prob2.5)) + 
  geom_line(aes(x=vac.dist.all/1000, y=pred.surv.prob3)) + 
  labs(x="Minimum distance to vaccine strains", y = "Presence post vaccine", colour = "") +
  scale_x_continuous(breaks=c(0,2,4, 6, 8, 10),
                     labels=c("0", "2k", "4k", "6k", "8k", "10k")) +
  theme_classic()
p

##  MODELLING POST-VACCINATION STRAIN REPLACEMENT


# genomes contextualised with additional large collection of chicken-associated genomes 
dist_large_col <- read.csv("large.coll.dist.csv", sep=";")

# for each sequence, choose the min vaccine distance
# --> one row per sequence

min.dist.seq <- dist_large_col %>%
  group_by(seq) %>%
  slice_min(order_by = dist) %>%
  ungroup()


# with theta corresponding to OR=2

theta2 <- 0.6931472
#intercept <- -2.442347
intercept <- log(bl.surv.prob)

# individual distance + related pred. surv. prob.
min.dist.seq$pred.surv.prob <- 1/(1+exp(-(intercept + theta2*(min.dist.seq$dist/1000))))

min.dist.seq <- min.dist.seq %>%
  group_by(CC) %>%
  mutate(ave.vac.dist.cc = mean(dist)) %>%
  mutate(min.vac.dist.cc = min(dist)) %>%
  mutate(med.vac.dist.cc = median(dist)) %>%
  mutate(med.surv.prob.cc = median(pred.surv.prob)) %>%
  mutate(vac.cc = if_else(CC %in% vac.cc, "vaccine CC", "non-vaccine CC")) %>%
  ungroup()


pre.vac.cc <- c("CC206", "CC257", "CC443", "ST-4430")
post.vac.cc <- c("CC21", "CC353", "CC828", "ST-7735")


# add variable pre/post vaccine CC and other
min.dist.seq$pre.post.other <- ifelse(min.dist.seq$CC %in% pre.vac.cc, "Pre-vaccine", 
                                      ifelse(min.dist.seq$CC %in% post.vac.cc, "Post-vaccine", "Other"))
min.dist.seq$pre.post.other <- factor(min.dist.seq$pre.post.other, levels = c("Pre-vaccine", "Post-vaccine", "Other"))



# use median vaccine distance for each CC
temp.med.cc <-  min.dist.seq %>%
  select(CC, med.vac.dist.cc) %>%
  distinct() %>%
  arrange(med.vac.dist.cc)

# arrange according to median vac distance
min.dist.seq$CC <- factor(min.dist.seq$CC, levels = temp.med.cc$CC)

# drop CCs with 3 isolates
names(which(sort(table(min.dist.seq$CC)) < 3))

pred.plot.cc <- min.dist.seq %>%
  filter(CC %in% names(which(sort(table(min.dist.seq$CC)) > 3)))


p1 <- ggplot(pred.plot.cc, aes(x = CC, y = pred.surv.prob, fill = pre.post.other)) +
  geom_boxplot() +
  #stat_summary(fun.y=mean, geom="point", shape=17, size=4, color="black", fill="black") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  labs(x="CC (ordered according to the median distance to the vaccine strains)", y = "Predicted survival probability", colour = "", fill = "") +
  scale_fill_manual(values=c("#F8766D","#00BFC4","grey")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1


# with theta corresponding to OR=2.5

theta2.5 <- 0.9162907
# intercept <- -2.442347
intercept <- log(bl.surv.prob)

# individual distance + related pred. surv. prob.
min.dist.seq$pred.surv.prob <- 1/(1+exp(-(intercept + theta2.5*(min.dist.seq$dist/1000))))

min.dist.seq <- min.dist.seq %>%
  group_by(CC) %>%
  mutate(ave.vac.dist.cc = mean(dist)) %>%
  mutate(min.vac.dist.cc = min(dist)) %>%
  mutate(med.vac.dist.cc = median(dist)) %>%
  mutate(med.surv.prob.cc = median(pred.surv.prob)) %>%
  mutate(vac.cc = if_else(CC %in% vac.cc, "vaccine CC", "non-vaccine CC")) %>%
  ungroup()

pre.vac.cc <- c("CC206", "CC257", "CC443", "ST-4430")
post.vac.cc <- c("CC21", "CC353", "CC828", "ST-7735")


# add variable pre/post vaccine CC and other
min.dist.seq$pre.post.other <- ifelse(min.dist.seq$CC %in% pre.vac.cc, "Pre-vaccine", 
                                      ifelse(min.dist.seq$CC %in% post.vac.cc, "Post-vaccine", "Other"))
min.dist.seq$pre.post.other <- factor(min.dist.seq$pre.post.other, levels = c("Pre-vaccine", "Post-vaccine", "Other"))


# use median vaccine distance for each CC
temp.med.cc <-  min.dist.seq %>%
  select(CC, med.vac.dist.cc) %>%
  distinct() %>%
  arrange(med.vac.dist.cc)

# arrange according to median vac distance
min.dist.seq$CC <- factor(min.dist.seq$CC, levels = temp.med.cc$CC)

# drop CCs with 3 isolates
names(which(sort(table(min.dist.seq$CC)) < 3))

pred.plot.cc <- min.dist.seq %>%
  filter(CC %in% names(which(sort(table(min.dist.seq$CC)) > 3)))

# boxplots
# individual distances instead of average, and then pred. surv. prob instead of average

p2 <- ggplot(pred.plot.cc, aes(x = CC, y = pred.surv.prob, fill = pre.post.other)) +
  geom_boxplot() +
  #stat_summary(fun.y=mean, geom="point", shape=17, size=4, color="black", fill="black") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  labs(x="CC (ordered according to the median distance to the vaccine strains)", y = "Predicted survival probability", colour = "", fill = "") +
  scale_fill_manual(values=c("#F8766D","#00BFC4","grey")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2


# with theta corresponding to OR=3

theta3 <- 1.098612
# intercept <- -2.442347
intercept <- log(bl.surv.prob)

# individual distance + related pred. surv. prob.
min.dist.seq$pred.surv.prob <- 1/(1+exp(-(intercept + theta3*(min.dist.seq$dist/1000))))

min.dist.seq <- min.dist.seq %>%
  group_by(CC) %>%
  mutate(ave.vac.dist.cc = mean(dist)) %>%
  mutate(min.vac.dist.cc = min(dist)) %>%
  mutate(med.vac.dist.cc = median(dist)) %>%
  mutate(med.surv.prob.cc = median(pred.surv.prob)) %>%
  mutate(vac.cc = if_else(CC %in% vac.cc, "vaccine CC", "non-vaccine CC")) %>%
  ungroup()


pre.vac.cc <- c("CC206", "CC257", "CC443", "ST-4430")
post.vac.cc <- c("CC21", "CC353", "CC828", "ST-7735")


# add variable pre/post vaccine CC and other
min.dist.seq$pre.post.other <- ifelse(min.dist.seq$CC %in% pre.vac.cc, "Pre-vaccine", 
                                      ifelse(min.dist.seq$CC %in% post.vac.cc, "Post-vaccine", "Other"))
min.dist.seq$pre.post.other <- factor(min.dist.seq$pre.post.other, levels = c("Pre-vaccine", "Post-vaccine", "Other"))


# use median vaccine distance for each CC
temp.med.cc <-  min.dist.seq %>%
  select(CC, med.vac.dist.cc) %>%
  distinct() %>%
  arrange(med.vac.dist.cc)

# arrange according to median vac distance
min.dist.seq$CC <- factor(min.dist.seq$CC, levels = temp.med.cc$CC)

# drop CCs with 3 isolates
names(which(sort(table(min.dist.seq$CC)) < 3))

pred.plot.cc <- min.dist.seq %>%
  filter(CC %in% names(which(sort(table(min.dist.seq$CC)) > 3)))

# boxplots
# individual distances instead of average, and then pred. surv. prob instead of average

p3 <- ggplot(pred.plot.cc, aes(x = CC, y = pred.surv.prob, fill = pre.post.other)) +
  geom_boxplot() +
  #stat_summary(fun.y=mean, geom="point", shape=17, size=4, color="black", fill="black") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  labs(x="CC (ordered according to the median distance to the vaccine strains)", y = "Predicted survival probability", colour = "", fill = "") +
  scale_fill_manual(values=c("#F8766D","#00BFC4","grey")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3



