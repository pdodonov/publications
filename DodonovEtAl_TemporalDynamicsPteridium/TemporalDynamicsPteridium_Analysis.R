library(nlme)
library(bbmle)
library(lme4)

patch <- rep(c("1a","1b","2a", "2b", "3", "4", "5"),4)
year <- rep(c(2009, 2011,2013,2015), each=7)
size <- c(20069, NA, 6825, NA, 1758, 470, 95,  20360, 1121, 6550, 849, 0, 802, 235, 19174, 1146, 7102, 762, 12, 937, 200, 19616, 1178, 7424, 840, 0, 802, 224)

dados <- data.frame(patch,year,size)
dados$patch<-as.factor(dados$patch)

dados1 <- subset(dados, year>2009)
dados2 <- subset(dados, !is.na(size))

dados1b <- subset(dados1, patch!="3")
dados2b <- subset(dados2, patch!="3")

lme.year1b <- lme(size ~ year, random = ~1|patch, data=dados1b)
qqnorm(resid(lme.year1b));qqline(resid(lme.year1b))

# Non-normal

# Log
lme.year1b.log <- lme(log(size) ~ year, random = ~1|patch, data=dados1b)
qqnorm(resid(lme.year1b.log));qqline(resid(lme.year1b.log))

# Seems normal
summary(lme.year1b.log)
lme.year1b.log.null <- lme(log(size) ~ 1, random = ~1|patch, data=dados1b)
AICctab(lme.year1b.log, lme.year1b.log.null)
lme.year2b <- lme(size ~ year, random = ~1|patch, data=dados2b)

qqnorm(resid(lme.year2b)); qqline(resid(lme.year2b))
# Non-normal

# Log
lme.year2b.log <- lme(log(size) ~ year, random = ~1|patch, data=dados2b)
qqnorm(resid(lme.year2b.log)); qqline(resid(lme.year2b.log))
lme.year2b.log.null <- lme(log(size) ~ 1, random = ~1|patch, data=dados2b)

AICctab(lme.year2b.log, lme.year2b.log.null)

lme.year2b.null <- lme(size ~ 1, random = ~1|patch, data=dados2b)

# How about lme4?
lmer.year1b.log <- lmer(log(size) ~ year + (1|patch), data=dados1b)
lmer.year1b.log.null <- lmer(log(size) ~ (1|patch), data=dados1b)
AICctab(lmer.year1b.log, lmer.year1b.log.null)

lmer.year2b.log <- lmer(log(size) ~ year + (1|patch), data=dados2b)
lmer.year2b.log.null <- lmer(log(size) ~ (1|patch), data=dados2b)
AICctab(lmer.year2b.log, lmer.year2b.log.null)

lmer.year1.log <- lmer(log(size+1) ~ year + (1|patch), data=dados1)
lmer.year1.log.null <- lmer(log(size+1) ~ (1|patch), data=dados1)
qqnorm(resid(lmer.year1.log)); qqline(resid(lmer.year1.log))
AICctab(lmer.year1.log, lmer.year1.log.null)

lmer.year2.log <- lmer(log(size+1) ~ year + (1|patch), data=dados2)
lmer.year2.log.null <- lmer(log(size+1) ~ (1|patch), data=dados2)
qqnorm(resid(lmer.year2.log)); qqline(resid(lmer.year2.log))
AICctab(lmer.year2.log, lmer.year2.log.null)

plot(log(size) ~ year, col=as.numeric(patch), data=dados1b, pch=16)
plot(log(size) ~ year, col=as.numeric(patch), data=dados2b, pch=16)

par(mfrow=c(2,2))
qqnorm(resid(lme.year1)); qqline(resid(lme.year1))
qqnorm(resid(lme.year2)); qqline(resid(lme.year2))
qqnorm(resid(lme.year1b)); qqline(resid(lme.year1b))
qqnorm(resid(lme.year2b)); qqline(resid(lme.year2b))




