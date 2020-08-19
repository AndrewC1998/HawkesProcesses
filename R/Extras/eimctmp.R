eimctmp <- c(9246.115,  8941.131, 7180.219, 0, 0, 0, 0, 0, 0, 0,
             -3282014,         0,        0, 0, 0, 0, 0, 0, 0, 0,
             1059.976,  623.0362, 497.2070, 0, 0, 0, 0, 0, 0, 0,
            -80168.67, -29246.95,        0, 0, 0, 0, 0, 0, 0, 0, 
            -483588.4,         0,        0, 0, 0, 0, 0, 0, 0, 0,
               297.78,   123.646, 119.9172, 0, 0, 0, 0, 0, 0, 0,
            -12619.09, -4847.089,        0, 0, 0, 0, 0, 0, 0, 0,
            -44870.29, -23711.72,        0, 0, 0)         

stuff <- c(1,11,21,31,41,51,61,71)
newstuff <- cbind(mabtmp[stuff, 1:4], eimctmp[stuff])
colnames(newstuff) <- c("mu", "alpha", "beta", "intensity", "det")

library(ggplot2)

plot(newstuff$intensity, newstuff$det, type = "p", col = "red", xlab = "Branching Ratio",
     ylab = "Determinant of Empirical Fisher Information", lwd = 2)

ggplot(newstuff, aes(x=intensity, y=det, color=beta))  +
      geom_point(aes(colour = factor(beta)))


eimctmp <- c(9246.115, 8941.131, 3179.353, 3362.850, 3189.220, -342742.386, -78397.495, 0, 0, 0,
             -3282014, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1059.976, 623.0362, 0, 0, 0, 0, 0, 0, 0, 0,
             -80168.67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
             -483588.4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             297.78, 123.646, 0, 0, 0, 0, 0, 0, 0, 0,
             -12619.09, -4847.089, 0, 0, 0, 0, 0, 0, 0, 0,
             -44870.29, -23711.72, 0, 0, 0)

stuff <- c(1:7,11,21,31,41,51,61,71)
newstuff <- cbind(mabtmp[stuff, 1:4], eimctmp[stuff])
colnames(newstuff) <- c("mu", "alpha", "beta", "intensity", "det")

library(ggplot2)

plot(newstuff$intensity, newstuff$det, type = "p", col = "red", xlab = "Branching Ratio",
     ylab = "Determinant of Empirical Fisher Information", lwd = 2)

ggplot(newstuff, aes(x=intensity, y=det, color=mu))  +
  geom_point(aes(colour = factor(alpha)))

# 
# cluster run

eimctmp <- c(802.0905, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              74.18151, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             -38689.72, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              20.43103, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
             -3312.237, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             -304.5464, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             -1594.868, 0, 0, 0, 0, 0, 0, 0, 0, 0)         

mabtmp$Fisher <- eimctmp
mabtmp

eils
