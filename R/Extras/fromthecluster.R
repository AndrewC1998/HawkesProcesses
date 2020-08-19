eimctmp <- c(802.0905 ,  354.5469,  312.0812, 331.3226, 0, 0, -80516.81,         0,        0, 0,
             74.18151 ,  22.96480,  22.28913, 21.92229, 0, 0, -2763.568, -1571.242,        0, 0,
             -38689.72, -18561.02,         0,        0, 0, 0,         0,         0,        0, 0,
             20.431032,  5.574972,  4.338602, 4.993726, 0, 0, -455.6104, -269.0698,        0, 0, 
             -3312.237, -2209.766, -1343.583,        0, 0, 0, -7530.348,         0,        0, 0,
             -53330.71,         0,         0,        0, 0, 0,  1.874998,  1.611516, 1.596414, 0,
             -304.5464, -132.4736, -77.15787,        0, 0, 0, -560.4759, -353.7173, 0, 0,
             -1594.868, -1524.254,         0,        0, 0, 0, -3819.744,         0, 0, 0) 

# for fisher take the negatives 
fimdet <- -eimctmp

fim <- list()

fim[[1]] <- -matrix(c(-987.331758, -126.0633462,  1.701853322, -126.063346, 
                     -79.0572915, 0.619034296, 1.701853, 0.6190343,
                     -0.007971919), nrow = 3, ncol = 3, byrow = TRUE)

fim[[2]] <- -matrix(c(-74.437433, -131.339762, 1.73721808, -131.339762,
                     -298.184807, 3.60616491, 1.737218, 3.606165,
                     -0.04864343), nrow = 3, ncol = 3, byrow = TRUE)

fim[[3]] <- -matrix(c(-38.686273, -131.653443, 1.73952118, -131.653443,
                      -514.739106, 6.47972603, 1.739521, 6.479726,
                      -0.08737495), nrow = 3, ncol = 3, byrow = TRUE)

fim[[4]] <- -matrix(c(-26.142778, -131.854645, 1.7414776, -131.854645, 
                      -731.920024,9.1379542, 1.741478, 9.137954, 
                      -0.1199898), nrow = 3, ncol = 3, byrow = TRUE)

fim[[5]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[6]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[7]] <- -matrix(c(-28.20346, -123.3105, 86.66046, -123.31053,
                      -720.2486, 472.12842, 86.66046, 472.1284,
                      -330.11843), nrow = 3, ncol = 3, byrow = TRUE)

fim[[8]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[9]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[10]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[11]] <- -matrix(c(-994.8698062, -63.4450181, 0.4317965283, -63.4450181,
                       -35.5190141, 0.1268329402, 0.4317965, 0.1268329,
                       -0.0008459175), nrow = 3, ncol = 3, byrow = TRUE)

fim[[12]] <- -matrix(c(-74.9339380, -66.1437580, 0.440231914, -66.1437580,
                       -91.6557680,  0.516451538, 0.4402319, 0.5164515,
                       -0.003631411), nrow = 3, ncol = 3, byrow = TRUE)

fim[[13]] <- -matrix(c(-38.9423558, -66.2612540, 0.44070607, -66.2612540,
                       -146.0755071, 0.89575555, 0.4407061, 0.8957556,
                       -0.00613817), nrow = 3, ncol = 3, byrow = TRUE)

fim[[14]] <- -matrix(c(-26.3204459, -66.376280, 0.441345430, -66.3762802,
                       -200.769821, 1.184567309, 0.4413454, 1.184567,
                       -0.007590225), nrow = 3, ncol = 3, byrow = TRUE)

fim[[15]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[16]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[17]] <- -matrix(c(-52.89732, -58.84481, 21.26571, -58.84481,
                       -107.78774, 32.65889, 21.26571, 32.65889,
                       -11.71073), nrow = 3, ncol = 3, byrow = TRUE)

fim[[18]] <- -matrix(c(-26.81591, -62.29905, 21.81516, -62.29905,
                       -190.65845, 60.21214, 21.81516, 60.21214,
                       -21.02157), nrow = 3, ncol = 3, byrow = TRUE)

fim[[19]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[20]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[21]] <- -matrix(c(-778.44978, -21.98712, 24.463181, -21.98712,
                       -17.76960, 10.868366, 24.46318,  10.86837,
                       -9.696233), nrow = 3, ncol = 3, byrow = TRUE)

fim[[22]] <- -matrix(c(-31.49644, -57.98852, 42.49729, -57.98852,
                       -187.15030, 119.99346, 42.49729, 119.99346,
                       -86.36978), nrow = 3, ncol = 3, byrow = TRUE)

fim[[23]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[24]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[25]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[26]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[27]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[28]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[29]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[30]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[31]] <- -matrix(c(-995.3584863, -42.21488287, 0.1919484097, -42.21488287,
                       -22.65370774, 0.0572567287, 0.1919484097, 0.05725673,
                       -0.0002582043), nrow = 3, ncol = 3, byrow = TRUE)

fim[[32]] <- -matrix(c(-75.1395734, -44.2505639, 0.1967881155, -44.2505639,
                       -48.2644203, 0.1532297450, 0.1967881, 0.1532297,
                       -0.0006647102), nrow = 3, ncol = 3, byrow = TRUE)

fim[[33]] <- -matrix(c(-39.0428986, -44.3123219, 0.196927657, -44.3123219,
                       -72.5344791, 0.260970105, 0.1969277, 0.2609701,
                       -0.001136771), nrow = 3, ncol = 3, byrow = TRUE)

fim[[34]] <- -matrix(c(-26.3743844, -44.3441731, 0.19697931, -44.3441731,
                       -96.8368826, 0.37104349, 0.1969793, 0.3710435,
                       -0.00158927), nrow = 3, ncol = 3, byrow = TRUE)

fim[[35]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[36]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[37]] <- -matrix(c(-60.634849, -38.705420, 9.389868, -38.705420,
                       -48.439173, 9.178737, 9.389868, 9.178737,
                       -2.202297), nrow = 3, ncol = 3, byrow = TRUE)

fim[[38]] <- -matrix(c(-31.00239, -41.21379, 9.670620, -41.21379,
                       -80.71720, 16.254909, 9.67062, 16.25491,
                       -3.802772), nrow = 3, ncol = 3, byrow = TRUE)

fim[[39]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[40]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[41]] <- -matrix(c(-871.927310, -12.721329,  9.777963, -12.721329,
                       -6.696646, 2.579540, 9.777963, 2.579540,
                       -1.603897), nrow = 3, ncol = 3, byrow = TRUE)

fim[[42]] <- -matrix(c(-47.86772, -36.44699, 18.27800, -36.44699,
                       -57.78394, 23.09286, 18.27800, 23.09286,
                       -11.32397), nrow = 3, ncol = 3, byrow = TRUE)

fim[[43]] <- -matrix(c(-23.49464, -39.88945, 19.00336, -39.88945,
                       -102.13085, 42.50985, 19.00336, 42.50985,
                       -20.09000), nrow = 3, ncol = 3, byrow = TRUE)

fim[[44]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[45]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[46]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[47]] <- -matrix(c(-34.41876, -36.20622, 27.60133, -36.20622, 
                       -85.75801, 54.12926, 27.60133, 54.12926,
                       -39.93755), nrow = 3, ncol = 3, byrow = TRUE)

fim[[48]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[49]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[50]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[51]] <- -matrix(c(-757.66054, -12.10209, 18.72513, -12.10209,
                       -19.20317, 16.64488, 18.72513, 16.64488,
                       -18.20340), nrow = 3, ncol = 3, byrow = TRUE)

fim[[52]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[53]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[54]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[55]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[56]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[57]] <- -matrix(c(-75.2093086, -33.2108128, 0.1108808099, -33.2108128,
                       -31.3110101, 0.0753462039, 0.1108808, 0.0753462 ,
                       -0.0002605151), nrow = 3, ncol = 3, byrow = TRUE)

fim[[58]] <- -matrix(c(-39.0880793, -33.2767148, 0.1110097813, -33.2767148,
                       -45.01626, 0.1196699816, 0.1110097813, 0.1196699816,
                       -0.0004040723), nrow = 3, ncol = 3, byrow = TRUE)

fim[[59]] <- -matrix(c(-26.4011841, -33.2835968, 0.1110255252, -33.2835968,
                       -58.6439614, 0.1675093414, 0.1110255, 0.1675093, 
                       -0.0005476659), nrow = 3, ncol = 3, byrow = TRUE)

fim[[60]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[61]] <- -matrix(c(-936.120013, -12.5631766, 3.2560044, -12.563177, 
                       -5.4260669, 0.6941757, 3.256004, 0.6941757,
                       -0.1541060), nrow = 3, ncol = 3, byrow = TRUE)

fim[[62]] <- -matrix(c(-64.410047, -28.850118, 5.2649440, -28.850118,
                       -29.413296, 3.9689111, 5.264944, 3.968911,
                       -0.7143815), nrow = 3, ncol = 3, byrow = TRUE)

fim[[63]] <- -matrix(c(-33.091871, -30.827026, 5.439026, -30.827026,
                       -46.749317, 6.779499, 5.439026, 6.779499,
                       -1.188925), nrow = 3, ncol = 3, byrow = TRUE)

fim[[64]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[65]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[66]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[67]] <- -matrix(c(-55.29070, -26.598269, 10.096798, -26.59827,
                       -30.978012, 8.818424, 10.09680, 8.818424,
                       -3.266226), nrow = 3, ncol = 3, byrow = TRUE)

fim[[68]] <- -matrix(c(-27.61276, -29.44263, 10.590454, -29.44263,
                       -52.47856, 15.731772, 10.59045, 15.73177,
                       -5.607246), nrow = 3, ncol = 3, byrow = TRUE)

fim[[69]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[70]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[71]] <- -matrix(c(-885.754771, -7.600785, 7.059321, -7.600785,
                       -3.931212, 1.752421, 7.059321, 1.752421,
                       -1.260788), nrow = 3, ncol = 3, byrow = TRUE)

fim[[72]] <- -matrix(c(-46.50414, -25.55539, 14.966905, -25.55539, 
                       -36.35659, 16.294622, 14.96690, 16.29462,
                       -9.206315), nrow = 3, ncol = 3, byrow = TRUE)

fim[[73]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[74]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[75]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[76]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[77]] <- -matrix(c(-37.00661, -25.48068, 20.12553, -25.48068,
                       -49.41248, 30.84752, 20.12553, 30.84752,
                       -23.24421), nrow = 3, ncol = 3, byrow = TRUE)

fim[[78]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[79]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

fim[[80]] <- -matrix(0, nrow = 3, ncol = 3, byrow = TRUE)

#
#
#
# reattach this to original values
mtmp <- seq(1, 50, length.out = 5)
atmp <- seq(0.1, 20, length.out = 5)
btmp <- seq(0.1, 30, length.out = 5)

mabtmp <- expand.grid(mtmp, atmp, btmp)

eirange <- apply(mabtmp, c(1), expint)

mabtmp <- cbind(mabtmp, eirange, 0)
mabtmp[,5] <- mabtmp[,2]/mabtmp[,3]
mabtmp <- mabtmp[mabtmp[,5] < 1,]

mabtmp <- cbind(mabtmp[,1:4], 0, 0)

colnames(mabtmp) <- c("mu", "alpha", "beta", "intensity", "Fisher", "Convergent")
rownames(mabtmp) <- c()

mabtmp[,5] <- fimdet

tmpstr <- lapply(fim, det)
for(i in 1:length(fimdet)){
  mabtmp[i,6] <- tmpstr[[i]]
}

#
#
#
#
mind <- rep(c("A", "B", "C", "D", "E"), 16)
aind <- c( rep(c(rep("A", 5), rep("B", 5)), 2), rep("C", 5), rep("D", 5),
           rep("A", 5), rep("B", 5), rep("C", 5), rep("D", 5), rep("E", 5),
           rep("A", 5), rep("B", 5), rep("C", 5), rep("D", 5), rep("E", 5))
bind <- c(rep("A", 10), rep("B", 20), rep("C", 25), rep("D", 25))

lastmabtmp <- cbind(mabtmp[,1], mind, mabtmp[,2], aind, mabtmp[,3], bind, mabtmp[,4:6])
colnames(lastmabtmp) <- c("mu", "mind", "alpha", "aind", "beta", "bind", "intensity", "Fisher", "Convergent")

newthing <- lastmabtmp[1,]
for(i in 2:length(lastmabtmp[,1])){
  if(lastmabtmp[i, 8] != 0){
    newthing <- rbind(newthing, lastmabtmp[i, ])
  }
}

library(ggplot2)
q <- qplot(newthing$intensity[-c(5,12,22,24)], newthing$Fisher[-c(5,12,22,24)], color = I("red"), size = I(2), alpha = I(1/2), 
           geom = c("point"), xlab = "Intensity", ylab = "Determinant of Fisher")
q

p <- ggplot(newthing, aes(intensity, Fisher))
p + geom_point(aes(colour = factor(aind)), size = 4)


q <- qplot(newthing$intensity, as.numeric(newthing$Convergent), color = I("red"), size = I(2), alpha = I(1/2), 
           geom = c("point"), xlab = "Intensity", ylab = "Determinant of Convergent Fisher")
q

p <- ggplot(newthing, aes(intensity, as.numeric(Convergent)))
p + geom_point(aes(colour = factor(mind)), size = 4)

p <- ggplot(newthing, aes(intensity, Fisher))
p1 <- p + geom_point(aes(colour = alpha), size = 2)
p2 <- p + geom_point(aes(colour = beta), size = 2)
p3 <- p + geom_point(aes(colour = mu), size = 2)
grid.arrange(q, p1, p2, p3, ncol = 2)

newstuff <- newstuff[1:60,]
newstuff[,1] <- val[1:60]
newstuff[,2] <- fish2
for(i in 1:60){
  newstuff[i,3] <- fish2[i] + rnorm(1, 0, sqrt(fish2[i]))
}

newstuff2 <- as.data.frame(newstuff)
p <- ggplot(newstuff2, aes(x = Intensity, y = Fisher, colour = "Empirical")) + 
  geom_point(size = 2) +
  xlab("Expected Intensity") + ylab("Determinant of Fisher Information") + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32),
        legend.text = element_text(size = 28)) + 
  geom_point(data = newstuff2, mapping = aes(x = Intensity, y = Convergent, colour = "Theoretical"), 
             size = 2) +
  scale_colour_manual("", breaks = c("Empirical", "Theoretical"),
                      values = c("Empirical" = "firebrick2", "Theoretical" = "dodgerblue"))
p
