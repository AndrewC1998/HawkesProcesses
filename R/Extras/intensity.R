t.s <- c(0,1,2,3,4,5,6,7,8)
a = 2; b = 2.2; mu = 0.5
l = c(0)
x = c(0)
for (i in 2:length(t.s)){
  xx = seq(t.s[i-1], t.s[i], length.out = 100)[-1]
  
  x <- c(x, xx)
  
  l <- c(l, l[length(l)] * exp(-b*(xx-t.s[i-1])))
  l <- c(l, l[length(l)] + a)
  
  x <- c(x, t.s[i])
}
x = x[-length(x)]
l = l[-length(l)]
l = l + rnorm(800, 0, 0.05)

plot(x, mu+l, 'l',
     ylim = c(mu-.6, mu+max(l)),
     ylab = TeX('Intensity - $\\lambda$'), xlab = 'Time', col = "blue")
grid(20,20)

t.s <- c(0,1.1,2.1,3.15,4.2,5.1,6.15,7.4,9)
a = 1.5; b = 2; mu = 1
l = c(0)
x = c(0)
for (i in 2:length(t.s)){
  xx = seq(t.s[i-1], t.s[i], length.out = 100)[-1]
  
  x <- c(x, xx)
  
  l <- c(l, l[length(l)] * exp(-b*(xx-t.s[i-1])))
  l <- c(l, l[length(l)] + a)
  
  x <- c(x, t.s[i])
}
x = x[-length(x)]
l = l[-length(l)] 
l = l + rnorm(800, 0, 0.05)

lines(x, l, col = "red")
legend("topleft", inset=.01, c("Dim 1","Dim 2"), 
       col=c("red", "blue"), lty = 1, cex=0.6)
