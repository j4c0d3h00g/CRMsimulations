library(CRMwf)

r <- seq(-10, 10, 0.001)

# Plot of Hampel weight function
plot(r, HampelWeightFunction(r, q1 = 1.645, q2 = 1.960, q3 = 3.090), type = "l", col = "black", ylab = "w(r)")
curve(HampelWeightFunction(x, q1 = 1.487, q2 = 2.974, q3 = 5.948), type = "l", col = "red", ylab = "w(r)", add = TRUE)
legend(x = "topright", legend = c("(Q1, Q2, Q3) = \n(1.487, 2.974, 5.948)", "(Q1, Q2, Q3) = \n(1.645, 1.960, 3.090)"), 
       lty = 1, lwd = 2, cex = 1, y.intersp = 1.3,
       col = c("red", "black"))

# Plot of Tukey bisquare weight function
plot(r, TukeyWeightFunction(r, q = 4.685), type = "l", col = "black", ylab = "w(r)")
curve(TukeyWeightFunction(x, q = 3.5), type = "l", col = "red", ylab = "w(r)", add = TRUE)
curve(TukeyWeightFunction(x, q = 5.5), type = "l", col = "orange", ylab = "w(r)", add = TRUE)
curve(TukeyWeightFunction(x, q = 6.5), type = "l", col = "yellow", ylab = "w(r)", add = TRUE)
legend(x = "topright", legend = c("Q = 3.5", "Q = 4.685", "Q = 5.5", "Q = 6.5"), 
       lty = 1, lwd = 2, cex = 1, y.intersp = 1,
       col = c("red", "black", "orange", "yellow"))

# Plot of Huber weight function
plot(r, HuberWeightFunction(r, q = 1.345), type = "l", col = "black", ylab = "w(r)")
curve(HuberWeightFunction(x, q = 1), type = "l", col = "red", ylab = "w(r)", add = TRUE)
curve(HuberWeightFunction(x, q = 2), type = "l", col = "orange", ylab = "w(r)", add = TRUE)
curve(HuberWeightFunction(x, q = 3), type = "l", col = "yellow", ylab = "w(r)", add = TRUE)
legend(x = "topright", legend = c("Q = 1", "Q = 1.345", "Q = 2", "Q = 3"), 
       lty = 1, lwd = 2, cex = 1, y.intersp = 1,
       col = c("red", "black", "orange", "yellow"))

# Plot of Andrews-sine weight function
plot(r, AndrewsWeightFunction(r, q = 1.339), type = "l", col = "black", ylab = "w(r)")
curve(AndrewsWeightFunction(x, q = 0.75), type = "l", col = "red", ylab = "w(r)", add = TRUE)
curve(AndrewsWeightFunction(x, q = 1.75), type = "l", col = "orange", ylab = "w(r)", add = TRUE)
curve(AndrewsWeightFunction(x, q = 2.25), type = "l", col = "yellow", ylab = "w(r)", add = TRUE)
legend(x = "topright", legend = c("Q = 0.75", "Q = 1.339", "Q = 1.75", "Q = 2.25"), 
       lty = 1, lwd = 2, cex = 1, y.intersp = 1,
       col = c("red", "black", "orange", "yellow"))

# Plot of Generalized Gauss weight function
plot(r, GaussWeightFunction(r, q = 1.063, a = 1.387, b = 1.5), type = "l", col = "black", ylab = "w(r)")
curve(GaussWeightFunction(x, q = 0.75, a = 1.387, b = 1.5), type = "l", col = "red", ylab = "w(r)", add = TRUE)
curve(GaussWeightFunction(x, q = 1.25, a = 1.387, b = 1.5), type = "l", col = "orange", ylab = "w(r)", add = TRUE)
curve(GaussWeightFunction(x, q = 1.5, a = 1.387, b = 1.5), type = "l", col = "yellow", ylab = "w(r)", add = TRUE)
legend(x = "topright", legend = c("Q = 0.75", "Q = 1.063", "Q = 1.25", "Q = 1.5"), 
       lty = 1, lwd = 2, cex = 1, y.intersp = 1,
       col = c("red", "black", "orange", "yellow"))

# Plot of linear quadratic quadratic weight function
plot(r, QuadraticWeightFunction(r, q1 = 0.982, q2 = 1.473, s = 1.5), type = "l", col = "black", ylab = "w(r)")
curve(QuadraticWeightFunction(x, q1 = 0.5, q2 = 0.75, s = 1.5), type = "l", col = "red", ylab = "w(r)", add = TRUE)
curve(QuadraticWeightFunction(x, q1 = 1.5, q2 = 2.25, s = 1.5), type = "l", col = "orange", ylab = "w(r)", add = TRUE)
curve(QuadraticWeightFunction(x, q1 = 2, q2 = 3, s = 1.5), type = "l", col = "yellow", ylab = "w(r)", add = TRUE)
legend(x = "topright", legend = c("(Q1, Q2) = \n(0.5, 0.75)", "(Q1, Q2) = \n(0.982, 1.473)", "(Q1, Q2) = \n(1.5, 2.25)", "(Q1, Q2) = \n(2, 3)"), 
       lty = 1, lwd = 2, cex = 1, y.intersp = 1.3,
       col = c("red", "black", "orange", "yellow"))
