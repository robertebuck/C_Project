setwd('C:/Users/User/source/repos/WbsProject')


# Plotting paths

NSIM <- 100
NT <- 1000

#  <- as.matrix(read.csv('paths_euler.csv', header = F))
# paths_exact1 <- as.matrix(read.csv('paths_exact1.csv', header = F))
# paths_exact2 <- as.matrix(read.csv('paths_exact2.csv', header = F))

# matplot(x=c(0:NT), y=paths_euler, type = "l",pch=1,col = 1:NSIM, xlab = "Number of Time steps", ylab = "Price")
# matplot(x=c(0:NT), y=paths_exact1, type = "l",pch=1,col = 1:NSIM, xlab = "Number of Time steps", ylab = "Price")
# matplot(x=c(0:NT), y=paths_exact2, type = "l",pch=1,col = 1:NSIM, xlab = "Number of Time steps", ylab = "Price")


# Ploting price accross spot values

# price_euler <- read.csv('Euler_price_range.csv')
# price_exact1 <- read.csv('Exact1_price_range.csv')
# price_exact2 <- read.csv('Exact2_price_range.csv')

# Plot of European options
# matplot(x=price_euler[,1], y=cbind(price_euler[,2], price_euler[,6]), lwd = c(1,2), type = "l",pch=1,col = c(1,3), xlab = "Underlying price", ylab = "Option Price")
# matplot(x=price_exact1[,1], y=cbind(price_exact1[,2], price_exact1[,6]), lwd = c(1,2), type = "l",pch=1,col = c(1,3), xlab = "Underlying price", ylab = "Option Price")
# matplot(x=price_exact2[,1], y=cbind(price_exact2[,2], price_exact2[,6]), lwd = c(1,2), type = "l",pch=1,col = c(1,3), xlab = "Underlying price", ylab = "Option Price")

# legend("topleft",c("Simulated","Black Scholes"),col=c(1,3),lwd=2,lty=1)
# legend("topright",c("Simulated","Black Scholes"),col=c(1,3),lwd=2,lty=1)

# Plot of Asian options
# matplot(x=price_euler[,1], y=price_euler[,2], type = "l",pch=1, col = 1, xlab = "Underlying price", ylab = "Option Price")
# matplot(x=price_exact1[,1], y=price_exact1[,2], type = "l",pch=1,col = 1, xlab = "Underlying price", ylab = "Option Price")
# matplot(x=price_exact2[,1], y=price_exact2[,2], type = "l",pch=1,col = 1, xlab = "Underlying price", ylab = "Option Price")
