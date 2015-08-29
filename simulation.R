# The constants from Gibbons' book
# Returns 4 values of b, depending on the n, k, and shape parameter values.
# The four values given are for pcs=0.75, pcs=0.90, pcs=0.95, and pcs=0.99 respectively.
gibbonsConstant <- function(n=10, k=4, shapeparam=1) {
   gibbons <- dget("gibbons.out")
   return(gibbons[gibbons$n == n & gibbons$k == k & gibbons$shapeparam == shapeparam,]$b)
}

# Turned into a function for added flexibility from the R interface
simulation <- function(k=4, n=10, rep=1000, delta=0.005) {
   n0 <- 0
   t0 <- 0
   pcs <- 1:k
   integrand <- 1:k
   p <- 1:k
   t <- 1:k
   lam <- 1:k
   b <- 1:4
   PCSGupta <- rep(0,4)
   EGupta <- rep(0,4)
   guptaS <- rep(0,4)
   igupta <- rep(0,4)
   numberS <- rep(0, 250)
   indicator <- rep(0, 250)
   N <- rep(0, 250)
   NN <- rep(0, 4)
   II <- rep(0, 4)
   I <- rep(0, 5)
   ProbCS <- rep(0, 250)
   ExpectedS <- rep(0, 250)
   
   b <- gibbonsConstant(n=10, k=4, shapeparam=1) # Once gibbons.out is updated, delete the defaults
   
   tmax <- 0
   
   for(z in 1:rep) {
      lam[4] <- lam[3] <- lam[2] <- lam[1]*(1+delta)
      
      lambdamin <- min(lam)      
      t <- rep(0,k)
       
      for(i in 1:k) {
         t[i] <- rgamma(1, n, rate=lam[i]) # sample of 1 with shape parameter n
      }
      
      tmax <- max(t)
      for(i in 1:4) {
         for(j in 1:k) {
            if(t[j]>b[i]*tmax) {
               guptaS[i] <- guptaS[i]+1
               if(lam[j]==lambdamin) {
                  igupta[i] <- igupta[i]+1
               }
            }
         }
         NN[i] <- guptaS[i]
         II[i] <- igupta[i]
      }
          
      integrandt1 <- function(x) {
         (1-pgamma(x, n+n0+1, rate=t[2]+t0)) * 
         (1-pgamma(x, n+n0+1, rate=t[3]+t0)) * 
         (1-pgamma(x, n+n0+1, rate=t[4]+t0)) * (t[1]+t0)^(n+n0+1) * 
         x^(n+n0) * exp(-x * (t[1] + t0))/factorial(n+n0)
      }
      
      p <- integrate(integrandt1, lower=0, upper=Inf)
      pcs[1] <- p$value
      
      integrandt2 <- function(x) {
         (1-pgamma(x, n+n0+1, rate=t[1]+t0)) *
         (1-pgamma(x, n+n0+1, rate=t[3]+t0)) *
         (1-pgamma(x, n+n0+1, rate=t[4]+t0)) * (t[2]+t0)^(n+n0+1) *
         x^(n+n0) * exp(-x * (t[2] + t0))/factorial(n+n0)
      }
   
      p <- integrate(integrandt2, lower=0, upper=Inf)
      pcs[2] <- p$value
   
      integrandt3 <- function(x) {
         (1-pgamma(x, n+n0+1, rate=t[1]+t0)) *
         (1-pgamma(x, n+n0+1, rate=t[2]+t0)) *
         (1-pgamma(x, n+n0+1, rate=t[4]+t0)) * (t[3]+t0)^(n+n0+1) *
         x^(n+n0) * exp(-x * (t[3] + t0))/factorial(n+n0)
      } 
      
      p <- integrate(integrandt3, lower=0, upper=Inf)
      pcs[3] <- p$value
   
      integrandt4 <- function(x) {
         (1-pgamma(x, n+n0+1, rate=t[1]+t0)) *
         (1-pgamma(x, n+n0+1, rate=t[2]+t0)) *
         (1-pgamma(x, n+n0+1, rate=t[3]+t0)) * (t[4]+t0)^(n+n0+1) *    # This line is semi-shady, because there was a t[3] instead of the expected t[4] in the base term.
         x^(n+n0) * exp(-x * (t[4] + t0))/factorial(n+n0)
      }
      
      p <- integrate(integrandt4, lower=0, upper=Inf)
      pcs[4] <- p$value
      for(c in 1:250) {
         for(j in 1:k) {
            if(pcs[j]>(1/(1+c))) {
               numberS[c] <- numberS[c] + 1
               if(lam[j] == lambdamin) {
                  indicator[c] <- indicator[c] + 1
               }
            }
         }
         N[c] <- numberS[c]
         I[c] <- indicator[c]
         ProbCS[c] <- I[c]/rep
         ExpectedS[c] <- N[c]/rep
      }
   }
   
   for(j in 1:4) {
      PCSGupta[j] <- II[j]/rep
      EGupta[j] <- NN[j]/rep
   }
   
#   for(c in 1:250) {
#      cat("  c=", c)
#      cat("  P(CS)=", ProbCS[c])
#      cat("  E[S]=", ExpectedS[c], "\n")
#   }

   out1 <- character(4)
   outcount <- 1
   for(c in 1:250) {
      if(c==16 || c==56 || c==124 || c==200) {
         out1[outcount] <- paste("c=", c, "  P(CS)=", ProbCS[c], "  E[S]=", ExpectedS[c], sep="")
         outcount <- outcount + 1
      }
   }
   out2 <- character(4)
   for(j in 1:4) {
      out2[j] <- paste("j=", j, "  P(CSGupta)=", PCSGupta[j], "  E[S_Gupta]=", EGupta[j], sep="")
   }
   #pcs is included to confirm that it's evaluated to a non-1 value
   output <- list(out1=out1, out2=out2, pcs=pcs)
   print(output)
}
