
library(nimble)

xObs <- c( 1,  3,  5,  7, 10, 15, 20)    ## input
yObs <- c(10, 15, 18, 14, 20, 13, 19)    ## input
xPred <- 1:20                            ## input

nObs <- length(xObs)
nPred <- length(xPred)
f <- function(xi, xj) (xi-xj)^2
diffOO <- outer(xObs,  xObs,  f)
diffPP <- outer(xPred, xPred, f)
diffPO <- outer(xPred, xObs,  f)

gpPred <- nimbleFunction(
    setup = function(model, params) {
        calcNodes <- model$getDependencies(params, determOnly = TRUE)
        nPred <- dim(model$SigPP)[1]
        E <- array(0, c(nPred, 1))
        C <- array(0, c(nPred, nPred))
    },
    run = function(samples = double(2)) {
        E <<- E * 0
        C <<- C * 0
        nSamples <- dim(samples)[1]
        for(i in 1:nSamples) {
            values(model, params) <<- samples[i,]
            calculate(model, calcNodes)
            intermediate <- model$SigPO %*% inverse(model$SigOO)
            Etemp <-               intermediate %*% model$yObs
            Ctemp <- model$SigPP - intermediate %*% t(model$SigPO)
            E <<- E + Etemp
            C <<- C + Ctemp
        }
        E <<- E / nSamples
        C <<- C / nSamples
    },
    methods = list(
        getE = function() { returnType(double(1)); return(E[,1]) },
        getC = function() { returnType(double(2)); return(C[, ]) }
    )
)

## this is *not* as elegant as I had hoped.
## still requires repeating (essentially) the same code three times.
## I discovered some limitations of NIMBLE while trying other approaches.
code <- nimbleCode({
   rho ~ dgamma(10, 1)
   sigGP ~ dunif(0, 1e5)
   sigOE ~ dunif(0, 1e5)
   SigOO[1:nObs, 1:nObs ] <- sigGP^2 * exp(-1/2 * diffOO[1:nObs, 1:nObs ] / rho^2) + sigOE^2 * IOO[1:nObs, 1:nObs ]
   SigPP[1:nPred,1:nPred] <- sigGP^2 * exp(-1/2 * diffPP[1:nPred,1:nPred] / rho^2) + sigOE^2 * IPP[1:nPred,1:nPred]
   SigPO[1:nPred,1:nObs ] <- sigGP^2 * exp(-1/2 * diffPO[1:nPred,1:nObs ] / rho^2)
   yObs[1:nObs] ~ dmnorm(mu[1:nObs], cov = SigOO[1:nObs,1:nObs])
})
constants <- list(nObs=nObs, nPred=nPred, diffOO=diffOO, diffPP=diffPP, diffPO=diffPO, IOO=diag(nObs), IPP=diag(nPred), mu=rep(0,nObs))
data <- list(yObs = yObs)
inits <- list(rho = 1, sigGP = 1, sigOE = 1)
Rmodel <- nimbleModel(code, constants, data, inits)
params <- Rmodel$getNodeNames(topOnly = TRUE)

spec <- configureMCMC(Rmodel, nodes = NULL)
## NOTE: the next line shouldn't work for you, since the MCMC API has changed.
## you can rebuild from source off nimble branch 'devel'.
## otherwise, using last public release of NIMBLE (0.3-1), use this line instead:
## spec$addSampler('RW_block', control = list(targetNodes = params))
spec$addSampler(params, 'RW_block')
## we can debate about univariate vs. block sampling at some point
Rmcmc <- buildMCMC(spec)
Rpred <- gpPred(Rmodel, params)

Cmodel <- compileNimble(Rmodel)
Cmcmc  <- compileNimble(Rmcmc, project = Rmodel)
Cpred  <- compileNimble(Rpred, project = Rmodel)

set.seed(0)
system.time(Cmcmc$run(100000))
samples <- as.matrix(Cmcmc$mvSamples)
if(!identical(params, colnames(samples))) stop('problem')
system.time(Cpred$run(samples))
E <- Cpred$getE()
C <- Cpred$getC()
E
diag(C)


plot(xObs, yObs, type='b', xlim=range(c(xObs,xPred)), ylim=range(c(yObs,E)))
points(xPred, E, pch=19, col='red')


## wield this new weapon carefully!
## happy to discuss data & applications anytime.








