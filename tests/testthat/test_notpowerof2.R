test_that("not a power of 2 smash has similar results as a power of 2 case",{
    
  #' @title simulate gaussian or poisson distributed data
  #' @param n the number of data
  #' @param type a string that is either "gaus" or "poiss"
  simulate = function(n, type){
    t=1:n/n
    spike.f = function(x) (0.75*exp(-500*(x-0.23)^2) +
        1.5*exp(-2000*(x-0.33)^2) +
        3*exp(-8000*(x-0.47)^2) +
        2.25*exp(-16000*(x-0.69)^2) +
        0.5*exp(-32000*(x-0.83)^2))
    mu.s = spike.f(t)
    if (type=="gaus"){
      mu.t = (1+mu.s)/5
      var.fn = (0.0001 + 4*(exp(-550*(t-0.2)^2) +
                            exp(-200*(t-0.5)^2) +
                            exp(-950*(t-0.8)^2)))/1.35
      rsnr=sqrt(5)
      sigma.t=sqrt(var.fn)/mean(sqrt(var.fn))*sd(mu.t)/rsnr^2
      set.seed(1)
      X.s=rnorm(n,mu.t,sigma.t)
    } else if(type=="poiss"){
      mu.t <- 0.01 + mu.s
      set.seed(1)
      X.s <- rpois(n,mu.t)
    } else {
      stop("Type should be either gaus or poiss.")
    }
    return(X.s=X.s)
  }

  # data 1 simulates 1000 data points
  #
  # data 2 includes data 1 and then has extra 24 data points, i.e.,
  # the number of data 2 is 1024
  data1_gaus = simulate(1000, "gaus")
  data2_gaus = c(data1_gaus, simulate(24, "gaus"))
  data1_poiss = simulate(1000, "poiss")
  data2_poiss = c(data1_poiss, simulate(24, "poiss"))

  expect_equal(smash(data2_gaus)[1:1000], smash(data1_gaus), tolerance=0.02)
  expect_equal(smash(data2_poiss)[1:1000], smash(data1_poiss), tolerance=0.02)
})
