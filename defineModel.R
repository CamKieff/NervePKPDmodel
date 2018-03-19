#combine the different ODE defined models that can exist for the system

#simple vs. complex ACH model
#one NT vs two NTs vs more NTs(?)
#summation method of model if more than two neurotransmitters are used

#ODEinits for one drug vs two drugs
require(RxODE)

defineModel <- function(ACH_mod="simple", unk_mod="none", effect_mod = "oneNT"){
  ach_models <- list(
    none = "",
    simple = "d/dt(depot2) = -KA2*depot2;
            d/dt(centr2) = KA2*depot2 - KE2*centr2;
    ",
    complex = "d/dt(depot2) = -KA2*depot2 + KA2*KA2f*depot2*depot2/(depot2 + IC501);
             d/dt(centr2) = KA2*depot2 - KA2*KA2f*depot2*depot2/(depot2 + IC501) - KE2*centr2 + KE2*KE2f*centr2*centr2/(centr2 + IC502);
    ")

  unk_models <- list(
    none = "",
    simple = "d/dt(depot1) = -KA1*depot1;
              d/dt(centr1) = KA1*depot1 - KE1*centr1;
    "
  )

  effect_models <- list(
    none = "",
    oneNT = "eff2 = centr2/(EC502+centr2);",
    twoNT_1 = "effa = max1*centr1/(EC501 + centr1);
            effb = centr2/(EC502 + centr2);
            eff2 = effb - effa + (effa * effb);",
    twoNT_2 = "eff1 = max1*centr1/(EC501 + centr1);
            eff2 = (1-eff1)*centr2/(EC502 + centr2);"
  )
  if(unk_mod == "none"){
    effect_mod <- "oneNT"
  }

  ODE1 <- paste(unk_models[unk_mod], ach_models[ACH_mod], effect_models[effect_mod])

  mod1 <- RxODE(model = ODE1) #Compile Model

  if(unk_mod == "none"){
    ODEinits <- c(depot2 = 0, centr2 = 0) #initial conditions for single neurotransmitter
    NTnum <- 1
  } else {
    ODEinits <- c(depot1 = 0, centr1 = 0, depot2 = 0, centr2 = 0) #initial conditions for two neurotransmitters
    NTnum <- 2
  }

  return(list(mod1, ODEinits, NTnum))
}
