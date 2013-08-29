load("130827-01 - Workspace Modell Monster KM DIS Quantos0001.rData")

library(mBioSpec)
dyn.load("C:/micro-biolytics/Firmen und Projekte/Q-Food/QFood002C - Monster Energydrink/Modell DIS/mBioerx.dll")


setwd("C:/micro-biolytics/Firmen und Projekte/Q-Food/QFood002E - diverse Softdrinks/R-Daten")
load("130827 - DivDrinks EU Validierungsmessungen für RunDIS.rlist")

setwd("C:/micro-biolytics/Firmen und Projekte/Q-Food/QFood002D - Konservierungsmittel/R-Daten")
load("130827 - Monster EU Validierungsmessungen für RunDIS.rlist")






Results <- list()

## Run mit FI.Diff = as.WaveNum(X3001:X2791, X1900X991), Calc.pH = TRUE, Calc.pHDiff = FALSE und DI.Mod = DI.Mod
   nam  <- "Run01"

   for (i in names(Monster)) {
      cat(round(match(i, names(Monster)) / length(Monster), 2) * 100, "% ... calculating data",rep(" ", 5), "\r", sep = "") ; flush.console()
      Results$Monster[[nam]][[i]] <- Run.DIS(Monster[[i]], Result.Admin = FALSE, Calc.pH = TRUE, Calc.pHDiff = FALSE)
   }

   Results$Monster[[nam]][["21496678"]] <- NULL

   ## Result.Admin = TRUE -----------
      r1.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Result$Result1))
      r2.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Result$Result2))
      r3.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Result$Result3))
      pHd.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Result$pH.Diff))
      coeff.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Coef))
      error.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Error))

   ## Result.Admin = FALSE -----------
      r.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Result))
      error.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Error))
      schidi.run2 <- t(sapply(Results$Monster[[nam]], function(x) x$Schidi))
