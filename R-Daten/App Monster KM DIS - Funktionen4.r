################################################################################################################################
## App für Monster DIS Konservierungsmittel
################################################################################################################################

Run.DIS <- function (SPEC, Result.Admin = FALSE, Calc.pH = TRUE, Calc.pHDiff = FALSE) {

   ## SPEC-Liste neu ordnen ------
      DAT <-  list()
      for (i in sort(names(SPEC), decreasing = TRUE)) {
         nam <- paste("IDX", match(i, sort(names(SPEC), decreasing = TRUE)), sep = "")
         DAT[[nam]] <- SPEC[[i]]
      }
            
   ## Freigabe relevanter Analysenwerte -------------
      freigabe.ges <- "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure"
      Result <- rep(-1, length(grep(freigabe.ges, unique(comp.dat$CompName))) + 1) ; names(Result) <- cas[grep(freigabe.ges, unique(cas$CompName)), "AnalytID "]  
            
   ## Methode einlesen und Auswerteparameter zuordnen-----------
      method <- as.numeric(DAT$IDX1$METHOD)
      MethodPar <- Select.MethodPar(method)

   ## Spektrenvorbereitung und Checks -------------
      AB.list <- list()           
      for (x in names(DAT)) {
         AB.list[[x]] <- Data.PreTreat(DAT[[x]], MethodPar$BL, MethodPar$CompAnalyt)
      }     
                
   ## Datenauswertung AB-Spektrum Original ---------
      ## DOSIM Berechnung ---------
         freigabe1 <- MethodPar$Freigabe$IDX1
         comp.rel1 <- grep(freigabe1, unique(comp.dat$CompName), value = TRUE)
         
         if (AB.list[["IDX1"]]$Do.Calc) {
            result1 <- Calc.DOSIM_AB(AB.list[["IDX1"]]$AB, comp.rel1, MethodPar$FI, MethodPar$BL, MethodPar$CompAnalyt, MethodPar$DI, Calc.pH)
            pH.low <- result1["pH"]         
         } else {
            result1 <- rep(-1, length(comp.rel1) + 2)
            names(result1) <- c(comp.rel1, "pH", "Gesamtzucker")
         }
         
         if (Result.Admin) {
            Result.All <- list()
            Result.All[["Result1"]] <- result1
         }
            
      ## Überarbeitung relevanter Analysenwerte ------      
         result1 <- result1[grep("pH", names(result1), invert = TRUE)]
         result1.cas <- Mod.Results(result1, freigabe1, cas, AB.list[["IDX1"]]$Do.Calc, lod.analyt)
         
   ## Datenauswertung AB-Spektrum moduliert ---------
      if (any("IDX2" == names(AB.list)) && (method >= 200)) {
         if (Result.Admin) {
            comp.rel3 <- comp.rel1
         } else {
            comp.rel3 <- comp.rel1[1]
         }
         
         ## DOSIM Berechnung ---------         
            if (AB.list[["IDX2"]]$Do.Calc) {
               dat <- AB.list[["IDX2"]]$AB
               result3 <- Calc.DOSIM_AB(dat, comp.rel3, MethodPar$FI, MethodPar$BL, MethodPar$CompAnalyt, MethodPar$DI.Mod, Calc.pH)
               pH.high <- result3["pH"]
            } else {
               result3 <- rep(-1, length(comp.rel3) + 2)
               names(result3) <- c(comp.rel3, "pH", "Gesamtzucker")
            }

            if (Result.Admin) {
               Result.All[["Result3"]] <- result3
            }
                  
         ## Überarbeitung relevanter Analysenwerte ------      
            #result3 <- result3[grep("pH", names(result3), invert = TRUE)]
            #result3.cas <- Mod.Results(result3, freigabe1, cas, AB.list[["IDX1"]]$Do.Calc, lod.analyt) 
      }
                     
   ## Datenauswertung Diff-Spektrum ---------
      if (any("IDX2" == names(AB.list)) && (method >= 200)) {
         freigabe2 <- MethodPar$Freigabe$IDX2
         comp.rel2 <- grep(freigabe2, unique(comp.dat$CompName), value = TRUE)
      
         ## DOSIM Berechnung ---------
            if (AB.list[["IDX2"]]$Do.Calc) {
               temp <- Calc.DOSIM_DIFF(AB.list[["IDX2"]]$AB, AB.list[["IDX1"]]$AB, comp.rel2, MethodPar$FI.Diff, MethodPar$BL.Diff, MethodPar$CompAnalyt.Diff, MethodPar$DI.Diff, Calc.pHDiff, pH.low, pH.high)
               result2 <- temp$konz.ck
               dilution <- result2["Dilution"]               
            } else {
               result2 <- rep(-1, length(comp.rel2))
               names(result2) <- c(comp.rel2, "Dilution")
            }

            if (Result.Admin) {
               Result.All[["Result2"]] <- result2            
               Result.All[["pH.Diff"]] <- temp$pH
            }
                              
         ## Überarbeitung relevanter Analysenwerte ------      
            result2 <- result2[grep("Dilution", names(result2), invert = TRUE)]
            result2.cas <- Mod.Results(result2, freigabe2, cas, AB.list[["IDX2"]]$Do.Calc, lod.analyt.diff)
      }
                                         
   ## Zusammenfassung Ergebnisse ----------
      Schidi <- Result.Schidi(AB.list, DAT, SPEC)      
      Error <- apply(sapply(AB.list, function(x) x$Error), 1, max)
      Result[names(result1.cas)] <- result1.cas
      
      if (any("IDX2" == names(AB.list)) && (method >= 200)) {
         Result[names(result2.cas)] <- result2.cas
      } 

   ## Check auf richtige App-Auswahl für Probe ------------
      error.app <- Check.AppSelect(Result, MethodPar$Comp.Limits)
      Error <- c(Error, error.app)
      if (error.app == 60) Result[] <- -1      

   ## Check auf richtige Verdünnung der modulierten Probe ------------            
      Error["Dilution"] <- 0
      if (any("IDX2" == names(AB.list)) && (method >= 200)) {     
         if ((dilut.max.warn < dilution) || (dilut.min.warn > dilution))                  Error["Dilution"] <- 71
         if ((dilut.max.error < dilution) || (dilut.min.error > dilution))                Error["Dilution"] <- 70
      }
      
   ## Check auf richtigen pH der modulierten Probe und der Originalprobe ------------            
      Error["Toleranz"] <- 0
      if (any("IDX2" == names(AB.list)) && (method >= 200)) {
         if (Calc.pH) {
            if ((pH.low.max.warn < pH.low) || (pH.low.min.warn > pH.low))                  Error["Toleranz"] <- 81
            if ((pH.high.max.warn < pH.high) || (pH.high.min.warn > pH.high))                Error["Toleranz"] <- 83
            #if ((pH.low.max.error < pH.low) || (pH.low.min.error > pH.low))                  Error["Toleranz"] <- 80
            #if ((pH.high.max.error < pH.high) || (pH.high.min.error > pH.high))                Error["Toleranz"] <- 82
         }
      }
                                                     
   ## Output -------
      if (Result.Admin) {
         if (any("IDX2" == names(AB.list)) && (method >= 200)) {
            Result.All$Result3[grep("pH", names(Result.All$Result3), invert = TRUE)] <- round(Result.All$Result3[grep("pH", names(Result.All$Result3), invert = TRUE)] / dilution, 3)
         }
         Coef <- sapply(AB.list, function(x) round(x$Coef, 5))               
         return(list(Result = Result.All, Schidi = Schidi, Coef = Coef, Error = Error, Method = method, Version = VersID))
      } else {
         return(list(Result = Result, Schidi = Schidi, Error = Error, Method = method, Version = VersID))
      }            
   }
   
   
################################################################################################################################
## Funktion für Run.DIS
################################################################################################################################

## Auswerteparameter nach Methode zuordnen-----------
   #Select.MethodPar <- function(method, MethodPar.List) {
         #nam <- paste("M", method, sep = "")
         #MethodPar <- MethodPar.List[[nam]]
      
      #return(MethodPar)
   #}                        
 
## Spektrenvorbereitung und Checks -------------
   Data.PreTreat <- function (DAT, BL, comp.analyt) {
      ## Eingabespektren konvertieren in mbio.data-Format ---------------
         AB <- Convert.Data(DAT)
         
      ## Make Compatible -------------
         # noch einfügen
            
      ## Schichtdicke berechnen ------------
         Schidi <- Calc.Schidi(AB)
         AB <- Schidi[[1]]
         schidi <- Schidi[[2]]
         schidi.bkpvalue <- Schidi[[3]]
   
      ## Atmosphärische Kompensation ----------
         AB <- AtmoKorr.DIS(AB)
         
      ## Ableiten Probenspektrum nach SDB-Vorgaben -------------
         AB <- savgol(AB, 1, 5, 2)[, c(3050:950),, FALSE]         
   
      ## Check valides Spektrum (Luft und CO2) ----------- 
         Error <- Check.AB(AB, comp.check, BL, comp.analyt, max.luft.warn, max.luft.error, max.CO2.warn, max.SqS.warn, max.SqS.error, max.WV.warn, max.WV.error, schidi.bkpvalue)
         Do.Calc <- Error$Do.Calc
         Error$Do.Calc <- NULL
         Coef <- Error$Coef      
               
      ## Ableiten Probenspektrum nach App-Vorgaben ---------
         if (Do.Calc) {
            AB <- savgol(AB, Abl.App, Gl.App, 2)[, c(3020:2450, 2250:980),, FALSE]
         }
                  
      return(list(AB = AB, Do.Calc = Do.Calc, Error = Error$Error, Schidi = schidi, Coef = Coef))
   }
      
      
## Einlesen Proben- und Schidispektren -----------
   Convert.Data <- function(DAT) {
      AB <- list()
   
      ## Probenspektrum ---------
         AB$Sam <- ToMbioData(DAT$SAM, DAT$PROPERTY, "Sample")
   
      ## Schidispektren ------------         
         schidi <- sapply(DAT$SCHIDI, function(x) x)
         schidi <- apply(schidi, 1, mean)

         AB$Schidi <- ToMbioData(schidi, DAT$PROPERTY, "Schidi")
      return(AB)
   }

   ToMbioData <- function(Data, property, Name) {   
      PROPERTY <- vector()
     
      PROPERTY["NPT"] <- property[1]
      PROPERTY["FXV"] <- property[3]                           
      PROPERTY["LXV"] <- property[2]                           
      PROPERTY["RES"] <- (PROPERTY["FXV"] - PROPERTY["LXV"]) / (PROPERTY["NPT"] - 1) 
      PROPERTY["DXU"] <- 1
              
      PARAM <- as.data.frame(matrix(c(1, 1), ncol = 2))
      HISTORY <- list()
      NOTE <- ""
      TYPE <- "AB"
      Data <- matrix(Data, ncol = length(Data), dimnames = list(Name, as.character(1:length(Data))))
                    
      dat.mbio <- new("mbio.spec.data", Data, PARAM = PARAM, PROPERTY = PROPERTY, HISTORY = HISTORY, NOTE = NOTE, TYPE = TYPE)
      
      return(dat.mbio)
   }

## Schichtdicke berechnen ----------
   Calc.Schidi <- function(AB) {
      
      ## Parameter festlegen ---------
         Gl <- 25
         Abl <- 1
         Allow.SchidiBackUp <- TRUE
         schidi.bkpvalue <- 0
         value.backup <- 9.1
         
      ## Schidi über FFT-Fit-------------------
         waves.fft <- wavenumbers(AB$Sam[, 7450:3900,, FALSE])
         schidi.fft <- NULL
         
         dat.fft <- savgol(AB$Schidi, Abl, Gl, 2)
         dat.fft <- dat.fft[, 7450:3900,, FALSE]  
   
         fft.result <- find_thick(dat.fft@.Data, wavenumbers(dat.fft), find_thick.par = list(NumberOfZeros = 2.5e5, Refrcorrect = TRUE))
         schidi.fft <- c(fft.result$thick_found, fft.result$SNR, fft.result$MMR)
         names(schidi.fft) <- c("Schidi.FFT", "SNR", "MMR")
         
         value.schidi <- round(schidi.fft["Schidi.FFT"], 3)

      ## Backupwert falls Schichtdickenbestimmung nicht fitted ------------
         if (Allow.SchidiBackUp) {
            factor.min <- 0.9
            factor.max <- 1.1
            check.value <- (value.schidi > value.backup * factor.max) || (value.schidi < value.backup * factor.min)
            if (is.na(value.schidi) || check.value) {
               value.schidi <- value.backup
               schidi.bkpvalue <- 1
            }
         }

      ## Einheitsspektrum berechnen -----------
         AB$Sam@.Data <- AB$Sam@.Data / value.schidi
         
      return(list(AB$Sam, value.schidi, schidi.bkpvalue))
   }

## atmosphärische Kompensation ----------------
   AtmoKorr.DIS <- function(dat) { 
   
      ## Parameter --------------
         Abl.Korr = 3
         res <- paste("RES", round(dat@PROPERTY["RES"], 0) * 2, sep = "")
      
      ## Auswahl Korrekturspektren -----------
         korrspec.wd <- AtmoKorrSpecs$KorrSpecs[AtmoKorrSpecs$NumbersPC, ]
         
      ## Datenvorbereitung Korrekturfit ------------       
         korrspec.wd.der <- korrspec.wd            
         dat.der <- dat            
         
         for (x in 1:Abl.Korr) {
            korrspec.wd.der <- savgol(korrspec.wd.der, 1, 5, 2)
            dat.der <- savgol(dat.der, 1, 5, 2)
         }      
                        
      ## Korrektur-Fit -----------                 
         comp <- t(korrspec.wd.der[, c(3980:3800, 1950:1820),, FALSE]@.Data)
         dat.fit <- t(dat.der[, c(3980:3800, 1950:1820),, FALSE]@.Data)
      
         ## lm-fit -------------
            fit.lm <-  lm(dat.fit ~ comp)
            coef.fit <- t(as.data.frame(coef(fit.lm))[2, ])
            dat@.Data <- dat@.Data  - coef.fit %*% korrspec.wd
     
         dat@PARAM$Coef_WDKorr <- as.vector(coef.fit)
      
      return(dat)
   }
   
## Checks ------------
   Check.AB <- function(AB, comp.check, BL, comp.analyt, max.luft.warn, max.luft.error, max.CO2.warn, max.SqS.warn, max.SqS.error, max.WV.warn, max.WV.error, schidi.bkpvalue) {
      ## Error Vector erstellen ---------
         Do.Calc <- TRUE
         Error <- rep(0, 5)
         names(Error) <- c("Luft", "WaterVapor", "SqS", "CO2", "Schidi.Value")
                                 
      ## Berechnung Modell ------------
         dat <- AB[, c(3020:2450, 2250:980),, FALSE]      

         ckname <- "Glukose"                                            ## Glukose ohne Relevanz,nur um richtige Auswahl der Komponenten zu treffen
         sel.comp <- rownames(comp.dat)[comp.dat$CompName %in% names(comp.analyt)[comp.analyt[grep(ckname, comp.analyt$CompName), ] == "x"]]
         bl <- as.numeric(BL[, c("Intercept", sel.comp)])
         bu <- rep(Inf, length(bl))
         bu[grep("PCA.Loading_LCS.1", names(BL[, c("Intercept", sel.comp)]))] <- abs(bl[grep("PCA.Loading_LCS.1", names(BL[, c("Intercept", sel.comp)]))])

         comp.check <- comp.check[, sel.comp]

         fit.bvls <- bvls(cbind(rep(0.0005, nrow(comp.check)), comp.check), t(dat@.Data), bl, bu)
         names(fit.bvls$x) <- c("Intercept", colnames(comp.check))

      ## Luft ---------         
         coef.luft <- fit.bvls$x[grep("Wasserspektrum", names(fit.bvls$x))] 

      ## Wasserdampf -------------
         #coef.wv <- fit.bvls$x[grep("PCA.Loading_LCS.2", names(fit.bvls$x))]
         coef.wv <- 0
         
      ## SqS ----------
         sqs <- sum(residuals(fit.bvls)^2) * 1e9

      ## CO2 -------------
         dat <- savgol(AB, 1, 5, 2)        
         coef.CO2 <- apply(dat[, 2450:2250,, FALSE], 1, max)         
                                 
         #fit.lm <- lm(t(AB@.Data) ~ comp.check)
         #coef.lm <- coef(fit.lm)[grep("Wasserspektrum", names(coef(fit.lm)))]     
            
      ## Error zuweisen ------------
         if (coef.luft > max.luft.warn)               Error["Luft"] <- 11
         if (coef.luft < -max.luft.warn)              Error["Luft"] <- 13
         if (coef.luft < -max.luft.error)          {  Error["Luft"] <- 10 ; Do.Calc <- FALSE }
         if (coef.luft > max.luft.error)           {  Error["Luft"] <- 12 ; Do.Calc <- FALSE }

         #if (abs(coef.wv) > max.WV.warn)              Error["WaterVapor"] <- 41
         #if (abs(coef.wv) > max.WV.error)          {  Error["WaterVapor"] <- 40 ; Do.Calc <- FALSE }

         if (sqs > max.SqS.warn)                      Error["SqS"] <- 21
         if (sqs > max.SqS.error)                  {  Error["SqS"] <- 20 ; Do.Calc <- FALSE }
         
         if (coef.CO2 > max.CO2.warn)                 Error["CO2"] <- 31
         
         if (schidi.bkpvalue == 1)                    Error["Schidi.Value"] <- 51
                  
         out <- c(coef.luft, coef.wv, sqs, coef.CO2, schidi.bkpvalue)
         names(out) <- names(Error)
                                    
      return(list(Error = Error, Do.Calc = Do.Calc, Coef = out))
   }

## Datenauswertung mittels DOSIM ---------------------
   ## AB-Spektren ----------
      Calc.DOSIM_AB <- function(AB, comp.rel, FitInterval, BL, comp.analyt, DI, Calc.pH) {
         Result <- list()
         Temp <- list()
         
         ## Konvertieren AB-Spektrum für Fit ------------
            fit.dat <- as.data.frame(AB, TRUE)
   
         ## Parameter Fit allgemein ------------
            #FitInterval <- as.WaveNum(X3001:X2791, X1900:X991)
            FitInterval <- FitInterval[FitInterval %in% colnames(fit.dat)[wrange]]
            Gl <- 0
            Abl <- 0
            weight.pot = 1/3
            weight.max = 10
            weight.min = 0.01         
             
         ## BVLS-Fit ----------------          
            for (ckname in comp.rel) {               
               sel.comp <- rownames(comp.dat)[comp.dat$CompName %in% names(comp.analyt)[comp.analyt[grep(ckname, comp.analyt$CompName), ] == "x"]]
               bl <- as.numeric(BL[, c("Intercept", sel.comp)])
               bu <- rep(Inf, length(bl))
               bu[grep("PCA.Loading_LCS.1", names(BL[, c("Intercept", sel.comp)]))] <- abs(bl[grep("PCA.Loading_LCS.1", names(BL[, c("Intercept", sel.comp)]))])
                           
               comp.med <- comp.dat[sel.comp, ]
   
               dat.imp <- DI[[ckname]][[as.character(Abl + Abl.App)]]
               dat.imp <- SavGol(dat.imp, 3, 2, 0)
               dat.imp <- abs(as.vector(t(dat.imp[, FitInterval])))^(weight.pot)
               dat.imp[dat.imp < weight.min] <- weight.min
               dat.imp[dat.imp > weight.max] <- weight.max
               
               Weights <- sqrt(dat.imp)
                        
               ## Baseline-Spektren offset-korrigiert -----------
                  #bas.mwkorr <- bas[, FitInterval] - apply(bas[, FitInterval], 1, mean)
   
               Temp[[ckname]] <- DOSIM.BVLS_AB(fit.dat, Abl, Gl, FitInterval, Weights, bl, bu, comp.med, Calc.pH)             
               Result$BVLS[[ckname]] <- Temp[[ckname]]$Konz.BVLS[ckname]
            }
   
         ## Auslesen aus Result-Liste ---------
            konz.ck <- sapply(Result$BVLS, function(x) as.vector(t(x)))
            pH <- round(mean(sapply(Temp, function(x) x$pH), na.rm = TRUE), 2)
            konz.ck["pH"] <- pH
   
         ## Berechnung Gesamtzucker ------------
            konz.ck["Gesamtzucker"] <- sum(konz.ck[c("Glukose", "Fruktose", "Saccharose")])
   
            return(konz.ck)
      }

   ## Diff-Spektren ----------
      Calc.DOSIM_DIFF <- function(Diff, AB, comp.rel, FitInterval, BL.Diff, comp.analyt.diff, DI.Diff, Calc.pHDiff, pH.low, pH.high) {
         Result <- list()
         Temp <- list()
         PH <- list()
         
         ## Konvertieren Diff-Spektrum für Fit ------------
            fit.dat <- as.data.frame(Diff, TRUE)
            ref.dat <- as.data.frame(AB, TRUE)
            ref.dat$CompID <- "XXXORIG"
            ref.dat$CompName <- "Original"
            ref.dat$Konz_Komp_mg_mL <- 1
            ref.dat$pH <- 3.55
            
         ## Parameter Fit allgemein ------------
            #FitInterval <- as.WaveNum(X3001:X2791, X1900:X991)
            FitInterval <- FitInterval[FitInterval %in% colnames(fit.dat)[wrange]] 
            Gl <- 0
            Abl <- 0
            weight.pot = 1/3
            weight.max = 10
            weight.min = 0.01
             
         ## BVLS-Fit ----------------          
            for (ckname in comp.rel) {               
               sel.comp <- rownames(comp.diff)[comp.diff$CompName %in% names(comp.analyt.diff)[comp.analyt.diff[grep(ckname, comp.analyt.diff$CompName), ] == "x"]]
               bl <- as.numeric(BL.Diff[, c("Intercept", sel.comp)])
               bu <- rep(Inf, length(bl))
               bu[grep("PCA.Loading_LCS.1", names(BL.Diff[, c("Intercept", sel.comp)]))] <- abs(bl[grep("PCA.Loading_LCS.1", names(BL.Diff[, c("Intercept", sel.comp)]))])
               
               comp.med <- comp.diff[sel.comp, ]
             
               ## Zufügen ref.dat-Parameter ----------
                  bl <- c(0, bl)                                    
                  bu <- c(Inf, bu)                                    
                  sel.comp <- c(rownames(ref.dat), sel.comp)
                  comp.med <- rbind(ref.dat[, colnames(comp.med)], comp.med)         
                  
               dat.imp <- DI.Diff[[ckname]][[as.character(Abl + Abl.App)]]
               dat.imp <- SavGol(dat.imp, 3, 2, 0)
               dat.imp <- abs(as.vector(t(dat.imp[, FitInterval])))^(weight.pot)
               dat.imp[dat.imp < weight.min] <- weight.min
               dat.imp[dat.imp > weight.max] <- weight.max
               
               Weights <- sqrt(dat.imp)
                        
               ## Baseline-Spektren offset-korrigiert -----------
                  #bas.mwkorr <- bas[, FitInterval] - apply(bas[, FitInterval], 1, mean)
   
               Temp <- DOSIM.BVLS_DIFF(fit.dat, Abl, Gl, FitInterval, Weights, bl, bu, comp.med, Calc.pHDiff, pH.low, pH.high)
               Result$BVLS[[ckname]] <- Temp$Konz.BVLS[ckname]
               Result$BVLS[["Dilution"]][[ckname]] <- Temp$Konz.BVLS["Original"]
               PH[[ckname]] <- Temp$pH                            
            }
           
            Result$BVLS[["Dilution"]] <- mean(sapply(Result$BVLS[["Dilution"]], function(x) as.vector(t(x))), na.rm = TRUE)
            pH <- round(sapply(PH, function(x) as.vector(t(x))), 2)
               
         ## Auslesen aus Result-Liste ---------
            konz.ck <- sapply(Result$BVLS, function(x) as.vector(t(x)))
  
            return(list(konz.ck = konz.ck, pH = pH))
      }

## Freigabe, Justierung und Überprüfung Nachweisgrenze relevanter Analysenwerte und Zuweisung AnalytID-Nummern ---------
   Mod.Results <- function(result, freigabe, cas, Do.Calc, lod.analyt) {
      #result[grep(freigabe, names(result), invert = TRUE)] <- NA
      result <- result[grep(freigabe, names(result))]

      ## Justierung der Ergebnisse und Überprüfung der Nachweisgrenze pro Analyt ----------
         if (Do.Calc) {
            for (i in names(result)) {
               result[i] <- round(result[i] * lod.analyt$Justierfaktor[match(i, lod.analyt$CompName)], 3)
              
               if (is.na(result[i]) || is.infinite(result[i])) {
                  result[i] <- -1
               } else {
                  if (result[i] < lod.analyt$LOD[match(i, lod.analyt$CompName)]) {
                     result[i] <- 0
                  }
               }
            }
         }
      
      ## Zuweisung AnalytID -------------
         result.cas <- sapply(names(result), function(x) cas$AnalytID[match(x, cas$CompName)])
         names(result) <- result.cas
     
      return(result)
   }
      
## Zusammenfassung Schidi und Namensänderung zurück zu SPEC -------------
   Result.Schidi <- function(AB.list, DAT, SPEC) {      
      Schidi <- sapply(AB.list, function(x) x$Schidi)
      for (i in sort(names(DAT), decreasing = TRUE)) {
         nam <- paste("SPEC", match(i, sort(names(DAT), decreasing = TRUE)) - 1, sep = "")
         names(Schidi) <- gsub(i, nam, names(Schidi))
      }
      
      return(Schidi)
   }      
      
## Check auf richtige App-Auswahl für Probe ------------
   Check.AppSelect <- function(Result, comp.limits) {
      error <- rep(0, 1)
      names(error) <- c("App.Check")
      error.count <- 0
      
      for(i in names(comp.limits)) {
         value <- Result[cas[grep(i , cas$CompName), "AnalytID "]]           
         if ((value > comp.limits[[i]]$Max) || (value < comp.limits[[i]]$Min)) error.count <- error.count + 1
      }
      
      if (error.count == 1) error["App.Check"] <- 61
      if (error.count >= 2) error["App.Check"] <- 60         
      
      return(error)
   }   
      
                        
################################################################################################################################
## DOSIM
################################################################################################################################

DOSIM.BVLS_AB <- function(fit.dat, Abl, Gl, FitInterval, Weights, bl, bu, comp.med, Calc.pH) {
   ## Fit (alle Proben) -------------
      #CompKonz.BVLS <- Do.BVLS(fit.dat, Abl, Gl, FitInterval, Weights, bl, bu)

      ## Zusammenstellen Komponenten -----------
         comp.list <- list()
         spec.offset <- as.data.frame(t(rep(0.00005, ncol(comp.med)))) ; rownames(spec.offset) <- "Intercept" ; colnames(spec.offset) <- colnames(comp.med)
         comp.list[[1]] <- as.data.frame(rbind(spec.offset, comp.med))
   
      ## Berechnung Fit ----------
         fit.bvls <- Fit.BVLS(fit.dat, comp.list, bas = NULL, Abl, Gl, FitInterval, Weights, bl, bu)

      ## Erstellen Fitergebnisse ---------
         dat.res <- as.data.frame(matrix(nrow = nrow(fit.dat), ncol = nrow(comp.med)))
         colnames(dat.res) <- rownames(comp.med) ; rownames(dat.res) <- rownames(fit.dat)
         dat.res$Intercept <- NA ; dat.res$SqSum <- NA
         #dat.res[, "SqSum"] <- fit.bvls$SqSum
   
         dat.res[, colnames(fit.bvls)] <- fit.bvls
         
         for (g in rownames(comp.med)) {
            if (!grepl("MCR.Loading", g)) {
               dat.res[, g] <- dat.res[, g] * comp.med[g, "Konz_Komp_mg_mL"]
            }
         }
      #CompKonz.BVLS <- dat.res
 
   ## Berechnung der Konzentrationen ---------
      dat.calc <- dat.res[, c(rownames(comp.med), "SqSum")]
      
      ID.dat <- c(levels(factor(comp.med$CompID)), "SqSum")
      dat.temp <- as.data.frame(t(dat.calc))
      dat.temp$CompID <- c(comp.med$CompID, "SqSum")
      dat.temp$CompName <- c(comp.med$CompName, "SqSum")
      konz <- SUMvers9(dat.temp, ID.dat, addChar = TRUE, File.ID = dat.temp$CompID, Counter = FALSE)$SUM
      rownames(konz) <- konz$CompName 

      ## Berechnung pH-Werte Monster Original aus Monster-Anteil der reinen Zitronensäure ------------
         if (Calc.pH) {
            calccomp.ph <- "Zitronensaeure"
            sel <- rownames(comp.med)[grep(calccomp.ph , comp.med$CompName)]
            E <- as.data.frame(cit.svd, TRUE)
            pks <- c(3.13, 4.5, 5.8)
                              
            if (!Gl == 0) {
               comp.der <- SavGol(comp.med[sel, ], Gl, 2, Abl)
               E <- SavGol(E, Gl, 2, Abl)
            } else {
               comp.der <- comp.med[sel, ]
            }
                                             
            mcr.konz <- as.matrix(dat.res[, sel])
            spec.dat <- mcr.konz %*% as.matrix(comp.der[, FitInterval])                                
            C.cit <- spec.dat %*% pinv(E[, FitInterval])
            C.cit[C.cit < 0] <- 0    

            pH <- calc.pH(C.cit, pks)
         } else {
            pH <- NA
         }
              
      ## Rückfit MCR-Loadings mit Komponentenspektren ----------------         
         mcr.num <- grep("MCR.Loading", rownames(comp.med))
      
         if (length(mcr.num) > 0) {
            mcr.names <- unique(comp.med[mcr.num, "CompName"])
            konz.mcr <- NULL
            compkonz.mcr <- NULL
                         
            for (x in mcr.names) {
               #konz.mcr <- matrix(nrow = length(mcr.names), ncol = rownames(dat.res)) ; rownames(konz.mcr) <- mcr.names ; colnames(konz.mcr) <- 
               sdb.med <- as.data.frame(complist.sdb[[x]], TRUE)
               sel <- rownames(comp.med)[grep(x , comp.med$CompName)]
               koef.mcr <- list()
               
               if (!Gl == 0) {
                  comp.der <- SavGol(comp.med[sel, ], Gl, 2, Abl)
                  sdb.med <- SavGol(sdb.med, Gl, 2, Abl)
               } else {
                  comp.der <- comp.med[sel, ]
               }
               
               mcr.konz <- as.matrix(dat.res[, sel])
               spec.dat <- mcr.konz %*% as.matrix(comp.der[, FitInterval])
                              
               mcr.fit <- lm(t(spec.dat) ~ t(sdb.med[, FitInterval]))
               
               #koef.mcr <- sapply(rownames(spec.dat), function(x) t(coef(mcr.fit))[x, ])
               koef.mcr <- as.matrix(coef(mcr.fit))
               konz.mcr <- rbind(konz.mcr, matrix(sdb.med$Konz_Komp_mg_mL, nrow = 1) %*% koef.mcr[-1, ])                 
               compkonz.x <- matrix(sdb.med$Konz_Komp_mg_mL * koef.mcr[-1, ]) ; rownames(compkonz.x) <- rownames(sdb.med)
               compkonz.mcr <- rbind(compkonz.mcr, compkonz.x)                     
            }
            
            rownames(konz.mcr) <- mcr.names
            konz[rownames(konz.mcr), 1] <- konz.mcr
            dat.calc <- cbind(t(compkonz.mcr), dat.calc[, grep("MCR.Loading", colnames(dat.calc), invert = TRUE)]) 
         }        
      
      Konz.BVLS <- round(as.data.frame(t(konz[, !colnames(konz) %in% c("CompID", "CompName"), drop = FALSE])), 3)

   ## Output ----------------
      return(list(Konz.BVLS = Konz.BVLS, pH = pH))
}

DOSIM.BVLS_DIFF <- function(fit.dat, Abl, Gl, FitInterval, Weights, bl, bu, comp.med, Calc.pHDiff, pH.low, pH.high) {
   ## Fit (alle Proben) -------------
      #CompKonz.BVLS <- Do.BVLS(fit.dat, Abl, Gl, FitInterval, Weights, bl, bu)

      ## Zusammenstellen Komponenten -----------
         comp.list <- list()
         spec.offset <- as.data.frame(t(rep(0.00005, ncol(comp.med)))) ; rownames(spec.offset) <- "Intercept" ; colnames(spec.offset) <- colnames(comp.med)
         comp.list[[1]] <- as.data.frame(rbind(spec.offset, comp.med))
   
      ## Berechnung Fit ----------
         fit.bvls <- Fit.BVLS(fit.dat, comp.list, bas = NULL, Abl, Gl, FitInterval, Weights, bl, bu)

      ## Erstellen Fitergebnisse ---------
         dat.res <- as.data.frame(matrix(nrow = nrow(fit.dat), ncol = nrow(comp.med)))
         colnames(dat.res) <- rownames(comp.med) ; rownames(dat.res) <- rownames(fit.dat)
         dat.res$Intercept <- NA ; dat.res$SqSum <- NA
         #dat.res[, "SqSum"] <- fit.bvls$SqSum
   
         dat.res[, colnames(fit.bvls)] <- fit.bvls
         
         if (FALSE) {
            for (g in rownames(comp.med)) {
               if (!grepl("MCR.Loading", g)) {
                  dat.res[, g] <- dat.res[, g] * comp.med[g, "Konz_Komp_mg_mL"]
               }
            }
         }
      #CompKonz.BVLS <- dat.res
 
   ## Berechnung der Konzentrationen ---------
      dat.calc <- dat.res[, c(rownames(comp.med), "SqSum")]
      
      ID.dat <- c(levels(factor(comp.med$CompID)), "SqSum")
      dat.temp <- as.data.frame(t(dat.calc))
      dat.temp$CompID <- c(comp.med$CompID, "SqSum")
      dat.temp$CompName <- c(comp.med$CompName, "SqSum")
      konz <- SUMvers9(dat.temp, ID.dat, addChar = TRUE, File.ID = dat.temp$CompID, Counter = FALSE)$SUM
      rownames(konz) <- konz$CompName 
      pks <- c(3.13, 4.5, 5.8)
      
      ## Berechnung pH-Werte Monster Original und Monster pH-moduliert aus dem Diff-Spektrum der reinen Zitronensäure ------------
         if(is.na(pH.low)) pH.low <- pH.fix[1]
         if(is.na(pH.high)) pH.high <- pH.fix[2]
         
         if (Calc.pHDiff) {
               calccomp.ph <- "Zitronensaeure"
               sel <- rownames(comp.med)[grep(calccomp.ph , comp.med$CompName)]
               E <- as.data.frame(complist.svd[[calccomp.ph]], TRUE)
                                 
               if (!Gl == 0) {
                  comp.der <- SavGol(comp.med[sel, ], Gl, 2, Abl)
                  E <- SavGol(E, Gl, 2, Abl)
               } else {
                  comp.der <- comp.med[sel, ]
               }
                                                
               mcr.konz <- as.matrix(dat.res[, sel])
               spec.dat <- mcr.konz %*% as.matrix(comp.der[, FitInterval])                                
               C.cit <- spec.dat %*% pinv(E[, FitInterval])
   
               pH <- round(unlist(calc.pH(C.cit, pks)), 2)               
               if (any(is.na(pH))) {
                  pH <- c(pH.low, pH.high)
               }
         } else {
            pH <- c(pH.low, pH.high)
         }
                          
      ## Rückfit MCR-Loadings mit Komponentenspektren ----------------         
         mcr.num <- grep("MCR.Loading", rownames(comp.med))
         
         if (length(mcr.num) > 0) {
            mcr.names <- unique(comp.med[mcr.num, "CompName"])
            konz.mcr <- NULL
            compkonz.mcr <- NULL
            konz.mcr <- matrix(nrow = length(mcr.names), ncol = nrow(dat.res)) ; rownames(konz.mcr) <- mcr.names ; colnames(konz.mcr) <- rownames(dat.res)
                
            ## Berechnung Konzentrationswerte --------------                                                
               for (x in mcr.names) {
                  ## Reinspektren der Komponente --------
                     sdb.med <- as.data.frame(complist.svd[[x]], TRUE)
                     sel <- rownames(comp.med)[grep(x , comp.med$CompName)]
                     koef.mcr <- list()
                  
                     if (!Gl == 0) {
                        comp.der <- SavGol(comp.med[sel, ], Gl, 2, Abl)
                        sdb.med <- SavGol(sdb.med, Gl, 2, Abl)
                     } else {
                        comp.der <- comp.med[sel, ]
                     }
                     
                  ## Rückfit des pH-Diffspektrums pro Komponente in Monster aus MCR.Loadings -------
                     mcr.konz <- as.matrix(dat.res[, sel])
                     spec.dat <- mcr.konz %*% as.matrix(comp.der[, FitInterval])
                                    
                     mcr.fit <- lm(t(spec.dat) ~ t(sdb.med[, FitInterval]))                  
                     koef.mcr <- sapply(1:nrow(spec.dat), function(x) t(coef(mcr.fit))[x, ])
                  
                  ## Korrektur Verdünnung durch Zugabe von Lauge zu Monsterlösung ----------
                     koef.mcr <- t(t(koef.mcr[-1, ]) / as.vector(t(konz["Original", colnames(konz.mcr)])))
                  
                     if (any(x == c("Benzoesaeure", "Sorbinsaeure", "Zitronensaeure"))) {
                        ## Berücksichtigung relative Konzänderung innerhalb pH-Sprung -----------------
                           pks.x <- as.vector(t(unique(sdb.med[, grep("pks", colnames(sdb.med))])))
                           c.high <- component.fractions(pks.x, pH[2])     
                           c.low <- component.fractions(pks.x, pH[1])
                           c.diff <- apply(c.high - c.low, 1, function(x) sum(abs(x)) / 2)
                           
                           koef.mcr <- t(t(koef.mcr) / c.diff)                              
                           konz.mcr[x, ] <- apply(koef.mcr, 2, function(x) sum(abs(x)) / 2)
                           compkonz.x <- koef.mcr
                     } else {                               
                        konz.neg <- apply(koef.mcr, 2, function(x) sum(x[x < 0]))
                        konz.pos <- apply(koef.mcr, 2, function(x) sum(x[x > 0]))
                        konz.ab <- abs(konz.pos + konz.neg)
                        konz.diff <- abs(konz.neg)
                        konz.mcr[x, ] <- konz.diff + konz.ab
                        compkonz.x <- koef.mcr
                     }                  
                     
                  rownames(compkonz.x) <- rownames(sdb.med) ; colnames(compkonz.x) <- rownames(dat.res)
                  compkonz.mcr <- rbind(compkonz.mcr, compkonz.x)                     
               }
            
               ## Einfügen Konzentrationswerte aus Rückfit ------------                 
                  konz[rownames(konz.mcr), colnames(konz.mcr)] <- konz.mcr
                  dat.calc <- cbind(t(compkonz.mcr), dat.calc[, grep("MCR.Loading", colnames(dat.calc), invert = TRUE)]) 
         }        
      
      Konz.BVLS <- round(as.data.frame(t(konz[, !colnames(konz) %in% c("CompID", "CompName"), drop = FALSE])), 3)

   ## Output ----------------
      return(list(Konz.BVLS = Konz.BVLS, pH = pH))
}

Fit.BVLS <- function(fit.dat, comp.list, bas, Abl, Gl, FitInterval, Weights, bl, bu) {

   ## Datenvorbereitung ---------
      if (!Gl == 0) {
         comp.list[[1]][2:nrow(comp.list[[1]]), ] <- SavGol(comp.list[[1]][2:nrow(comp.list[[1]]), ], Gl, 2, Abl)
         fit.dat <- SavGol(fit.dat, Gl, 2, Abl)
         if(!is.null(bas)) {
            bas <- SavGol(bas, Gl, 2, Abl)
         }
      }
      
   ## Gewichten der Komponenten- und Fitdaten -----------
      comp.w <- apply(comp.list[[1]][, FitInterval], 1, function(x)  x * Weights)
      fit.w <- Weights * t(fit.dat[, FitInterval])
      #bas[, FitInterval] <- t(Weights * t(bas[, FitInterval]))

   ## BVLS-Fit --------
      fit.model.bvls <- bvls(comp.w, fit.w, bl, bu)

      #koef.bvls <- as.data.frame(t(coef(fit.model.bvls)))
      koef.bvls <- as.data.frame(t(fit.model.bvls$x))
      rownames(koef.bvls) <- rownames(fit.dat) ; colnames(koef.bvls) <- colnames(comp.w)

   ## Residuen ------------      
      if (FALSE) {
         specs.bvls <- matrix(ncol = length(FitInterval), nrow = ncol(fit.w))
         res.bvls <- matrix(ncol = length(FitInterval), nrow = ncol(fit.w))
         specs.bvls <- t(fitted.values(fit.model.bvls) / Weights)
         res.bvls <- t(residuals(fit.model.bvls) / Weights)
         fit.w <- fit.w / Weights      
         specs.bvls <- as.data.frame(specs.bvls) ; rownames(specs.bvls) <- rownames(fit.dat) ; colnames(specs.bvls) <- FitInterval
         res.bvls <- as.data.frame(res.bvls) ; rownames(res.bvls) <- rownames(fit.dat) ; colnames(res.bvls) <- FitInterval

         ## SquareSum Residuals -------
            if(!is.null(bas)) {
               sqsum.bas <- mean(apply(bas[, FitInterval]^2, 1, sum))
               sqsum.res <- apply(res.bvls^2, 1, sum)
               sqsum.res <- as.data.frame(sqsum.res / sqsum.bas)
            }
      
         ## Ausgabeliste ---------
            OUT <- list()
            OUT$Res <- res.bvls
            OUT$Koef <- koef.bvls
            if(is.null(bas)) {
               OUT$SqSum <- NULL
            } else {
               OUT$SqSum <- sqsum.res
            }
      }
      
      return(koef.bvls)
}

## Auswerteparameter nach Methode zuordnen -----------
   Select.MethodPar <- function(method) {
      nam <- paste("M", method, sep = "")
      
      MethodPar <- switch(nam,      
                  M101 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M101, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure")
                                 ),

                  M102 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M102, 
                                       DI = DI$G, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure")
                                 ),

                  M103 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M103, 
                                       DI = DI$GP, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure")
                                 ),

                  M104 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M104, 
                                       DI = DI$GSP, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure")
                                 ),

                  M105 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M105, 
                                       DI = DI$GP, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.05, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 0.2, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                               
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure")
                                 ),

                  M107 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M107, 
                                       DI = DI$GSP, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.1, Max = 0.8),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 4, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure")
                                 ),
                                 
                  M110 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M110, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Gesamtzucker   = list(Min = 70, Max = 110),
                                                               Zitronensaeure = list(Min = 0, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Glukose|Saccharose|Fruktose|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M120 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M120, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Gesamtzucker   = list(Min = 70, Max = 100),
                                                               Zitronensaeure = list(Min = 0, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Glukose|Saccharose|Fruktose|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M130 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M130, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Gesamtzucker   = list(Min = 15, Max = 35),
                                                               Zitronensaeure = list(Min = 0, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Glukose|Saccharose|Fruktose|Zitronensaeure|Gesamtzucker|Ethanol|Benzoesaeure")
                                 ),

                  M131 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M131, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Gesamtzucker   = list(Min = 10, Max = 30),
                                                               Zitronensaeure = list(Min = 0, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Glukose|Saccharose|Fruktose|Zitronensaeure|Gesamtzucker|Ethanol|Benzoesaeure")
                                 ),
                                    
                  M150 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M150, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Gesamtzucker   = list(Min = 50, Max = 90),
                                                               Zitronensaeure = list(Min = 0, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Glukose|Saccharose|Fruktose|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M160 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M160, 
                                       DI = DI$Basic, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Gesamtzucker   = list(Min = 70, Max = 110),
                                                               Zitronensaeure = list(Min = 0, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Glukose|Saccharose|Fruktose|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M170 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M170, 
                                       DI = DI$G, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0, Max = 0.8),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 1, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M180 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M180, 
                                       DI = DI$G, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0, Max = 0.8),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 1, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M190 =   list(       FI = FI, 
                                       CompAnalyt = CompAnalyt$M190, 
                                       DI = DI$G, 
                                       BL = BL, 
                                       Comp.Limits = list   (  Koffein        = list(Min = 0, Max = 0.8),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 1, Max = 10)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol")
                                 ),

                  M201 =   list(       FI = FI, 
                                       FI.Diff = FI.Diff,
                                       CompAnalyt = CompAnalyt$M101, 
                                       CompAnalyt.Diff = CompAnalyt$Diff,
                                       DI = DI$Basic, 
                                       DI.Diff = DI$Basic, 
                                       DI.Mod = DI.Mod, 
                                       BL = BL, 
                                       BL.Diff = BL.Diff,
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure",
                                                            IDX2 = "Sorbinsaeure|Benzoesaeure")
                                 ),

                  M202 =   list(       FI = FI, 
                                       FI.Diff = FI.Diff,
                                       CompAnalyt = CompAnalyt$M102, 
                                       CompAnalyt.Diff = CompAnalyt$Diff,
                                       DI = DI$G, 
                                       DI.Diff = DI$Basic, 
                                       DI.Mod = DI.Mod, 
                                       BL = BL, 
                                       BL.Diff = BL.Diff,
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure",
                                                            IDX2 = "Sorbinsaeure|Benzoesaeure")
                                 ),

                  M203 =   list(       FI = FI, 
                                       FI.Diff = FI.Diff,
                                       CompAnalyt = CompAnalyt$M103, 
                                       CompAnalyt.Diff = CompAnalyt$Diff,
                                       DI = DI$GP, 
                                       DI.Diff = DI$Basic, 
                                       DI.Mod = DI.Mod, 
                                       BL = BL, 
                                       BL.Diff = BL.Diff,
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure",
                                                            IDX2 = "Sorbinsaeure|Benzoesaeure")
                                 ),

                  M204 =   list(       FI = FI, 
                                       FI.Diff = FI.Diff,
                                       CompAnalyt = CompAnalyt$M104, 
                                       CompAnalyt.Diff = CompAnalyt$Diff,
                                       DI = DI$GSP, 
                                       DI.Diff = DI$Basic, 
                                       DI.Mod = DI.Mod, 
                                       BL = BL, 
                                       BL.Diff = BL.Diff,
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.2, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 3, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure",
                                                            IDX2 = "Sorbinsaeure|Benzoesaeure")
                                 ),

                  M205 =   list(       FI = FI, 
                                       FI.Diff = FI.Diff,
                                       CompAnalyt = CompAnalyt$M105, 
                                       CompAnalyt.Diff = CompAnalyt$Diff,
                                       DI = DI$GP, 
                                       DI.Diff = DI$Basic, 
                                       DI.Mod = DI.Mod, 
                                       BL = BL, 
                                       BL.Diff = BL.Diff,
                                       Comp.Limits = list   (  Koffein        = list(Min = 0.05, Max = 0.4),
                                                               Gesamtzucker   = list(Min = 90, Max = 130),
                                                               Taurin         = list(Min = 0.2, Max = 5),
                                                               Zitronensaeure = list(Min = 6, Max = 9)
                                                            ),                     
                                       Freigabe = list   (  IDX1 = "Koffein|Glukose|Saccharose|Fruktose|Taurin|Zitronensaeure|Gesamtzucker|Ethanol|Sorbinsaeure|Benzoesaeure",
                                                            IDX2 = "Sorbinsaeure|Benzoesaeure")
                                 )
      )
      
      return(MethodPar)
   }     

################################################################################
## Ergänzungsfuktionen zu mbioSpec 1.2
################################################################################

.mbioCodeToProperty <- function(x)
{
  x <- x[which(substr(x, 1, 9) == "mbioCode:")[1]];
  x <- substr(x, 10, 500);
  x <- paste("0x",substring(x,seq(1,5*16,2),seq(2,5*16,2)),sep="");
  x <- t(array(as.raw(as.integer(x)),c(8,5)));
  x <- unlist(lapply(1:5, function(y) readBin(x[y,],double())));
  names(x) = c("NPT", "FXV", "LXV", "RES", "DXU");
  return(x);
}

as.WaveNum <-  function (...)
{
    wn <- substr(WaveNums, 2, 5)
    get.WaveIndex  <- function (Value)
    {
        Value   <-  as.numeric(Value);
        idx <-  match(c(Value,Value + 1,Value - 1), wn)
        if (!is.na(idx[1])){
            return (idx[1])
        }
        else if (!is.na(idx[2])){
            cat("X", Value, " is shiftet by one to: X", Value + 1, "\n", sep=""); flush.console();
            return (idx[2])
        }
        else if (!is.na(idx[3])){
            cat("X", Value, " is shiftet by one to: X", Value - 1, "\n", sep=""); flush.console();
            return (idx[3])
        }
        else{
            stop(paste("WaveNr. X", Value, " not found !!\n", sep=""))
        }
    }
    parameters  <-  as.list(match.call())[c(-1)];
    WNrRegEx    <-  "^[ ]*[Xx][0123456789]{1,5}[ ]*$";
    WaveVector  <-  c()

    for (param in parameters){
        class(param)    <-  "character"
        is.WaveNr   <-  regexpr(WNrRegEx, param) > -1
        is.area <-  (length(param) == 3 && param[[1]] == ":" && all(is.WaveNr[2:3]))
        if (is.area){
            WaveVector  <-  c(WaveVector, paste("X", wn[get.WaveIndex(substr(param[2], 2, 500)):get.WaveIndex(substr(param[3], 2, 500))], sep = ""))
        }
        else if (is.WaveNr[1]){
            WaveVector  <-  c(WaveVector, paste("X", wn[get.WaveIndex(substr(param[1], 2, 500))], sep = ""))
        }
        else{
            if (param[1] == ":") is.WaveNr[1] <-  TRUE;
            stop(paste("WaveNr. ", param[!is.WaveNr], " not found !!\n", sep = ""))
        }
    }
    return (WaveVector)
}


################################################################################
## Excel-Files laden
################################################################################

Read.XLS <- function(Dir.XLS, xlsfile, method = c("csv", "tab"), sheet = 1, pattern, pattern.header, SkipRows = 1, SkipHeaderRows = 0, Dir.Perlxls = "C:/Programs/Perl/perl/site/lib", Dir.Perlexe = "C:/Programs/Perl/perl/bin") {

   ## Beschreibung Variablen -----------
      #Dir.XLS          vollständiger Pfad, in dem Excel-File vorliegt
      #xlsfile:         Filename des Excel-Files
      #method:          Auswahl, ob Konvertierung des Excel-Files in csv- oder tab-Format durchgeführt werden soll
      #sheet:           Auswahl des einzulesenden Tabellenblattes in Excel-File
      #pattern:         wenn ausgewählt, liest Zeilen in Excel-File erst ab angegebenem Pattern ein
      #SkipRows:        überspringt Zeilen in Excel-Tabelle und liest Daten erst ab nachfolgender Zeilennummer ein (derzeit nur für csv-Format). Wird von pattern - falls angeben - überschrieben
      #SkipHeaderRows:  überspringt Zeilen in Excel-Tabelle und liest Header erst ab nachfolgender Zeilennummer ein (derzeit nur für csv-Format)
      #Dir.Perlxls:     Verzeichnis, in dem Perlprogramme "xls2csv.pl" oder "xls2tab.pl" zu finden sind
      
   ## interne Funktionen ------
      ## Einlesen Ascii -------
         dQuote.ascii <- function(x)
         {
             paste('"', x,'"', sep = '')
         }

      ## Einlesen csv-Datei mit Header -----------      
         Read.csv <- function(con, ...) {
            param <- read.csv(con, header = FALSE, skip = SkipRows, sep = ",", dec = ".", stringsAsFactors = FALSE)
            param.header <- read.csv(con, header = FALSE, skip = SkipHeaderRows, nrows = 1, sep = ",", stringsAsFactors = FALSE)
            colnames(param) <- param.header
            
            return(param)
         }
      
   ## Konvertieren xls-File ---------
      xls <- paste(Dir.XLS, xlsfile, sep = "")
      if (file.access(xls, 4) != 0) {
        #stop("Unable to read xls file '", xls, "'.")
      }
      if (method == "csv") {
        script <- file.path(Dir.Perlxls, "xls2csv.pl")
        con <- paste(tempfile(), "csv", sep = ".")
      } else {
         if (method == "tab") {
            script <- file.path(Dir.Perlxls, "xls2tab.pl")
            con <- paste(tempfile(), "tab", sep = ".")
         } else {
            stop("Unknown method", method)
         }
      }
      cat("Converting xls file to", method, "file ...\n")
      #cmd <- paste("perl", script, dQuote.ascii(xls), dQuote.ascii(con), sheet, sep = " ")
      Dir.Prog <- paste(Dir.Perlexe, "/perl.exe", sep = "")      
      cmd <- paste(Dir.Prog, script, dQuote.ascii(xls), dQuote.ascii(con), sheet, sep = " ")
      results <- system(cmd)
            
   ## Einlesen konvertiertes File ------------
      if (!missing(pattern)) {
         cat("Searching for lines containing pattern ", pattern, "... ")
         idx <- grep(pattern, readLines(con))
         SkipRows <- idx[1] - 1
         if (length(idx) == 0) {
            warning("\npattern not found")
            return(NULL)
         } else {
            cat("Done.\n")
         }
      }
      if (!missing(pattern.header)) {
         cat("Searching for lines containing header pattern ", pattern.header, "... ")
         idx <- grep(pattern.header, readLines(con))
         SkipHeaderRows <- idx[1] - 1
         if (length(idx) == 0) {
            warning("\npattern not found")
            return(NULL)
         } else {
            cat("Done.\n")
         }
      }         
      cat("Reading", method, "file ... ")
      if (method == "csv") {
         param <- Read.csv(con)
      } else {
         if (method == "tab") {
            param <- read.delim(con)
         } else {
            stop("Unknown method", method)
         }
      }
      cat("Done.\n")
      
   ## Löschen temporäres File ---------
      file.remove(con)
      cat("Deleting temporary", method, "file ... ")
      cat("Done.\n")
                   
   ## Output ------
      return(param)
}


###################################################################################################
## Schichtdicke FFT
###################################################################################################

##
## periodtothick : calculates the layer thickness corresponding to the
##               period of the sine function which describes the baseline
##
## example :  periodtothick(100, 4, 1.35)
##               the thickness (9.26um) is calculated, which corresponds
##               to the period of 100 data points collected at the resolution
##               of 4 cm^-1, refractive index being 1.35
##
periodtothick <- function(periodnpoints, resolution, refraction)
{
  return( 10000/2/periodnpoints/refraction/resolution )
}  
  
##
## thicktoperiod : calculates the  period of the sine function (the number of
##               data points) which would describe the baseline emerging when
##               working with a sample with the given layer thickness.
##
## example :  thicktoperiod(9.26, 4, 1.35)
##               the period (100 data points) is calculated, which corresponds
##               to the thickness of 9.26um when acquiring the spectrum with
##               the resolution of 4 cm^-1, refractive index being 1.35
##
thicktoperiod <- function(thickness, resolution, refraction)
{
  return( 10000/2/thickness/refraction/resolution )
}
    
##
## freqspect : calculates the Fourier-spectrum of a data set
##
## returns a list with the components
##  $periods.npoints - a range of periods (expressed as numbers of data points),
##  $fourier.weights - the weights of the corresponding Fourier-frequencies. 
## 
freqspect <- function(Y)
{ 
  N <- length(Y)
  X1 <- 2/seq(0, 1, length.out = floor(N/2) )
  Y1 <- fft(Y, inverse = FALSE)
  Y1 <- abs(Y1[1:floor(N/2)])
  return(list(fourier.weights = Y1, periods.npoints = X1))
}

## 
##  find_thick : recovers the layer thickness for the given IR-spectrum
##
##  syntax : find_thick(Spectr, Wavenumbers, find_thick.par = NULL)                                   
##           
##  where
##               Spectr - a vector containing the spectral data;
##               Wavenumbers - a vector containing the wavenumbers;
##               find_thick.par - a list of optional parameters 
##
## returns a list with the components:
##  $thick_found - the calculated layer thickness [um],
##  $SNR - the signal-to-noise ratio for the peak of interest in the Fourier-spectrum. 
##  $SNR2 - the ratio peak of interest to the second tallest peak
##  $X and $Y - the x- and y-coordinates of the Fourier-spectrum.
##
## returns only the value of the thickness if find_thick.par$ListOutput = FALSE
##
## examples : 
##    find_thick(Spectr, Wavenumbers)
##    find_thick(Spectr, Wavenumbers, list(NumberOfZeros=1e5, ListOutput=FALSE))
##
find_thick <- function(Spectr, Wavenumbers, find_thick.par = NULL)
{ 
   ## Parsing input ------------
      if (!is.list(find_thick.par)) find_thick.par <- list()
      if (is.null(find_thick.par$min_thick)) find_thick.par$min_thick <- 5 # minimum thickness
      if (is.null(find_thick.par$max_thick)) find_thick.par$max_thick <- 20 # max. thickness
      if (is.null(find_thick.par$absmin_thick)) find_thick.par$absmin_thick <- 4 # lower thickness values are disregarded when predicting SNR
      if (is.null(find_thick.par$absmax_thick)) find_thick.par$absmax_thick <- 150 # higher thickness values are disregarded when predicting SNR
      if (is.null(find_thick.par$Refrcorrect)) find_thick.par$Refrcorrect <- TRUE # whether to correct for non-constant refractive ind.
      if (is.null(find_thick.par$BigValOff)) find_thick.par$BigValOff <- FALSE # whether to set the intensive signals to 0
      if (is.null(find_thick.par$Calc.SNR)) find_thick.par$Calc.SNR <- TRUE # whether to calculate SNR
      if (is.null(find_thick.par$Calc.MMR)) find_thick.par$Calc.MMR <- TRUE # whether to calculate the ratio of maxima
      if (is.null(find_thick.par$Apodisation)) find_thick.par$Apodisation <- TRUE # whether to apply Hamming's apodisation
      if (is.null(find_thick.par$NumberOfZeros)) find_thick.par$NumberOfZeros <- 1e5 # number of zeros for zero-filling
      if (is.null(find_thick.par$ListOutput)) find_thick.par$ListOutput <- TRUE # output as a list
   
   ## Calculating the Fourier spectrum and determining the thickness ----------
   if (find_thick.par$BigValOff) Spectr[abs(Spectr) > quantile(abs(Spectr), .9)] <- 0
      #if (find_thick.par$Refrcorrect) Spectr <- correct_refr(Spectr, Wavenumbers) 
   if (find_thick.par$Refrcorrect) {
      Spectr <- correct_refr(Spectr, Wavenumbers, ListOutput = TRUE)
      Wavenumbers <- Spectr$x
      Spectr <- Spectr$y
   } 
   resolution <- Wavenumbers[1] - Wavenumbers[2];
   refraction <- n_H2O(Wavenumbers[1]) # The refractive index for the 1st wavenumber is used for the calculations  
   N <- length(Spectr)
   if (find_thick.par$Apodisation) { 
      Win_coeff <- 1 - cos( 2*pi*(0:(N-1))/(N-1) )
      Spectr <- Spectr * Win_coeff 
   }
   NZ <- max( round(log2(N + find_thick.par$NumberOfZeros)), ceiling(log2(N)) )
   NZ <- (2^NZ - N) / 2    
   XY <- freqspect(c(rep(0, floor(NZ)), Spectr, rep(0, ceiling(NZ))))  # Two-sided zero-filling
   X <- periodtothick(XY$periods.npoints, resolution, refraction)
   Y <- XY$fourier.weights 
   Pointsofinterest <- which( (X > find_thick.par$min_thick ) & (X < find_thick.par$max_thick ) )
   if (length(Pointsofinterest)==0) Pointsofinterest <- 1:length(X)
   thick0 <- findallmax( Y[Pointsofinterest])  
   MMR <- thick0$x[1] / thick0$x[2]
   thick0 <- Pointsofinterest[thick0$ix[1]]
   thick_found <- sum( X[(thick0-1):(thick0+1)] * Y[(thick0-1):(thick0+1)] ) /
      sum( Y[ (thick0-1) : (thick0+1) ] )
      
   ## calculating S/N-ratio ------------
      if (find_thick.par$Calc.SNR & find_thick.par$ListOutput) {
         Noise <- Y[- c( 
            which((X > find_thick.par$min_thick) & (X < find_thick.par$max_thick)), 
            which(X < find_thick.par$absmin_thick) , which(X > find_thick.par$absmax_thick) 
            ) ]
         if (is.null(Noise)) Noise <- Y[- c( which((X > thick_found-1) & (X < thick_found+1)), which(X < find_thick.par$absmin_thick) )] 
         SNR <- Y[thick0] / mean(Noise^2)^.5
         attributes(SNR)$names <- NULL
      } else {
         SNR <- NA
      }

   ## ratio of maxima ------------      
      if (!find_thick.par$Calc.MMR ) MMR <- NA

   Output <- list(thick_found = thick_found, SNR = SNR, MMR = MMR, X = X, Y = Y)
   if (find_thick.par$ListOutput) return(Output) else return(thick_found)
}

##
## findallmax: determining the 1st and 2nd biggest local maxima.
## 
## Used by find_thick
##
findallmax <- function(Y)
{   
   n <- 2
   maxY <- which.max(Y)
   max.pos <- c(cummax( Y[1 : (maxY + n - 1)]), cummax(Y[length(Y) : maxY]))
   max.pos <- table(max.pos)
   max.pos <- as.numeric(names(max.pos[max.pos >= n]))
   max.pos <- sort(max.pos, decreasing = TRUE)
   allmax <- list(x = max.pos, ix = maxY  )
   return(allmax)      
}

##
## findallmax_: looking for all the local maxima on the given segment
##              (can be used in place of findallmax).
## 
findallmax_ <- function(Y)
{   
   x0 <- diff(c(0, Y))
   x01 <- diff(c(Y, 0))
   max.pos <- which((x0>=0) & (x01<0))
   allmax <- sort(Y[max.pos], decreasing = TRUE, index.return = TRUE)
   allmax$ix <- max.pos[allmax$ix]
   return(allmax)      
}


##
## correct_refr: calculates the corrected spectrum (baseline) when acquiring
##               it experimentally within a broad range of wavenumbers so that
##               the refractive indices of water are different. 
##
## correct_refr(Spectr, Wavenumbers = NULL, maxrefr = NULL, minrefr = NULL, ListOutput = FALSE)
##
## The optional parameters maxref, minref are used when only the maximum and
## minimum refractive indices are given, assuming they exhibit a linear trend.
## Otherwise the wavenumbers are used, the refractive indices being interpolated
## from the reference data:
## http://refractiveindex.info/?group=LIQUIDS&material=Water
##
## If ListOutput is set to FALSE, the function returns the corrected counts corresponding
## to the initial wavenumbers, otherwise a list is calculated, its components $x and $y being
## certain "new" evenly spaced wavenumbers and the corresponding signal values.
##
correct_refr <- function(Spectr, Wavenumbers = NULL, maxrefr = NULL, minrefr = NULL, ListOutput = FALSE)
{
   N <- length(Spectr)
   if (is.null(Wavenumbers)) {
      Wavenumbers <- N:1
      if (is.null(maxrefr)) maxrefr <- 1.375
      if (is.null(minrefr)) minrefr <- 1.2
      Refraction <- seq(maxrefr, minrefr, len=N)
   } else {
      Refraction <- n_H2O(Wavenumbers)
   }
   
   Wave.seem <- Wavenumbers * Refraction / Refraction[1]
   if (!ListOutput) {
      out <- spline(Wave.seem, Spectr, N, "fmm", xout = Wavenumbers)$y
   } else {
      out <- spline(Wave.seem, Spectr, N, "fmm")
      out$x = rev(out$x)
      out$y = rev(out$y)
   }        
   return(out)
   # The corrected baseline is calculated for evenly spaced wavenumbers!!!
} 

##
## n_CaF2: calculates the rerfractive index of CaF2 for the given
##         wavenumber [cm^-1].
##
## see:
## http://refractiveindex.info/?group=CRYSTALS&material=CaF2
##
n_CaF2 <- function(wavenumber)
{
   L2 <- (10000/wavenumber)^2
   n <- sqrt(1 + 0.5675888*L2/(L2-0.050263605^2) + 0.4710914*L2/(L2-0.1003909^2) 
      + 3.8484723*L2/(L2-34.649040^2))
   return(n)
}

##
## n_H2O: calculates the rerfractive index of H2O for the given
##         wavenumber [cm^-1].
##
## see:
## http://refractiveindex.info/?group=LIQUIDS&material=Water
##
n_H2O <- function(wavenumber)
{
   Wasser <- matrix(c(
   1.0,	1.327,	2.89E-6,
   1.2,	1.324,	9.89E-6,
   1.4,	1.321,	1.38E-4,
   1.6,	1.317,	8.55E-5,
   1.8,	1.312,	1.15E-4,
   2.0,	1.306,	1.10E-3,
   2.2,	1.296,	2.89E-4,
   2.4,	1.279,	9.56E-4,
   2.6,	1.242,	3.17E-3,
   2.65,	1.219,	6.70E-3,
   2.70,	1.188,	0.019,
   2.75,	1.15,	   0.059,
   2.80,	1.142,	0.115,
   2.85,	1.149,	0.185,
   2.90,	1.201,	0.268,
   2.95,	1.292,	0.298,
   3.00,	1.371,	0.272,
   3.05,	1.426,	0.240,
   3.10,	1.467,	0.192,
   3.15,	1.483,	0.135,
   3.20,	1.478,	0.0924,
   3.25,	1.467,	0.0610,
   3.30,	1.450,	0.0368,
   3.35,	1.432,	0.0261 ,
   3.40,	1.420,	0.0195 ,
   3.45,	1.410,	0.0132,
   3.50,	1.400,	0.0094,
   3.6,	1.385,	0.00515,
   3.7,	1.374,	0.00360,
   3.8,	1.364,	0.00340,
   3.9,	1.357,	0.00380,
   4.0,	1.351,	0.00460,
   4.1,	1.346,	0.00562,
   4.2,	1.342,	0.00688,
   4.3,	1.338,	0.00845,
   4.4,	1.334,	0.0103 ,
   4.5,	1.332,	0.0134 ,
   4.6,	1.330,	0.0147 ,
   4.7,	1.330,	0.0157 ,
   4.8,	1.330,	0.0150 ,
   4.9,	1.328,	0.0137 ,
   5.0,	1.325,	0.0124 ,
   5.1,	1.322,	0.0111 ,
   5.2,	1.317,	0.0101 ,
   5.3,	1.312,	0.0098,
   5.4,	1.305,	0.0103 ,
   5.5,	1.298,	0.0116 ,
   5.6,	1.289,	0.0142,
   5.7,	1.277,	0.0203 ,
   5.8,	1.262,	0.0330 ,
   5.9,	1.248,	0.0622 ,
   6.0,	1.265,	0.107  ,
   6.1,	1.319,	0.131  ,
   6.2,	1.363,	0.0880 ,
   6.3,	1.357,	0.0570 ,
   6.4,	1.347,	0.0449 ,
   6.5,	1.339,	0.0392 ,
   6.6,	1.334,	0.0356 ,
   6.7,	1.329,	0.0337 ,
   6.8,	1.324,	0.0327 ,
   6.9,	1.321,	0.0322 ,
   7.0,	1.317,	0.0320 ,
   7.1,	1.314,	0.0320 ,
   7.2,	1.312,	0.0321 ,
   7.3,	1.309,	0.0322 ,
   7.4,	1.307,	0.0324 ,
   7.5,	1.304,	0.0326 ,
   7.6,	1.302,	0.0328 ,
   7.7,	1.299,   0.0331 ,
   7.8,	1.297,	0.0335 ,
   7.9,	1.294,	0.0339 ,
   8.0,	1.291,	0.0343 ,
   8.2,	1.286,	0.0351 ,
   8.4,	1.281,	0.0361 ,
   8.6,	1.275,	0.0372 ,
   8.8,	1.269,	0.0385 ,
   9.0,	1.262,	0.0399 ,
   9.2,	1.255,	0.0415 ,
   9.4,	1.247,	0.0433 ,
   9.6,	1.239,	0.0454 ,
   9.8,	1.229,	0.0479 ), 
      ncol = 3, byrow = TRUE )

   Wasserwave <- 10000/Wasser[,1]  # cm^-1
   Wasser_n <- Wasser[,2]
   return(spline(Wasserwave, Wasser_n, xout = wavenumber)$y)
}

##########################################################################################################
## Berechnung relative Anteile verschiedener pH-Formen in Substanzen mit einem oder mehreren pKs-Werten
##########################################################################################################

## Anteile von verschiedenen protonierten Formen ----------------
##
## component.fractions(pK, pH)
##
## pK - pKs der Sauere
## pH - pH-Wert(e)
## Berechnet werden (als ein Vektor) die auf 1 normierten Anteile aller (de)protonierten Formen
## dieser Saeure bei dem angegebenen pH

component.fractions <- function(pK, pH) {
   pK <- cumsum( rev(pK) ) # pKs -> p(Beta_H) 
   m <- length(pK)
   pK <- c(0, pK)
   Cfr <- NULL
   for (i.pH in pH) {
    pH.j <- -c(0:m) * i.pH
    numerator <- 10^(pH.j + pK)
    denominator <- sum( numerator )
    Cfr <- rbind( Cfr, numerator / denominator )
   }
   colnames( Cfr ) <- paste("H", as.character(c(0:m)), sep = "")
   rownames( Cfr ) <- paste("pH", as.character(pH), sep = "_" )
   return( Cfr )
}

##########################################################################################################
## Pseudoinversion 
##########################################################################################################
 
pinv <- function(X, r = NULL) {
   if (is.null(r)) {
      s <- svd(X) 
   } else {
      s <- svd(X, r, r)
      s$d <- s$d[1:r]
   }
   
   if (length(s$d) > 1) {
      D <- diag(s$d)
   } else {
      D <- s$d
   }
   
   dat <- s$v %*% solve(D) %*% t(s$u)                                      ##  V D^-1 U'
   
   return(dat) 
}
  

##############################################################################################################
## pH-Suche 
##############################################################################################################
## calc.pH( delta.C , pKs = c(3.05, 4.5, 5.88) )
## find.pH( delta.C , pKs = c(3.05, 4.5, 5.88) )
## 
## delta.C - Vektor von 4 Konzentrationsdifferenzen für Zitrat (Cit3-, HCit2-, H2Cit-, H3Cit)
##           bei der pH-Änderung pH0 -> pH1;
## pKs - pKs-Werte Zitrat.
##
## Es wird eine Liste mit den Elementen (je 1 Zahl) ph0, ph1 berechnet.
## find.pH liest alle möglichen pH-Kombinationen aus, 
## calc.pH berechnet es für Zitrat.
             
calc.pH <- function( delta.C , pKs = c(3.05, 4.5, 5.88) ) {
  if ( any( delta.C < 0 ) ) {
    if ( delta.C[ 1 ] == 0 | delta.C[ 2 ] == 0 ) return ( NULL )
    if ( delta.C[ 2 ] > 0 ) {
      ph1 <- log10( delta.C[ 1 ] ) - log10( delta.C[ 2 ] ) + pKs[ 3 ]
      C.predict <- component.fractions( pKs , ph1 )
      conc <- delta.C[ 1 ] / C.predict[ 1 ]
      C.predict <- C.predict * conc 
      C.predict <- C.predict - delta.C
      ph0 <- pKs[ 1 ] - log10( C.predict[4] ) + log10( C.predict[3] )
      C.predict <- component.fractions( pKs , ph0 )
      C.predict <- C.predict * conc 
      C.predict <- C.predict + delta.C
      ph1 <- log10( C.predict[ 1 ] ) - log10( C.predict[ 2 ] ) + pKs[ 3 ]  
    } else {
      ph0 <- log10( -delta.C[ 3 ] ) - log10( -delta.C[ 4 ] ) + pKs[ 1 ]
      C.predict <- component.fractions( pKs , ph0 )
      conc <- -( delta.C[ 4 ] + delta.C[3] ) / ( C.predict[ 4 ] + C.predict[ 3 ] )
      C.predict <- C.predict * conc + delta.C
      ph1 <- log10( C.predict[ 1 ] ) - log10( C.predict[ 2 ] ) + pKs[ 3 ]
    }
    return( list( ph0=ph0, ph1=ph1 ) )
  } else { ## Erweiterung 23.08.2013
    pKs <- rev( pKs )
    x <- sort( delta.C, decreasing = TRUE, index.return = TRUE )
    ind.pK <- floor( mean( x$ix[ 1:2 ] ) )
    pH <- pKs[ ind.pK ] + log10( x$x[1] / x$x[2]  ) * (-1) ^ as.numeric( x$ix[1] > x$ix[2] )
    return( pH )
  }
}