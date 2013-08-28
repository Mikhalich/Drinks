##############################################################################################
## Modell DIS zur Bestimmung von Monster Energydrinks auf Quantos0001 mit konservierungsmittel
##############################################################################################

## Daten einladen --------------
   Dir <- "C:/"
   if (Dir == "C:/") Dir <- paste(Dir, "micro-biolytics/", sep = "") 
   Dir.Main <- paste(Dir, "Firmen und Projekte/Q-Food/QFood002C - Monster Energydrink/Modell DIS/", sep = "")
   Dir.Diff <- paste(Dir, "Firmen und Projekte/Q-Food/QFood002D - Konservierungsmittel/Modell DIS/", sep = "")
   Dir.Div <- paste(Dir, "Firmen und Projekte/Q-Food/QFood002E - diverse Softdrinks/Modell DIS/", sep = "")
   Dir.Model <- paste(Dir.Main, "R-Daten/", sep = "")
   Dir.ModelDiff <- paste(Dir.Diff, "R-Daten/", sep = "")
   Dir.SDB <- paste(Dir, "Spektrendatenbank Quantos0001/Datenbank/", sep = "")
   Dir.AtmoKorr <- paste(Dir, "Spektrendatenbank Quantos0001/Atmosphärische Kompensation/R-Daten/", sep = "")
         
   ## Funktionen laden ----------
      library(mBioSpec)
      #library(bvls)

      #dyn.load("C:/micro-biolytics/Firmen und Projekte/Q-Food/QFood2C - Monster Energydrink/Modell DIS/bvls.dll")              
      dyn.load(paste(Dir, "Firmen und Projekte/Q-Food/QFood002C - Monster Energydrink/Modell DIS/mBioerx.dll", sep = ""))              
      
      bvls <- function (A, b, bl, bu, key = 0, istate = rep(0, ncol(A) + 1)) {       
         M <- nrow(A)
         N <- ncol(A)
         X <- rep(0, N)
         if (!is.numeric(A)) {
           cat(" 'A' must be numeric \n")
           return()
         }
         if (!is.numeric(b)) {
           cat(" 'b' must be numeric \n")
           return()
         }
         if (length(bu) != N || length(bl) != N) {
           cat("bounds not correct length \n")
           return()
         }
         if (!(is.numeric(bu) && is.numeric(bl))) {
           cat("bounds contain non-numeric values \n")
           return()
         }
         W <- rep(0, N)
         mm <- min(M, N)
         act <- rep(0, M * (mm + 2))
         zz <- rep(0, M)
         loopA <- 0
         sol <- .Fortran("bvls", key = key, m = M, n = N, a = A, b = b,
           bl = bl, bu = bu, x = X, w = W, act = act, zz = zz, istate = istate,
           loopA = loopA, PACKAGE = "mBioerx", NAOK = TRUE)
         fitted <- A %*% sol$x
         resid <- b - fitted
         bvls.out <- list(x = sol$x, deviance = sum(resid^2), residuals = resid,
           fitted = fitted)
         class(bvls.out) <- "bvls"
         bvls.out
      }
      
      #source(paste(Dir, 'micro-biolytics/R-Workspace/Funktionen/MINmBioDataStruc.R', sep = ""))
      source(paste(Dir, 'R-Workspace/Funktionen/MWvers9.R', sep = ""))
      source(paste(Dir, 'R-Workspace/Funktionen/mBio_MathFunctions_V2.1.R', sep = ""))
      source(paste(Dir, 'R-Workspace/Funktionen/QFood/App Monster DIS/App Monster KM DIS - Funktionen4.r', sep = ""))
      source(paste(Dir, 'R-Workspace/Funktionen/mBio_Korrekturen mbioMVE_1.5.r', sep = ""))

   ## Komponentenspektren SDB laden ----------
      setwd(Dir.SDB)
      load(file = "130817 - Spektrendatenbank Quantos0001 Monster.rList")
      spec.case <- "AtmoKorr_NormMass"
      res.sel <- "RES4"
      pH.vers <- "pH_3.15_10.5"
      
      complist.sdb <- list()
      complist.svd <- list()
      
      ## Ableiten der CompSpectra nach SDB-Vorgaben ------------
         sel <- c("Glukose", "Fruktose", "Saccharose", "Zitronensaeure", "Phosphorsaeure")
         
         for (i in sel) {
            complist.sdb[[i]] <- savgol(sdblist.mw$CompSpectra[[spec.case]][[res.sel]][[i]][, c(3050:950),, FALSE], 1, 5, 2)
            complist.sdb[[i]] <- complist.sdb[[i]][, c(3020:2450, 2250:980),, FALSE]
         }
         
      complist.sdb[["Zitronensaeure"]] <- complist.sdb[["Zitronensaeure"]][grep("130410", complist.sdb[["Zitronensaeure"]]@PARAM$Datum_Teilansatz), ]
      
      sel.phos <- rownames(complist.sdb[["Phosphorsaeure"]])[(complist.sdb[["Phosphorsaeure"]]@PARAM$pH >= 3.15) & (complist.sdb[["Phosphorsaeure"]]@PARAM$pH <= 10.5)]
      complist.sdb[["Phosphorsaeure"]] <- complist.sdb[["Phosphorsaeure"]][sel.phos, ]

      complist.svd[["Phosphorsaeure"]] <- sdblist.mw$SingleComp[[spec.case]][[res.sel]][["Phosphorsaeure"]][[pH.vers]][["130806"]]               
      complist.svd[["Zitronensaeure"]] <- sdblist.mw$SingleComp[[spec.case]][[res.sel]][["Zitronensaeure"]][[pH.vers]][["130806"]]      
      complist.svd[["Benzoesaeure"]] <- sdblist.mw$SingleComp[[spec.case]][[res.sel]][["Benzoesaeure"]][[pH.vers]][["130807"]]      
      complist.svd[["Sorbinsaeure"]] <- sdblist.mw$SingleComp[[spec.case]][[res.sel]][["Sorbinsaeure"]][[pH.vers]][["130807"]]
      complist.svd[["Kohlensaeure"]] <- sdblist.mw$SingleComp[[spec.case]][[res.sel]][["Kohlensaeure"]][[pH.vers]][["130805"]]      
      
      cit.svd <- sdblist.mw$SingleComp[[spec.case]][[res.sel]][["Zitronensaeure"]][[pH.vers]][["130806"]]      
      
      rm(sdblist.mw)
              
   ## Spektrendatenbank laden --------------------
      setwd(Dir.SDB)
            
      ## Diff-Spektren ----------
         load("Qfood002D - Spektrendatenbank DIS Monster Konservierungsmittel Quantos0001 MCR_PCA.rData")
         comp.diff <- comp.dat
         rm(comp.dat)

      ## AB-Spektren ----------
         #load("Qfood2C - Spektrendatenbank DIS Monster Quantos0001 MCR_PCA ohne Monster Ingredients.rData")
         #load("Qfood002C - Spektrendatenbank DIS Monster Quantos0001 MCR_PCA Cit_MCR.rData")
         load("Qfood002C - Spektrendatenbank DIS Monster Quantos0001 MCR_PCA.rData")
                        
   ## Laden Wasserdamp- und CO2-Korrekturspektren ------------
      setwd(Dir.AtmoKorr)
      load("130801 - Korrekturspektren atmosphärische Kompensation Quantos0001.rList")
      AtmoKorrSpecs <- AtmoKorrSpecs$WD$PCA[["RES4"]][["Abl0"]]

   ## Laden Vorlagespektrum mbio.data ------------------
      setwd(Dir.Model)
      load("Vorlage mbiodata Spektrum.rData")

   ## Lade Tabellen 
      ## relevante Komponenten pro Analyt-Fit ------------
         #comp.analyt <- Read.XLS(Dir.Main, "Monster KM DIS -  relevante Komponenten pro Analyt-Fit.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt <- list()             
         CompAnalyt$M101 <- Read.XLS(Dir.Main, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 101.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M102 <- Read.XLS(Dir.Main, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 102.xlsx", sheet = 1, method = "csv", SkipRows = 1)         
         CompAnalyt$M103 <- Read.XLS(Dir.Main, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 103.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M104 <- Read.XLS(Dir.Main, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 104.xlsx", sheet = 1, method = "csv", SkipRows = 1)         
         CompAnalyt$M105 <- Read.XLS(Dir.Main, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 105.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M107 <- Read.XLS(Dir.Main, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 107.xlsx", sheet = 1, method = "csv", SkipRows = 1)         
         CompAnalyt$M110 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 110.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M120 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 120.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M130 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 130.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M131 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 131.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M150 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 150.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M160 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 160.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M170 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 170.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M180 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 180.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         CompAnalyt$M190 <- Read.XLS(Dir.Div, "Monster DIS - relevante Komponenten pro Analyt-Fit Methode 190.xlsx", sheet = 1, method = "csv", SkipRows = 1)    
         CompAnalyt$Diff <- Read.XLS(Dir.Diff, "Monster KM DIS -  relevante pH-aktive Komponenten pro Analyt-Fit.xlsx", sheet = 1, method = "csv", SkipRows = 1)

         ## Korrektur Leerzeichen ---------------
            for (x in names(CompAnalyt)) {
               for (i in 1:ncol(CompAnalyt[[x]])) {
                  CompAnalyt[[x]][, i] <- gsub("x ", "x", CompAnalyt[[x]][, i])
               }
               
               colnames(CompAnalyt[[x]]) <- gsub(" ", "", colnames(CompAnalyt[[x]]))
            }
                                    
      ## Justierfaktor und LOD pro Analyt ------------
         lod.analyt <- Read.XLS(Dir.Main, "Justierfaktor und LOD relevanter Analyten Monster KM DIS.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         lod.analyt.diff <- Read.XLS(Dir.Diff, "Justierfaktor und LOD relevanter pH-aktiver Analyten Monster DIS.xlsx", sheet = 1, method = "csv", SkipRows = 1)
         
      ## Lade AnalytID Tabelle (CAS Nummern) ------------
         cas <- Read.XLS(paste(Dir, "mbio Produkte/DOSIM/", sep = ""), "AnalytID für DOSIM-Apps.xls", sheet = 1, method = "csv", SkipRows = 1)      
         
   ## Lade RF Importancewerte für SDB-Komponenten in Monster Energydrink ----------
      setwd(Dir.Model)
      load(paste("RF Importancewerte für SDB-Komponenten in Monster Energydrink DIS RES4 Manipulation .rList"))
      DI <- Dat.Imp

      setwd(Dir.ModelDiff)
      load(paste("RF Importancewerte für pH-aktive DiffSDB-Komponenten in Monster Energydrink DIS RES4 AtmoKorr_NormMass.rList"))
      DI.Diff <- Diff.Imp
      
      DI.Mod <- DI$Basic
      for (i in names(DI.Mod)) {
         for (x in names(DI.Mod[[i]])) {
            DI.Mod[[i]][[x]] <- DI$Basic[[i]][[x]]
            DI.Mod[[i]][[x]][,] <- 1
         }
      } 
       
      #load("RF Importancewerte für SDB-Komponenten in Monster Energydrink DIS RES4 Manipulation ohne Monster Ingredients.rData") 
      #DI.oMI <- Dat.Imp                           
      #load("RF Importancewerte für SDB-Komponenten in Monster Energydrink DIS RES4 Manipulation mit Monster Ingredients.rData") 
      #DI.mMI <- Dat.Imp   

   ## Lade Constraints für SDB-Komponenten in Monster Energydrink ----------
      #setwd(Dir.Model)
      #load("Uservorgabe BVLS-Constraints Cmin Monster All_Comp3.rData") 
      #BL <- Cmin.Range
      
      BL <- as.data.frame(matrix(0, 1, nrow(comp.dat) + 1)) ; colnames(BL) <- c("Intercept", rownames(comp.dat))
      BL[, grep("Glukose", colnames(BL))] <- -1
      BL[, grep("Intercept", colnames(BL))] <- -1
      BL[, grep("130405-07_Sorbinsaeure_pH7,1", colnames(BL))] <- -1
      BL[, grep("MFET", colnames(BL))] <- -1
      BL[, grep("MonsterFlavour", colnames(BL))] <- -1
      BL[, grep("Wasserspektrum", colnames(BL))] <- -150
      BL[, grep("PCA.Loading_LCS", colnames(BL))] <- -0.3
      BL[, grep("Natriumchlorid", colnames(BL))] <- -10
            
      setwd(Dir.ModelDiff)
      #load("Uservorgabe BVLS-Constraints Cmin pH-Diff Monster.rData") 
      #BL.Diff <- Cmin.Range.Diff

      BL.Diff <- as.data.frame(matrix(0, 1, nrow(comp.diff) + 1)) ; colnames(BL.Diff) <- c("Intercept", rownames(comp.diff))
      BL.Diff[, grep("Wasserspektrum", colnames(BL.Diff))] <- -150
      BL.Diff[, grep("PCA.Loading_LCS", colnames(BL.Diff))] <- -0.3
      BL.Diff[, grep("Natriumchlorid", colnames(BL.Diff))] <- -10
                  
      #load("Uservorgabe BVLS-Constraints Cmin Monster ohne Ingredients.rData") 
      #BL.oMI <- as.numeric(Cmin.Range)
      #load("Uservorgabe BVLS-Constraints Cmin Monster mit Ingredients.rData") 
      #BL.mMI <- Cmin.Range
                     
## allgemeine Parameter ---------
   wrange <- which("X" == substr(colnames(as.data.frame(vorlage, TRUE)), 1, 1))
   WaveNums <- colnames(as.data.frame(vorlage, TRUE))[wrange]    
   FI <- as.WaveNum(X3001:X2791, X1900:X991)
   FI.Diff <- as.WaveNum(X3001:X2791, X1900:X991)
   Do.Weight <- TRUE
   bas <- NULL
   pH.fix <- c(3.55, 5.55)            
   max.CO2.warn <- 0.00075
   max.luft.warn <- 0.7   
   max.luft.error <- 1.5   
   #max.WV.warn <- 0.7
   #max.WV.error <- 3
   max.SqS.warn <- 10   
   max.SqS.error <- 40
   dilut.min.warn <- 0.933
   dilut.max.warn <- 0.947
   dilut.min.error <- 0.91
   dilut.max.error <- 0.965
   pH.low.min.warn <- 3.3
   pH.low.max.warn <- 3.7   
   pH.high.min.warn <- 5.3
   pH.high.max.warn <- 5.85 
   pH.low.min.error <- 3
   pH.low.max.error <- 4   
   pH.high.min.error <- 4.8
   pH.high.max.error <- 6.5 
      
## Comp-Daten für Checks ------------
   comp.check <- t(comp.dat[, which("X" == substr(colnames(comp.dat), 1, 1))])
      
## Ableiten comp-Daten nach App-Vorgaben ---------
   Abl.App <- 1
   Gl.App <- 5
   for (i in names(complist.sdb)) {
      complist.sdb[[i]] <- savgol(complist.sdb[[i]], Abl.App, Gl.App, 2)
   }

   for (i in names(complist.svd)) {
      complist.svd[[i]] <- savgol(complist.svd[[i]], Abl.App, Gl.App, 2)
   }
   
   cit.svd <- savgol(cit.svd, Abl.App, Gl.App, 2)      
   comp.dat <- SavGol(comp.dat, Gl.App, 2, Abl.App)
   comp.diff <- SavGol(comp.diff, Gl.App, 2, Abl.App)
                                   
## Workspace speichern --------
   setwd(Dir.ModelDiff)
   VersID <- "130827-01"
   VersionWS <- paste(VersID, "- Workspace Modell Monster KM DIS Quantos0001")
   save.image(file = paste(VersionWS, ".rData", sep = ""))