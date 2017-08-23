#' Functional information tool
#'
#' Makes raw MiDAS function data compatible with ampvis format. Internal function, not exported.
#'
#' @usage amp_cleanMiF(data)
#'
#' @param data (required) A data frame with MiDAS functions.
#' @import tidyverse
#' @return A data frame.
#' 
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_cleanMiF <- function(data){
  
  MiF <- mutate(data,
                MiDAS = "POS",
                
                FIL = paste(Filamentous.Other,Filamentous.In.situ),
                FIL = ifelse(FIL %in% c("POS POS","NEG POS","NT POS","POS NT", "POS VAR", "VAR POS"), "POS", FIL),
                FIL = ifelse(FIL %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", FIL),
                FIL = ifelse(FIL %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", FIL),
                FIL = ifelse(FIL == "NT NT", "NT", FIL),
                
                AOB = paste(AOB.Other,AOB.In.situ),
                AOB = ifelse(AOB %in% c("POS POS","NEG POS","NT POS","POS NT", "POS VAR", "VAR POS"), "POS", AOB),
                AOB = ifelse(AOB %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", AOB),
                AOB = ifelse(AOB %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", AOB),
                AOB = ifelse(AOB == "NT NT", "NT", AOB),
                
                NOB = paste(NOB.Other,NOB.In.situ),
                NOB = ifelse(NOB %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", NOB),
                NOB = ifelse(NOB %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", NOB),
                NOB = ifelse(NOB %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", NOB),               
                NOB = ifelse(NOB == "NT NT", "NT", NOB),
                
                Anammox = paste(Anammox.Other,Anammox.In.situ),
                Anammox = ifelse(Anammox %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", Anammox),
                Anammox = ifelse(Anammox %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", Anammox),
                Anammox = ifelse(Anammox %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", Anammox),               
                Anammox = ifelse(Anammox == "NT NT", "NT", Anammox),
                
                AU.MIX = paste(Autotroph.Mixotroph.Other,Autotroph.Mixotroph.In.situ),
                AU.MIX = ifelse(AU.MIX %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", AU.MIX),
                AU.MIX = ifelse(AU.MIX %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", AU.MIX),
                AU.MIX = ifelse(AU.MIX %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", AU.MIX),               
                AU.MIX = ifelse(AU.MIX == "NT NT", "NT", AU.MIX),
                
                PAO = paste(PAO.Other,PAO.In.situ),
                PAO = ifelse(PAO %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", PAO),
                PAO = ifelse(PAO %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", PAO),
                PAO = ifelse(PAO %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", PAO),               
                PAO = ifelse(PAO == "NT NT", "NT", PAO),
                
                GAO = paste(GAO.Other,GAO.In.situ),
                GAO = ifelse(GAO %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", GAO),
                GAO = ifelse(GAO %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", GAO),
                GAO = ifelse(GAO %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", GAO),               
                GAO = ifelse(GAO == "NT NT", "NT", GAO),
                
                HET = paste(Aerobic.heterotroph.Other,Aerobic.heterotroph.In.situ),
                HET = ifelse(HET %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", HET),
                HET = ifelse(HET %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", HET),
                HET = ifelse(HET %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", HET),               
                HET = ifelse(HET == "NT NT", "NT", HET),
                
                DN = paste(Nitrite.reduction.Other,Nitrite.reduction.In.situ),
                DN = ifelse(DN %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", DN),
                DN = ifelse(DN %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", DN),
                DN = ifelse(DN %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", DN),               
                DN = ifelse(DN == "NT NT", "NT", DN),   
                
                FER = paste(Fermentation.Other,Fermentation.In.situ),
                FER = ifelse(FER %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", FER),
                FER = ifelse(FER %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", FER),
                FER = ifelse(FER %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", FER),               
                FER = ifelse(FER == "NT NT", "NT", FER),  
                
                SUL = paste(Sulphate.reduction.Other,Sulphate.reduction.In.situ),
                SUL = ifelse(SUL %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", SUL),
                SUL = ifelse(SUL %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", SUL),
                SUL = ifelse(SUL %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", SUL),               
                SUL = ifelse(SUL == "NT NT", "NT", SUL),
                
                ACE = paste(Acetogen.Other,Acetogen.In.situ),
                ACE = ifelse(ACE %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", ACE),
                ACE = ifelse(ACE %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", ACE),
                ACE = ifelse(ACE %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", ACE),               
                ACE = ifelse(ACE == "NT NT", "NT", ACE),
                
                MET = paste(Methanogen.Other,Methanogen.In.situ),
                MET = ifelse(MET %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", MET),
                MET = ifelse(MET %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", MET),
                MET = ifelse(MET %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", MET),               
                MET = ifelse(MET == "NT NT", "NT", MET),
                
                FA = paste(Fatty.acids.Other,Fatty.acids.In.situ),
                FA = ifelse(FA %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", FA),
                FA = ifelse(FA %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", FA),
                FA = ifelse(FA %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", FA),               
                FA = ifelse(FA == "NT NT", "NT", FA),
                
                SUG = paste(Sugars.Other,Sugars.In.situ),
                SUG = ifelse(SUG %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", SUG),
                SUG = ifelse(SUG %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", SUG),
                SUG = ifelse(SUG %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", SUG),               
                SUG = ifelse(SUG == "NT NT", "NT", SUG),
                
                PRO = paste(Proteins.Amino.acids.Other,Proteins.Amino.acids.In.situ),
                PRO = ifelse(PRO %in% c("POS POS","NEG POS","NT POS","POS NT" , "POS VAR" , "VAR POS"), "POS", PRO),
                PRO = ifelse(PRO %in% c("NEG NEG","POS NEG","NT NEG","NEG NT", "VAR NEG"), "NEG", PRO),
                PRO = ifelse(PRO %in% c("VAR NT","NT VAR","NEG VAR", "VAR VAR"), "VAR", PRO),               
                PRO = ifelse(PRO == "NT NT", "NT", PRO)
  )
  return(MiF)
}
