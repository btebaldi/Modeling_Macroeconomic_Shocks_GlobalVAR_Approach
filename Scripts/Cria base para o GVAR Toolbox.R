#' Autor: Bruno Tebaldi de Queiroz Barbosa
#' 
#' Data: 2025-06-16
#' 
#' Criacao de base de dados para ser utilizado no GVAR
#' 

# Setup -------------------------------------------------------------------

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(lubridate)

# User defined Values -----------------------------------------------------

mDirPath <- "Base de dados"

mFile <- "CAGED - Admt e Deslig - 2004-2024.xlsx"
mAdmOutputFile <- "Adm_OutputFile.xlsx"
mDesOutputFile <- "Des_OutputFile.xlsx"
mPIMOutputFile <- "PIM_OutputFile.xlsx"
mPIMForecastOutputFile <- "PIM_Forecast_OutputFile.xlsx"


mPimFile1 <- "PIM ipeadata 2020-07-30.xlsx"
mPimFile2 <- "PIM ipeadata 2025-06-17.xlsx"


mFile.fullPath <- file.path(mDirPath, mFile)
mPimFile1.fullPath <- file.path(mDirPath, mPimFile1)
mPimFile2.fullPath <- file.path(mDirPath, mPimFile2)

Cadastro_Municipios.path <- file.path(mDirPath, "Cadastro Municipios.xlsx")

# User defined Functions --------------------------------------------------
parseDateForGvar <- function(x){
  ret <- stringr::str_replace(x, pattern = "((Adm)|(Des))_", "") %>%
    lubridate::as_date() 
  ret <- sprintf("%dM%02d", lubridate::year(ret), lubridate::month(ret))
  return(ret)  
}


# Data Load ---------------------------------------------------------------

Data <- list("ALL" = NA, "Adm" = NA, "Des" = NA, PIM_2020 = NA, PIM_2025 = NA, PIM_ALL = NA)

Data$ALL <- readxl::read_xlsx(mFile.fullPath, 
                              sheet = "Admit-Deslig",
                              range = "A2:SM5587")
# range = cell_limits(ul = c(NA, NA), lr = c(NA, NA)))
Data$Adm <- Data$ALL[, c(1:254)]
dim(Data$Adm)
Data$Des <- Data$ALL[, c(1, 2, 256:507)]
dim(Data$Des)

De_Para <- readxl::read_xlsx(file.path(mDirPath, "De-Para IBGE CAGED(e RAIS).xlsx"))
CadastroMunicipios <- list("GrandeRegioes" = NA,
                           "Estados" = NA,
                           "Mesoregioes" = NA,
                           "Microregioes" = NA,
                           "Municipios" = NA,
                           "TabelaCompleta" = NA)

CadastroMunicipios$GrandeRegioes <- readxl::read_xlsx(Cadastro_Municipios.path, sheet = "GrandeRegioes")
CadastroMunicipios$Estados <- readxl::read_xlsx(Cadastro_Municipios.path, sheet = "Estados")
CadastroMunicipios$Mesoregioes <- readxl::read_xlsx(Cadastro_Municipios.path, sheet = "Mesoregioes")
CadastroMunicipios$Microregioes <- readxl::read_xlsx(Cadastro_Municipios.path, sheet = "Microregioes")
CadastroMunicipios$Municipios <- readxl::read_xlsx(Cadastro_Municipios.path, sheet = "Municipios",
                                                   col_types = c("numeric", "numeric", "numeric", "skip", "numeric", "text", "numeric", "text"))
CadastroMunicipios$TabelaCompleta <-  dplyr::left_join(CadastroMunicipios$Municipios,
                                                       CadastroMunicipios$Microregioes,
                                                       by = c("ID_Micro"="ID_Micro") ) %>%
  dplyr::left_join(CadastroMunicipios$Mesoregioes,
                   by = c("ID_Meso"="ID_Meso") ) %>% 
  dplyr::left_join(CadastroMunicipios$Estados,
                   by = c("ID_UF"="ID_UF") ) %>%
  dplyr::left_join(CadastroMunicipios$GrandeRegioes,
                   by = c("ID_GR"="ID_GR") ) %>%
  dplyr::select(ID_GR, Sigla_GR, Nome_GR,
                ID_UF, Sigla_UF, Nome_UF,
                ID_Meso, Nome_Mesorregiao, UF_NomeMesoregiao,
                ID_Micro, Nome_Microrregiao, UF_NomeMicroregiao,
                ID_Municipio, 
                ID_MunicipioNoDigit, Capital, Nome_Municipio, everything())


ShortNameOrder <- c("Ro1101", "Ro1102", "Ac1201", "Ac1202", "Am1301", "Am1302", "Am1303", "Am1304",
                    "Rr1401", "Rr1402", "Pa1501", "Pa1502", "Pa1503", "Pa1504", "Pa1505", "Pa1506",
                    "Ap1601", "Ap1602", "To1701", "To1702", "Ma2101", "Ma2102", "Ma2103", "Ma2104",
                    "Ma2105", "Pi2201", "Pi2202", "Pi2203", "Pi2204", "Ce2301", "Ce2302", "Ce2303",
                    "Ce2304", "Ce2305", "Ce2306", "Ce2307", "Rn2401", "Rn2402", "Rn2403", "Rn2404",
                    "Pb2501", "Pb2502", "Pb2503", "Pb2504", "Pe2601", "Pe2602", "Pe2603", "Pe2604",
                    "Pe2605", "Al2701", "Al2702", "Al2703", "Se2801", "Se2802", "Se2803", "Ba2901",
                    "Ba2902", "Ba2903", "Ba2904", "Ba2905", "Ba2906", "Ba2907", "Mg3101", "Mg3102",
                    "Mg3103", "Mg3104", "Mg3105", "Mg3106", "Mg3107", "Mg3108", "Mg3109", "Mg3110",
                    "Mg3111", "Mg3112", "Es3201", "Es3202", "Es3203", "Es3204", "Rj3301", "Rj3302",
                    "Rj3303", "Rj3304", "Rj3305", "Rj3306", "Sp3501", "Sp3502", "Sp3503", "Sp3504",
                    "Sp3505", "Sp3506", "Sp3507", "Sp3508", "Sp3509", "Sp3510", "Sp3511", "Sp3512",
                    "Sp3513", "Sp3514", "Sp3515", "Pr4101", "Pr4102", "Pr4103", "Pr4104", "Pr4105",
                    "Pr4106", "Pr4107", "Pr4108", "Pr4109", "Pr4110", "Sc4201", "Sc4202", "Sc4203",
                    "Sc4204", "Sc4205", "Sc4206", "Rs4301", "Rs4302", "Rs4303", "Rs4304", "Rs4305",
                    "Rs4306", "Rs4307", "Ms5001", "Ms5002", "Ms5003", "Ms5004", "Mt5101", "Mt5102",
                    "Mt5103", "Mt5104", "Mt5105", "Go5201", "Go5202", "Go5203", "Go5204", "Go5205",
                    "Go5301")


Data$PIM_2020 <- readxl::read_xlsx(mPimFile1.fullPath, sheet = "Séries",
                                   range = cell_limits(ul = c(2, NA), lr = c(NA, NA)))
Data$PIM_2025 <- readxl::read_xlsx(mPimFile2.fullPath, sheet = "Séries")

colnames(Data$PIM_2020) <- c("Data", "PIM")
colnames(Data$PIM_2025) <- c("Data", "PIM")

Data$PIM_2020 %>% tail()

Data$PIM_2020$Data <- lubridate::as_date(Data$PIM_2020$Data, format="%Y.%m")
Data$PIM_2025$Data <- lubridate::as_date(Data$PIM_2025$Data, format="%Y.%m")


# Data Generation - ADM/DES -----------------------------------------------

Data$Adm %>% 
  dplyr::filter( is.na(stringr::str_match(CAGED, "Ignorado")) ) %>% 
  tidyr::pivot_longer(cols = starts_with("Adm_"), names_to = "Periodo", values_to = "Valores") %>% 
  dplyr::mutate(date = parseDateForGvar(Periodo)) %>% 
  dplyr::left_join(De_Para, by =c("CAGED" = "municipio")) %>% 
  dplyr::inner_join(CadastroMunicipios$TabelaCompleta, by=c("codigo" = "ID_Municipio")) %>% 
  dplyr::group_by(Short_name, date) %>% 
  dplyr::summarise(Valor = sum(Valores)) %>% 
  tidyr::pivot_wider(id_cols = date, names_from = "Short_name", values_from = "Valor") %>% 
  dplyr::select(c("date", ShortNameOrder)) %>% 
  writexl::write_xlsx(path = file.path(mDirPath, mAdmOutputFile) )

Data$Des %>% 
  dplyr::filter( is.na(stringr::str_match(CAGED, "Ignorado")) ) %>% 
  tidyr::pivot_longer(cols = starts_with("Des_"), names_to = "Periodo", values_to = "Valores") %>% 
  dplyr::mutate(date = parseDateForGvar(Periodo)) %>% 
  dplyr::left_join(De_Para, by =c("CAGED" = "municipio")) %>% 
  dplyr::inner_join(CadastroMunicipios$TabelaCompleta, by=c("codigo" = "ID_Municipio")) %>% 
  dplyr::group_by(Short_name, date) %>% 
  dplyr::summarise(Valor = sum(Valores)) %>% 
  tidyr::pivot_wider(id_cols = date, names_from = "Short_name", values_from = "Valor") %>% 
  dplyr::select(c("date", ShortNameOrder)) %>% 
  writexl::write_xlsx(path = file.path(mDirPath, mDesOutputFile) )

# Data Generation - PIM ---------------------------------------------------


Data$PIM_ALL <- Data$PIM_2025 %>% left_join(Data$PIM_2020, by = c("Data" = "Data")) %>% mutate(PIM = PIM.y)
# Data$PIM_ALL %>% View()

for(i in 2:nrow(Data$PIM_ALL)){
  if(is.na(Data$PIM_ALL$PIM[i])){
    g <- Data$PIM_ALL$PIM.x[i] / Data$PIM_ALL$PIM.x[i-1]
    Data$PIM_ALL$PIM[i] <- Data$PIM_ALL$PIM[i-1] * g
  }
}

Data$PIM_ALL %>% 
  dplyr::select(Data, PIM) %>% 
  dplyr::mutate(lpim_BR = log(PIM), Date2 = parseDateForGvar(Data)) %>% 
  dplyr::filter(Data >= "2004-01-01") %>% 
  dplyr::select(Date2, lpim_BR) %>% 
  dplyr::rename(date = Date2) %>% 
  writexl::write_xlsx(path = file.path(mDirPath, mPIMOutputFile) )


# Forecasting PIM ---------------------------------------------------------

# Tabela com os forecasts condicionais
PIM_fcast <- tibble(h = -11:96,
                      date = as.Date(NA),
                      lPIM_0pc = NA,
                      lPIM_2pc = NA,
                      lPIM_4pc = NA,
                      lPIM_6pc = NA,
                      lPIM_8pc = NA)

# define growth rates
tx_2pc <- log(1.02)
tx_4pc <- log(1.04)
tx_6pc <- log(1.06)
tx_8pc <- log(1.08)

# Set current date
cur_date <- as.Date("2016-12-01")

# calculate the forecast dates for every horizon h
PIM_fcast$date <- cur_date
PIM_fcast <- PIM_fcast %>%
  dplyr::mutate(date = lubridate::add_with_rollback(date, months(h)))

# get the previous values
idx_ini <- which(Data$PIM_ALL$Data == lubridate::add_with_rollback(cur_date, months(-11)))
idx_fim <- which(Data$PIM_ALL$Data == cur_date)
PIM_fcast$lPIM_0pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast$lPIM_2pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast$lPIM_4pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast$lPIM_6pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast$lPIM_8pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])

# Initialize the 0% forecast scenario
for(i in 13:nrow(PIM_fcast)){
  PIM_fcast$lPIM_0pc[i] <- PIM_fcast$lPIM_0pc[i-12]
  PIM_fcast$lPIM_2pc[i] <- PIM_fcast$lPIM_2pc[i-12] + tx_2pc
  PIM_fcast$lPIM_4pc[i] <- PIM_fcast$lPIM_4pc[i-12] + tx_4pc
  PIM_fcast$lPIM_6pc[i] <- PIM_fcast$lPIM_6pc[i-12] + tx_6pc
  PIM_fcast$lPIM_8pc[i] <- PIM_fcast$lPIM_8pc[i-12] + tx_8pc
}

PIM_fcast <- PIM_fcast %>% dplyr::filter(h >0)




# Fazer previsoes fora da amostraa ----------------------------------------


# Tabela com os forecasts condicionais
PIM_fcast_2024 <- tibble(h = -11:96,
                    date = as.Date(NA),
                    lPIM_0pc = NA,
                    lPIM_2pc = NA,
                    lPIM_4pc = NA,
                    lPIM_6pc = NA,
                    lPIM_8pc = NA)

# define growth rates
tx_2pc <- log(1.02)
tx_4pc <- log(1.04)
tx_6pc <- log(1.06)
tx_8pc <- log(1.08)

# Set current date
cur_date <- as.Date("2024-12-01")

# calculate the forecast dates for every horizon h
PIM_fcast_2024$date <- cur_date
PIM_fcast_2024 <- PIM_fcast_2024 %>%
  dplyr::mutate(date = lubridate::add_with_rollback(date, months(h)))

# get the previous values
idx_ini <- which(Data$PIM_ALL$Data == lubridate::add_with_rollback(cur_date, months(-11)))
idx_fim <- which(Data$PIM_ALL$Data == cur_date)
PIM_fcast_2024$lPIM_0pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast_2024$lPIM_2pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast_2024$lPIM_4pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast_2024$lPIM_6pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])
PIM_fcast_2024$lPIM_8pc[1:12] <- log(Data$PIM_ALL$PIM[idx_ini:idx_fim])

# Initialize the 0% forecast scenario
for(i in 13:nrow(PIM_fcast_2024)){
  PIM_fcast_2024$lPIM_0pc[i] <- PIM_fcast_2024$lPIM_0pc[i-12]
  PIM_fcast_2024$lPIM_2pc[i] <- PIM_fcast_2024$lPIM_2pc[i-12] + tx_2pc
  PIM_fcast_2024$lPIM_4pc[i] <- PIM_fcast_2024$lPIM_4pc[i-12] + tx_4pc
  PIM_fcast_2024$lPIM_6pc[i] <- PIM_fcast_2024$lPIM_6pc[i-12] + tx_6pc
  PIM_fcast_2024$lPIM_8pc[i] <- PIM_fcast_2024$lPIM_8pc[i-12] + tx_8pc
}

PIM_fcast_2024 <- PIM_fcast_2024 %>% dplyr::filter(h >0)




tail(PIM_fcast)
tail(PIM_fcast_2024)

writexl::write_xlsx(x = list("Fcast2016" = PIM_fcast, "Fcast2024" = PIM_fcast_2024),
                    path = file.path(mDirPath, mPIMForecastOutputFile) )
