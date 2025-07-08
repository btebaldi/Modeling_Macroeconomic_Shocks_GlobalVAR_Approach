# Setup -------------------------------------------------------------------
rm(list = ls())


library(readxl)
library(tidyverse)


# Internal Variabels ------------------------------------------------------
# Dir and files of results
mDirMask <- "GVAR_Toolbox2.0/Output/Meso17_FullSample_%dperc"
mFileMask <- "m17_full_%dperc_output.xlsx"

# Actual Values
# mActualDir <- "Base de dados"
# mAdmFile <- "Adm_OutputFile.xlsx"
# mDesFile <- "Des_OutputFile.xlsx"


mOutputGraphDir <- "Export/Graph/Results Full Sample"


# Data Load ---------------------------------------------------------------

# load Actual values
# tblAdmActual <- readxl::read_excel(path = file.path(mActualDir, mAdmFile))
# tblDeslActual <- readxl::read_excel(path = file.path(mActualDir, mDesFile))

percent <- c(0,2,4,6,8)
# percent <- c(2, 8)
list_of_models <- list()

for(i in seq_along(percent)){
  percent_current <- percent[i]
  cat("Lendo resultados",  percent_current, "%\n")
  
  # Abretura dos resultados da previsao condicional
  tbl <- readxl::read_excel(path = file.path(sprintf(mDirMask, percent_current), sprintf(mFileMask, percent_current)),
                            sheet = "conditional forecasts",
                            range = cell_limits(ul = c(5, 1), lr = c(NA, NA)))
  
  # Regularizacao de colunas faltantes
  colnames(tbl)[1] <- "Region"
  colnames(tbl)[2] <- "Series"
  colnames(tbl)[3] <- "F_A"
  
  # Regularizacao de colunas faltantes
  tbl <- tbl %>% 
    dplyr::select(- F_A) %>% 
    dplyr::filter(!is.na(Region)) %>% 
    tidyr::pivot_longer(cols = starts_with("20"), names_to = "Periodo") %>% 
    dplyr::mutate(Date = lubridate::ymd(Periodo, truncated = 1))
  
  list_of_models[[i]] <- tbl 
}

# list_of_models[[length(list_of_models)+1]] <- tblAdmActual %>% 
#   pivot_longer(cols = -date) %>% 
#   dplyr::mutate(Series = "adm", 
#                 date = lubridate::ymd(date, truncated = 1)) %>%
#   dplyr::rename(Actual = value)
# 
# 
# list_of_models[[length(list_of_models)+1]] <- tblDeslActual %>% 
#   pivot_longer(cols = -date) %>% 
#   dplyr::mutate(Series = "desl", 
#                 date = lubridate::ymd(date, truncated = 1)) %>%
#   dplyr::rename(Actual = value)


# names(list_of_models) <- c(sprintf("Perc%d", percent), c("ActualAdm", "ActualDesl"))
names(list_of_models) <- sprintf("Perc%d", percent)



# Data interpretation -----------------------------------------------------

i=1
for(i in seq_along(percent)){
  percent_current <- percent[i]
  cat("Lendo resultados",  percent_current, "%\n")
  
  tbl <- list_of_models[[i]]
  newNameCol <- sprintf("Perc_%d", percent_current)
  tbl <- tbl %>% dplyr::rename(!!newNameCol := value) %>% dplyr::select(-Periodo) %>% 
    dplyr::filter(Date < as.Date("2029-01-01"))
  
  
  tbl$Region[tbl$Series == "Sp3515A"] <- "Sp3515"
  tbl$Series[tbl$Series == "Sp3515A"] <- "adm"
  
  tbl$Region[tbl$Series == "Sp3515D"] <- "Sp3515"
  tbl$Series[tbl$Series == "Sp3515D"] <- "desl"
  
  
  if(i == 1){
    tblfull <- tbl 
  } else {
    tblfull <- inner_join(tblfull, tbl, by = c("Region"="Region",
                                               "Series"="Series",
                                               "Date" = "Date") )
  }
}


# tblfull <- dplyr::left_join(tblfull, 
#                             dplyr::bind_rows(list_of_models[["ActualAdm"]],list_of_models[["ActualDesl"]]),
#                             by = c("Region" = "name", 
#                                    "Date"="date",
#                                    "Series"="Series"))


# # list of reagions
mListRegions <- c( "RJ" = "Rj3306",
                   "BH" = "Mg3107",
                   "Salvador" = "Ba2905",
                   "SP" = "Sp3515",
                   "PortoAlegre" = "Rs4305",
                   "Jequitinhonha" = "Mg3103",
                   "Borborema" = "Pb2502",
                   "CentralGoias" = "Go5203")
i=1
for(i in seq_along(mListRegions)){
  
  region <- mListRegions[i]
  
  tbl2 <- tblfull %>% 
    dplyr::filter(Region %in% region) %>%
    tidyr::pivot_longer(cols = sprintf("Perc_%d", percent), names_to = "Tipo") %>% 
    tidyr::pivot_wider(id_cols = c("Region", "Date", "Tipo"),
                       names_from = "Series",
                       values_from = "value") %>% 
    dplyr::arrange(Date) %>% 
    dplyr::group_by(Region, Tipo) %>% 
    dplyr::mutate(Net = adm - desl,
                  NetCumSum = cumsum(Net))
  
  graph <- tbl2 %>% 
    ggplot() + 
    geom_line(aes(x = Date, y = NetCumSum, colour = Tipo, linetype = Tipo)) + 
    labs(title = "Cumulative employment rate",
         subtitle = sprintf("%s regional level",  names(mListRegions)[i]),
         x = NULL, y = NULL, colour = NULL, linetype= NULL) +
    theme_bw() + 
    theme(legend.position = "bottom") +
    scale_color_manual(breaks = c("Actual", sprintf("Perc_%d", c(0:4)*2)),
                       values = c("#000000", "#F8766D", "#B79F00",
                                  "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))
  
  print(graph)
  ggsave(filename = sprintf("Region %s.png", region),
         path = mOutputGraphDir,
         scale = 1,
         width = 8, height = 6, units = "in", dpi = 200)
  
}



# Agregacao nivel Brasil --------------------------------------------------


tbl2 <- tblfull %>% 
  tidyr::pivot_longer(cols = sprintf("Perc_%d", percent), names_to = "Tipo") %>% 
  dplyr::group_by(Series, Date, Tipo) %>%
  dplyr::summarise(value = sum(value)) %>% 
  tidyr::pivot_wider(id_cols = c("Date", "Tipo"),
                     names_from = "Series",
                     values_from = "value") %>% 
  dplyr::arrange(Date) %>% 
  dplyr::group_by(Tipo) %>% 
  dplyr::mutate(Net = adm - desl,
                NetCumSum = cumsum(Net))

graph <- tbl2 %>% 
  ggplot() + 
  geom_line(aes(x = Date, y = NetCumSum, colour = Tipo, linetype = Tipo)) + 
  labs(title = "Cumulative employment forecast",
       subtitle = "Brazilian level",
       y=NULL, x=NULL, colour = NULL, linetype = NULL) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  scale_color_manual(breaks = c("Actual", sprintf("Perc_%d", c(0:4)*2)),
                     values = c("#000000", "#F8766D", "#B79F00",
                                "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))

print(graph)
ggsave(filename = "Brazil net employment forecast.png",
       path = mOutputGraphDir,
       scale = 1,
       width = 8, height = 6, units = "in", dpi = 200)


tail(tbl2)


# Big Region Classification ----------------------------------------------------


# Agregacao NO=North; NE=Northeast; MW=Midwest; SE=Southeast; SO=South
my_label_n_level <- c(
  "Ro1101" = "NO",
  "Ro1102" = "NO",
  "Ac1201" = "NO",
  "Ac1202" = "NO",
  "Am1301" = "NO",
  "Am1302" = "NO",
  "Am1303" = "NO",
  "Am1304" = "NO",
  "Rr1401" = "NO",
  "Rr1402" = "NO",
  "Pa1501" = "NO",
  "Pa1502" = "NO",
  "Pa1503" = "NO",
  "Pa1504" = "NO",
  "Pa1505" = "NO",
  "Pa1506" = "NO",
  "Ap1601" = "NO",
  "Ap1602" = "NO",
  "To1701" = "NO",
  "To1702" = "NO",
  # 
  "Ma2101" = "NE",
  "Ma2102" = "NE",
  "Ma2103" = "NE",
  "Ma2104" = "NE",
  "Ma2105" = "NE",
  "Pi2201" = "NE",
  "Pi2202" = "NE",
  "Pi2203" = "NE",
  "Pi2204" = "NE",
  "Ce2301" = "NE",
  "Ce2302" = "NE",
  "Ce2303" = "NE",
  "Ce2304" = "NE",
  "Ce2305" = "NE",
  "Ce2306" = "NE",
  "Ce2307" = "NE",
  "Rn2401" = "NE",
  "Rn2402" = "NE",
  "Rn2403" = "NE",
  "Rn2404" = "NE",
  # 
  "Pb2501" = "NE",
  "Pb2502" = "NE",
  "Pb2503" = "NE",
  "Pb2504" = "NE",
  "Pe2601" = "NE",
  "Pe2602" = "NE",
  "Pe2603" = "NE",
  "Pe2604" = "NE",
  "Pe2605" = "NE",
  "Al2701" = "NE",
  
  "Al2702" = "NE",
  "Al2703" = "NE",
  "Se2801" = "NE",
  "Se2802" = "NE",
  "Se2803" = "NE",
  "Ba2901" = "NE",
  "Ba2902" = "NE",
  "Ba2903" = "NE",
  "Ba2904" = "NE",
  "Ba2905" = "NE",
  # 
  "Ba2906" = "NE",
  "Ba2907" = "NE",
  # 
  "Mg3101" = "SE",
  "Mg3102" = "SE",
  "Mg3103" = "SE",
  "Mg3104" = "SE",
  "Mg3105" = "SE",
  "Mg3106" = "SE",
  "Mg3107" = "SE",
  "Mg3108" = "SE",
  # 
  "Mg3109" = "SE",
  "Mg3110" = "SE",
  "Mg3111" = "SE",
  "Mg3112" = "SE",
  "Es3201" = "SE",
  "Es3202" = "SE",
  "Es3203" = "SE",
  "Es3204" = "SE",
  "Rj3301" = "SE",
  "Rj3302" = "SE",
  # 
  "Rj3303" = "SE",
  "Rj3304" = "SE",
  "Rj3305" = "SE",
  "Rj3306" = "SE",
  "Sp3501" = "SE",
  "Sp3502" = "SE",
  "Sp3503" = "SE",
  "Sp3504" = "SE",
  "Sp3505" = "SE",
  "Sp3506" = "SE",
  # 
  "Sp3507" = "SE",
  "Sp3508" = "SE",
  "Sp3509" = "SE",
  "Sp3510" = "SE",
  "Sp3511" = "SE",
  "Sp3512" = "SE",
  "Sp3513" = "SE",
  "Sp3514" = "SE",
  "Sp3515" = "SE",
  "Pr4101" = "SO",
  # 
  "Pr4102" = "SO",
  "Pr4103" = "SO",
  "Pr4104" = "SO",
  "Pr4105" = "SO",
  "Pr4106" = "SO",
  "Pr4107" = "SO",
  "Pr4108" = "SO",
  "Pr4109" = "SO",
  "Pr4110" = "SO",
  "Sc4201" = "SO",
  # 
  "Sc4202" = "SO",
  "Sc4203" = "SO",
  "Sc4204" = "SO",
  "Sc4205" = "SO",
  "Sc4206" = "SO",
  "Rs4301" = "SO",
  "Rs4302" = "SO",
  "Rs4303" = "SO",
  "Rs4304" = "SO",
  "Rs4305" = "SO",
  # 
  "Rs4306" = "SO",
  "Rs4307" = "SO",
  "Ms5001" = "MW",
  "Ms5002" = "MW",
  "Ms5003" = "MW",
  "Ms5004" = "MW",
  "Mt5101" = "MW",
  "Mt5102" = "MW",
  "Mt5103" = "MW",
  "Mt5104" = "MW",
  # 
  "Mt5105" = "MW",
  "Go5201" = "MW",
  "Go5202" = "MW",
  "Go5203" = "MW",
  "Go5204" = "MW",
  "Go5205" = "MW",
  "Go5301" = "MW",
  "du_model" = "DU")


# Aggregation by Big Region ----------------------------------------------------

graph <- tblfull %>% 
  dplyr::mutate(MajorRegion = factor(x = Region, levels = names(my_label_n_level),  labels = my_label_n_level)) %>%
  tidyr::pivot_longer(cols = c("Perc_0", "Perc_2", "Perc_4", "Perc_6", "Perc_8"),
                      names_to = "Tipo", values_to = "value") %>% 
  dplyr::group_by(MajorRegion, Series, Date, Tipo) %>%
  dplyr::summarise(value = sum(value),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(id_cols = c("MajorRegion", "Date", "Tipo"), 
                     names_from = Series,
                     values_from = value) %>% 
  dplyr::arrange(Date) %>% 
  dplyr::group_by(MajorRegion, Tipo) %>% 
  dplyr::mutate(Net = adm - desl,
                NetCumSum = cumsum(Net)) %>% 
  dplyr::filter(MajorRegion != "DU") %>% 
  ggplot() + 
  geom_line(aes(x = Date, y = NetCumSum, colour = Tipo)) + 
  facet_wrap(~MajorRegion, scales="free_y") +
  # facet_grid(.~MajorRegion, scales="free_y") +
  labs(title = "Cumulative employment forecast",
       subtitle = "Major Region level",
       x = NULL, y = NULL, colour = NULL) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  scale_color_manual(breaks = c("Actual", sprintf("Perc_%d", c(0:4)*2)),
                     values = c("#000000", "#F8766D", "#B79F00",
                                "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))


print(graph)
ggsave(filename = "Major Regions net employment forecast.png",
       path = mOutputGraphDir,
       scale = 1,
       width = 8, height = 6, units = "in", dpi = 200)


# State Classification----------------------------------------------------------

my_label_n_level <- c(
  "Ro1101" = "Rondônia",
  "Ro1102" = "Rondônia",
  "Ac1201" = "Acre",
  "Ac1202" = "Acre",
  "Am1301" = "Amazonas",
  "Am1302" = "Amazonas",
  "Am1303" = "Amazonas",
  "Am1304" = "Amazonas",
  "Rr1401" = "Roraima",
  "Rr1402" = "Roraima",
  "Pa1501" = "Pará",
  "Pa1502" = "Pará",
  "Pa1503" = "Pará",
  "Pa1504" = "Pará",
  "Pa1505" = "Pará",
  "Pa1506" = "Pará",
  "Ap1601" = "Amapá",
  "Ap1602" = "Amapá",
  "To1701" = "Tocantins",
  "To1702" = "Tocantins",
  "Ma2101" = "Maranhão",
  "Ma2102" = "Maranhão",
  "Ma2103" = "Maranhão",
  "Ma2104" = "Maranhão",
  "Ma2105" = "Maranhão",
  "Pi2201" = "Piauí",
  "Pi2202" = "Piauí",
  "Pi2203" = "Piauí",
  "Pi2204" = "Piauí",
  "Ce2301" = "Ceará",
  "Ce2302" = "Ceará",
  "Ce2303" = "Ceará",
  "Ce2304" = "Ceará",
  "Ce2305" = "Ceará",
  "Ce2306" = "Ceará",
  "Ce2307" = "Ceará",
  "Rn2401" = "Rio Grande do Norte",
  "Rn2402" = "Rio Grande do Norte",
  "Rn2403" = "Rio Grande do Norte",
  "Rn2404" = "Rio Grande do Norte",
  "Pb2501" = "Paraíba",
  "Pb2502" = "Paraíba",
  "Pb2503" = "Paraíba",
  "Pb2504" = "Paraíba",
  "Pe2601" = "Pernambuco",
  "Pe2602" = "Pernambuco",
  "Pe2603" = "Pernambuco",
  "Pe2604" = "Pernambuco",
  "Pe2605" = "Pernambuco",
  "Al2701" = "Alagoas",
  "Al2702" = "Alagoas",
  "Al2703" = "Alagoas",
  "Se2801" = "Sergipe",
  "Se2802" = "Sergipe",
  "Se2803" = "Sergipe",
  "Ba2901" = "Bahia",
  "Ba2902" = "Bahia",
  "Ba2903" = "Bahia",
  "Ba2904" = "Bahia",
  "Ba2905" = "Bahia",
  "Ba2906" = "Bahia",
  "Ba2907" = "Bahia",
  "Mg3101" = "Minas Gerais",
  "Mg3102" = "Minas Gerais",
  "Mg3103" = "Minas Gerais",
  "Mg3104" = "Minas Gerais",
  "Mg3105" = "Minas Gerais",
  "Mg3106" = "Minas Gerais",
  "Mg3107" = "Minas Gerais",
  "Mg3108" = "Minas Gerais",
  "Mg3109" = "Minas Gerais",
  "Mg3110" = "Minas Gerais",
  "Mg3111" = "Minas Gerais",
  "Mg3112" = "Minas Gerais",
  "Es3201" = "Espírito Santo",
  "Es3202" = "Espírito Santo",
  "Es3203" = "Espírito Santo",
  "Es3204" = "Espírito Santo",
  "Rj3301" = "Rio de Janeiro",
  "Rj3302" = "Rio de Janeiro",
  "Rj3303" = "Rio de Janeiro",
  "Rj3304" = "Rio de Janeiro",
  "Rj3305" = "Rio de Janeiro",
  "Rj3306" = "Rio de Janeiro",
  "Sp3501" = "São Paulo",
  "Sp3502" = "São Paulo",
  "Sp3503" = "São Paulo",
  "Sp3504" = "São Paulo",
  "Sp3505" = "São Paulo",
  "Sp3506" = "São Paulo",
  "Sp3507" = "São Paulo",
  "Sp3508" = "São Paulo",
  "Sp3509" = "São Paulo",
  "Sp3510" = "São Paulo",
  "Sp3511" = "São Paulo",
  "Sp3512" = "São Paulo",
  "Sp3513" = "São Paulo",
  "Sp3514" = "São Paulo",
  "Sp3515" = "São Paulo",
  "Pr4101" = "Paraná",
  "Pr4102" = "Paraná",
  "Pr4103" = "Paraná",
  "Pr4104" = "Paraná",
  "Pr4105" = "Paraná",
  "Pr4106" = "Paraná",
  "Pr4107" = "Paraná",
  "Pr4108" = "Paraná",
  "Pr4109" = "Paraná",
  "Pr4110" = "Paraná",
  "Sc4201" = "Santa Catarina",
  "Sc4202" = "Santa Catarina",
  "Sc4203" = "Santa Catarina",
  "Sc4204" = "Santa Catarina",
  "Sc4205" = "Santa Catarina",
  "Sc4206" = "Santa Catarina",
  "Rs4301" = "Rio Grande do Sul",
  "Rs4302" = "Rio Grande do Sul",
  "Rs4303" = "Rio Grande do Sul",
  "Rs4304" = "Rio Grande do Sul",
  "Rs4305" = "Rio Grande do Sul",
  "Rs4306" = "Rio Grande do Sul",
  "Rs4307" = "Rio Grande do Sul",
  "Ms5001" = "Mato Grosso do Sul",
  "Ms5002" = "Mato Grosso do Sul",
  "Ms5003" = "Mato Grosso do Sul",
  "Ms5004" = "Mato Grosso do Sul",
  "Mt5101" = "Mato Grosso",
  "Mt5102" = "Mato Grosso",
  "Mt5104" = "Mato Grosso",
  "Mt5103" = "Mato Grosso",
  "Mt5105" = "Mato Grosso",
  "Go5201" = "Goiás",
  "Go5202" = "Goiás",
  "Go5203" = "Goiás",
  "Go5204" = "Goiás",
  "Go5205" = "Goiás",
  "Go5301" = "Goiás", # Distrito Federal
  "du_model" = "Dominat Unit")

# Aggregation by State ----------------------------------------------------


tbl2 <- tblfull %>% 
  dplyr::mutate(MajorRegion = factor(x = Region, levels = names(my_label_n_level),  labels = my_label_n_level)) %>%
  tidyr::pivot_longer(cols = c("Perc_0", "Perc_2", "Perc_4", "Perc_6", "Perc_8"),
                      names_to = "Tipo", values_to = "value") %>% 
  dplyr::group_by(MajorRegion, Series, Date, Tipo) %>%
  dplyr::summarise(value = sum(value),
                   .groups = "drop") %>% 
  tidyr::pivot_wider(id_cols = c("MajorRegion", "Date", "Tipo"), 
                     names_from = Series,
                     values_from = value) %>% 
  dplyr::arrange(Date) %>% 
  dplyr::group_by(MajorRegion, Tipo) %>% 
  dplyr::mutate(Net = adm - desl,
                NetCumSum = cumsum(Net))

for(estado in c("São Paulo", "Minas Gerais", "Espírito Santo", "Rio de Janeiro")){
  
  graph <- tbl2 %>% 
    dplyr::filter(MajorRegion == estado) %>% 
    ggplot() + 
    geom_line(aes(x = Date, y = NetCumSum, colour = Tipo)) + 
    # facet_wrap(~MajorRegion, scales="free_y") +
    # facet_grid(.~MajorRegion, scales="free_y") +
    labs(title = "Cumulative employment forecast",
         subtitle = sprintf("%s State Region level", estado),
         x = NULL, y = NULL, colour = NULL) +
    theme_bw() + 
    theme(legend.position = "bottom") +
    scale_color_manual(breaks = c("Actual", sprintf("Perc_%d", c(0:4)*2)),
                       values = c("#000000", "#F8766D", "#B79F00",
                                  "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))
  
  
  print(graph)
  ggsave(filename = sprintf("State %s Region net employment forecast.png", estado),
         path = mOutputGraphDir,
         scale = 1,
         width = 8, height = 6, units = "in", dpi = 200)
}
