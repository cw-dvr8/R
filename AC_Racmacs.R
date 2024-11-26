library(tidyverse)
library(dplyr)

get(load("/Users/cindy/projects/antigenic_cartography/data/input_files/An_Dataset_BL.Rdata"))
get(load("/Users/cindy/projects/antigenic_cartography/data/input_files/An_Dataset_D29.Rdata"))
get(load("/Users/cindy/projects/antigenic_cartography/data/input_files/An_Dataset_PtidData.Rdata"))

# biggab_0 <- biggab_0 %>% arrange(ptid, antigen)
biggab_1 <- biggab_1 %>% arrange(ptid, antigen)

biggab_1_filtered <- biggab_1 %>% filter(ptid != 808991795)

titer_ap_df <- biggab_1_filtered %>%
               select(antigen, ptid, result) %>% 
               pivot_wider(names_from="ptid", values_from="result")

titer_ap_df <- titer_ap_df %>% remove_rownames %>% column_to_rownames(var="antigen")
titer_ap_rfiltered_df <- titer_ap_df %>% drop_na()
titer_ap_cfiltered_df <- titer_ap_df %>% select_if(~ !any(is.na(.)))

titer_pa_df <- biggab_0 %>%
               select(antigen, ptid, result) %>% 
               pivot_wider(names_from="antigen", values_from="result")
titer_pa_df$ptid <- as.character(titer_pa_df$ptid)

titer_pa_df <- titer_pa_df %>% remove_rownames %>% column_to_rownames(var="ptid")
titer_pa_rfiltered_df <- titer_pa_df %>% drop_na()
titer_pa_cfiltered_df <- titer_pa_df %>% select_if(~ !any(is.na(.)))

ptid_boost_df <- ptdat %>% 
    mutate(vac_group=paste(Vac_Prime, Vac_Boost, sep="_")) %>%
    select(ptid, vac_group)

jan_jan_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Jan_Jan"]
jan_mod_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Jan_Mod"]
jan_pfi_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Jan_Pfi"]
jan_mod211_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Jan_Mod.211"]
jan_nvx_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Jan_Nvx"]
mod_jan_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Mod_Jan"]
mod_mod_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Mod_Mod"]
mod_pfi_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Mod_Pfi"]
mod_nvx_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Mod_Nvx"]
pfi_jan_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Pfi_Jan"]
pfi_mod_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Pfi_Mod"]
pfi_pfi_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Pfi_Pfi"]
pfi_mod211_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Pfi_Mod.211"]
pfi_nvx_vt <- ptid_boost_df$ptid[ptid_boost_df$vac_group == "Pfi_Nvx"]
# mod_vt <- ptid_boost_df$ptid[ptid_boost_df$Vac_Boost == "Mod"]
# pfizer_vt <- ptid_boost_df$ptid[ptid_boost_df$Vac_Boost == "Pfi"]
# mod211_vt <- ptid_boost_df$ptid[ptid_boost_df$Vac_Boost == "Mod.211"]
# janssen_vt <- ptid_boost_df$ptid[ptid_boost_df$Vac_Boost == "Jan"]
# nvx_vt <- ptid_boost_df$ptid[ptid_boost_df$Vac_Boost == "Nvx"]

omicron_vt <- c("SARSCOV2SP-B.1.1.529")

library(Racmacs)

ap_rfiltered_map <- acmap(titer_table = titer_ap_rfiltered_df)
ap_rfiltered_map <- optimizeMap(
    map                     = ap_rfiltered_map,
    number_of_dimensions    = 2,
    number_of_optimizations = 1500,
    minimum_column_basis    = "none"
)

jj_ptids <- srNames(ap_rfiltered_map) %in% jan_jan_vt
jm_ptids <- srNames(ap_rfiltered_map) %in% jan_mod_vt
jm211_ptids <- srNames(ap_rfiltered_map) %in% jan_mod211_vt
jn_ptids <- srNames(ap_rfiltered_map) %in% jan_nvx_vt
jp_ptids <- srNames(ap_rfiltered_map) %in% jan_pfi_vt
mj_ptids <- srNames(ap_rfiltered_map) %in% mod_jan_vt
mm_ptids <- srNames(ap_rfiltered_map) %in% mod_mod_vt
mn_ptids <- srNames(ap_rfiltered_map) %in% mod_nvx_vt
mp_ptids <- srNames(ap_rfiltered_map) %in% mod_pfi_vt
pj_ptids <- srNames(ap_rfiltered_map) %in% pfi_jan_vt
pm_ptids <- srNames(ap_rfiltered_map) %in% pfi_mod_vt
pm211_ptids <- srNames(ap_rfiltered_map) %in% pfi_mod211_vt
pn_ptids <- srNames(ap_rfiltered_map) %in% pfi_nvx_vt
pp_ptids <- srNames(ap_rfiltered_map) %in% pfi_pfi_vt
# mod_ptids <- srNames(ap_rfiltered_map) %in% mod_vt
# mod211_ptids <- srNames(ap_rfiltered_map) %in% mod211_vt
# pfizer_ptids <- srNames(ap_rfiltered_map) %in% pfizer_vt
# janssen_ptids <- srNames(ap_rfiltered_map) %in% janssen_vt
# nvx_ptids <- srNames(ap_rfiltered_map) %in% nvx_vt
srFill(ap_rfiltered_map)[jj_ptids] <- "skyblue"
srFill(ap_rfiltered_map)[jm_ptids] <- "darkgreen"
srFill(ap_rfiltered_map)[jm211_ptids] <- "maroon4"
srFill(ap_rfiltered_map)[jn_ptids] <- "orange"
srFill(ap_rfiltered_map)[jp_ptids] <- "hotpink"
srFill(ap_rfiltered_map)[mj_ptids] <- "wheat"
srFill(ap_rfiltered_map)[mm_ptids] <- "red"
srFill(ap_rfiltered_map)[mn_ptids] <- "lawngreen"
srFill(ap_rfiltered_map)[mp_ptids] <- "navy"
srFill(ap_rfiltered_map)[pj_ptids] <- "yellow"
srFill(ap_rfiltered_map)[pm_ptids] <- "plum3"
srFill(ap_rfiltered_map)[pm211_ptids] <- "goldenrod"
srFill(ap_rfiltered_map)[pn_ptids] <- "turquoise1"
srFill(ap_rfiltered_map)[pp_ptids] <- "blue"

ptDrawingOrder(ap_rfiltered_map) <- c(
  seq_len(numAntigens(ap_rfiltered_map)) + numSera(ap_rfiltered_map),
  which(mod211_ptids),
  which(pfizer_ptids),
  which(janssen_ptids),
  which(nvx_ptids),
  which(mod_ptids)
)

par(mar=rep(0,4))
plot(ap_rfiltered_map)
legend(x="topright",text.font=1,
       legend=c("Jan/Jan", "Jan/Mod", "Jan/Mod2.11", "Jan/Nvx", "Jan/Pfi",
                "Mod/Jan", "Mod/Mod", "Mod/Nvx", "Mod/Pfi",
                "Pfi/Jan", "Pfi/Mod", "Pfi/Mod2.11", "Pfi/Nvx", "Pfi/Pfi"),
       fill=c("skyblue", "darkgreen", "maroon4", "orange", "hotpink", "wheat",
              "red", "lawngreen", "navy", "yellow", "plum3", "goldenrod",
              "turquoise1", "blue"))

pa_rfiltered_map <- acmap(titer_table = titer_pa_rfiltered_df)
pa_rfiltered_map <- optimizeMap(
  map                     = pa_rfiltered_map,
  number_of_dimensions    = 2,
  number_of_optimizations = 500,
  minimum_column_basis    = "none"
)

mod_ptids <- agNames(pa_rfiltered_map) %in% mod_vt
mod211_ptids <- agNames(pa_rfiltered_map) %in% mod211_vt
pfizer_ptids <- agNames(pa_rfiltered_map) %in% pfizer_vt
janssen_ptids <- agNames(pa_rfiltered_map) %in% janssen_vt
nvx_ptids <- agNames(pa_rfiltered_map) %in% nvx_vt
agFill(pa_rfiltered_map)[mod_ptids] <- "red"
agFill(pa_rfiltered_map)[mod211_ptids] <- "purple"
agFill(pa_rfiltered_map)[pfizer_ptids] <- "blue"
agFill(pa_rfiltered_map)[janssen_ptids] <- "skyblue"
agFill(pa_rfiltered_map)[nvx_ptids] <- "yellow"

ptDrawingOrder(pa_rfiltered_map) <- c(
  seq_len(numSera(pa_rfiltered_map)) + numAntigens(pa_rfiltered_map),
  which(mod211_ptids),
  which(pfizer_ptids),
  which(janssen_ptids),
  which(nvx_ptids),
  which(mod_ptids)
)

plot(pa_rfiltered_map)

ap_cfiltered_map <- acmap(titer_table = titer_ap_cfiltered_df)
ap_cfiltered_map <- optimizeMap(
  map                     = ap_cfiltered_map,
  number_of_dimensions    = 2,
  number_of_optimizations = 1500,
  minimum_column_basis    = "none"
)

jj_ptids <- srNames(ap_cfiltered_map) %in% jan_jan_vt
jm_ptids <- srNames(ap_cfiltered_map) %in% jan_mod_vt
jm211_ptids <- srNames(ap_cfiltered_map) %in% jan_mod211_vt
jn_ptids <- srNames(ap_cfiltered_map) %in% jan_nvx_vt
jp_ptids <- srNames(ap_cfiltered_map) %in% jan_pfi_vt
mj_ptids <- srNames(ap_cfiltered_map) %in% mod_jan_vt
mm_ptids <- srNames(ap_cfiltered_map) %in% mod_mod_vt
mn_ptids <- srNames(ap_cfiltered_map) %in% mod_nvx_vt
mp_ptids <- srNames(ap_cfiltered_map) %in% mod_pfi_vt
pj_ptids <- srNames(ap_cfiltered_map) %in% pfi_jan_vt
pm_ptids <- srNames(ap_cfiltered_map) %in% pfi_mod_vt
pm211_ptids <- srNames(ap_cfiltered_map) %in% pfi_mod211_vt
pn_ptids <- srNames(ap_cfiltered_map) %in% pfi_nvx_vt
pp_ptids <- srNames(ap_cfiltered_map) %in% pfi_pfi_vt
# mod_ptids <- srNames(ap_cfiltered_map) %in% mod_vt
# mod211_ptids <- srNames(ap_cfiltered_map) %in% mod211_vt
# pfizer_ptids <- srNames(ap_cfiltered_map) %in% pfizer_vt
# janssen_ptids <- srNames(ap_cfiltered_map) %in% janssen_vt
# nvx_ptids <- srNames(ap_cfiltered_map) %in% nvx_vt
# omicron_antigen <- agNames(ap_cfiltered_map) %in% omicron_vt
srFill(ap_cfiltered_map)[jj_ptids] <- "skyblue"
srFill(ap_cfiltered_map)[jm_ptids] <- "darkgreen"
srFill(ap_cfiltered_map)[jj_ptids] <- "hotpink"
srFill(ap_cfiltered_map)[mj_ptids] <- "wheat"
srFill(ap_cfiltered_map)[mm_ptids] <- "red"
srFill(ap_cfiltered_map)[mp_ptids] <- "navy"
srFill(ap_cfiltered_map)[pj_ptids] <- "yellow"
srFill(ap_cfiltered_map)[pm_ptids] <- "plum3"
srFill(ap_cfiltered_map)[pp_ptids] <- "blue"
# agFill(ap_cfiltered_map)[omicron_antigen] <- "black"

ptDrawingOrder(ap_cfiltered_map) <- c(
  seq_len(numAntigens(ap_cfiltered_map)) + numSera(ap_cfiltered_map),
  which(mod211_ptids),
  which(pfizer_ptids),
  which(janssen_ptids),
  which(nvx_ptids),
  which(mod_ptids)
)

# setLegend is only useful with view(map), not plot(map)
legend_map <- setLegend(ap_cfiltered_map,
                        c("Moderna Boost", "Pfizer Boost", "Janssen Boost"),
                        fill=c("red", "blue", "skyblue"))

par(mar=rep(0,4))
plot(ap_cfiltered_map)
legend(x="topright",text.font=1, legend=c("Moderna Boost", "Pfizer Boost", "Janssen Boost"),
       fill=c("red","blue","skyblue"))

pa_cfiltered_map <- acmap(titer_table = titer_pa_cfiltered_df)
pa_cfiltered_map <- optimizeMap(
  map                     = pa_cfiltered_map,
  number_of_dimensions    = 2,
  number_of_optimizations = 500,
  minimum_column_basis    = "none"
)

mod_ptids <- agNames(pa_cfiltered_map) %in% mod_vt
mod211_ptids <- agNames(pa_cfiltered_map) %in% mod211_vt
pfizer_ptids <- agNames(pa_cfiltered_map) %in% pfizer_vt
janssen_ptids <- agNames(pa_cfiltered_map) %in% janssen_vt
nvx_ptids <- agNames(pa_cfiltered_map) %in% nvx_vt
agFill(pa_cfiltered_map)[mod_ptids] <- "red"
agFill(pa_cfiltered_map)[mod211_ptids] <- "purple"
agFill(pa_cfiltered_map)[pfizer_ptids] <- "blue"
agFill(pa_cfiltered_map)[janssen_ptids] <- "skyblue"
agFill(pa_cfiltered_map)[nvx_ptids] <- "yellow"

ptDrawingOrder(pa_cfiltered_map) <- c(
  seq_len(numSera(pa_cfiltered_map)) + numAntigens(pa_cfiltered_map),
  which(mod211_ptids),
  which(pfizer_ptids),
  which(janssen_ptids),
  which(nvx_ptids),
  which(mod_ptids)
)
plot(pa_cfiltered_map)
