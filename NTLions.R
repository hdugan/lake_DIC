library(tidyverse)
source('lake_DIC.R') # bicarbonate function

# Lake locations
lakelocations = data.frame(lakeid = c('AL','BM','CB','CL','SP','TB','TL','ME','MO','WI','FI'),
                           Lake = c("Allequash Lake", "Big Muskellunge Lake", 
                                    "Crystal Bog", "Crystal Lake", "Sparkling Lake", "Trout Bog", 
                                    "Trout Lake", "Lake Mendota", "Lake Monona", "Lake Wingra", "Fish Lake"),
                           Lat = c(46.038317, 46.021067, 46.007583, 46.00275, 46.007733, 
                                   46.04125, 46.029267, 43.09885, 43.06337, 43.05258, 43.28733), 
                           Long = c(89.620617, -89.611783, -89.606183, -89.612233, -89.701183, 
                                    -89.686283, -89.665017, -89.40545, -89.36086, -89.42499, 
                                    -89.65173)) 


###################### DOWNLOAD NTL-LTER DATA FROM EDI ######################
# doi:10.6073/pasta/a457e305538a0d8e669b58bb6f35721f
# Package ID: knb-lter-ntl.2.34 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Chemical Limnology of Primary Study Lakes: Major Ions 1981 - current.
infile1 <- tempfile()
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/2/34/3f740d0b77b3caf6930a8ce9cca4306a" 
download.file(inUrl1, infile1, method="curl")
ions <- read_csv(infile1)


# doi:10.6073/pasta/8359d27bbd91028f222d923a7936077d
# Package ID: knb-lter-ntl.1.52 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Chemical Limnology of Primary Study Lakes: Nutrients, pH and Carbon 1981 - current.
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/1/52/802d63a4c35050b09ef6d1e7da3efd3f" 
infile1 <- tempfile()
download.file(inUrl1, infile1, method="curl")
nutrients <- read_csv(infile1)


# Package ID: knb-lter-ntl.29.8 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER:
# Physical Limnology of Primary Study Lakes 1981 - current
inUrl3  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/29/28/03e232a1b362900e0f059859abe8eb97"
infile1 <- tempfile()
download.file(inUrl3, infile1, method="curl")
temp <- read_csv(infile1)

#############################################################################

# Clean data by removing flagged data 
LTERtemp = 
  temp %>% select(lakeid, sampledate, depth, rep, wtemp) |> 
  group_by(lakeid, sampledate, depth) |> 
  summarise(temp = mean(wtemp, na.rm = T))

ions.long = ions %>% select(lakeid, sampledate, depth, rep, cl:k, cond, flagcl:flagk, flagcond) |> 
  mutate(across(everything(), ~replace(., .<0 , NA))) %>%
  rename_at(vars(cl:cond), ~ str_c("value_",.)) %>%
  rename_at(vars(flagcl:flagcond), ~ str_c("error_",.)) %>%
  rename_all(~str_replace_all(.,"flag","")) %>%
  pivot_longer(-(lakeid:rep), names_to = c('.value','item'), names_sep = '_') %>%
  filter(!is.na(value) & value>= 0) %>%
  filter(!str_detect(error,'A|K|L|H') | is.na(error)) %>%
  dplyr::select(-error)

ions.wide = ions.long |> pivot_wider(names_from = item, values_from = value, values_fn = mean) |> 
  select(-cond, cond) |> 
  filter(if_all(ca:so4, ~ !is.na(.))) |> 
  group_by(lakeid, sampledate, depth) |> 
  summarise(across(ca:cond, mean, na.rm = T))

ph.long = nutrients |> select(lakeid, sampledate, depth, rep, ph, dic, alk, flagph, flagdic, flagalk) |> 
  mutate(across(everything(), ~replace(., .<0 , NA))) %>%
  rename_at(vars(ph:alk), ~ str_c("value_",.)) %>%
  rename_at(vars(flagph:flagalk), ~ str_c("error_",.)) %>%
  rename_all(~str_replace_all(.,"flag","")) %>%
  pivot_longer(-(lakeid:rep), names_to = c('.value','item'), names_sep = '_') %>%
  filter(!is.na(value) & value>= 0) %>%
  filter(!str_detect(error,'A|K|L|H') | is.na(error)) %>%
  dplyr::select(-error)

ph.wide = ph.long |> pivot_wider(names_from = item, values_from = value, values_fn = mean) |> 
  filter(if_all(ph:dic, ~ !is.na(.)))  |> 
  group_by(lakeid, sampledate, depth) |> 
  summarise(across(ph:dic, mean, na.rm = T))

# Join dataframes together
df = ions.wide |> inner_join(ph.wide) |> 
  inner_join(LTERtemp)

################################# Calculate bicarbonate #################################
df2 = df %>% 
  rowwise() |>
  mutate(bicarbonate = carbonate(TEMP = temp, DIC = dic, pH = ph, output = 'mg') |> pull(bicarbonate_mgkg)) |> 
  mutate(carbonate = carbonate(TEMP = temp, DIC = dic, pH = ph, output = 'mg') |> pull(carbonate_mgkg)) |> 
  mutate(alkalinity_calculated = carbonate(TEMP = temp, DIC = dic, pH = ph, output = 'mg') |> pull(alkalinity_uEqkg))

# Plot check on calculated alkalinity  
ggplot(df2) +
  geom_point(aes(x = alk, y = alkalinity_calculated)) +
  geom_abline() +
  facet_wrap(~lakeid, scales = 'free')

# Check ion balance
df.meq = df2 |> 
  mutate(hplus = 1000 * 10^(-ph)) |> #millimoles of H3O+
  mutate(ca = 2 * ca/40.078) |> 
  mutate(mg = 2 * mg/24.305) |> 
  mutate(na = 1 * na/22.989) |> 
  mutate(k = 1 * k/39.098) |> 
  mutate(cl = 1 * cl/35.453) |> 
  mutate(so4 = 2 * so4/96.06) |> 
  mutate(bicarbonate = 1 * bicarbonate/61.0168) |> 
  mutate(carbonate = 2 * carbonate/60.009) |> 
  mutate(cations =  ca + mg + na + k) |> 
  mutate(anions = cl + so4 + bicarbonate + carbonate) 

# Plot check on ion balance
ggplot(df.meq) +
  geom_point(aes(x = cations, y = anions)) +
  geom_abline() +
  facet_wrap(~lakeid, scales = 'free')

################################# Output datasets #################################
output.mg = df2 |> left_join(lakelocations)
write_csv(output.mg, 'NTLions_mg.csv')

output.meq = df.meq |> left_join(lakelocations)
write_csv(output.meq, 'NTLions_meq.csv')
