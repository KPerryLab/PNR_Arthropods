# ---- Powdermill Nature Preserve Community Data ----

library(tidyverse)
library(dplyr)

# ---- 2015 data ----
pnr2015 <- read.csv("Data/Processed/PNR2015_InvertebrateCommunity.csv")
pnr15 <- na.omit(pnr2015)
pnr15 <- pnr15 %>% select(-c(Canopy, Veg, Sorter))
pnr15$Treatment <- as.factor(pnr15$Treatment)
levels(pnr15$Treatment)


pnr15 <- pnr15 %>% mutate(Treatment = dplyr::recode(Treatment, 
                                                    "B" = "Windthrow", 
                                                    "BS" = "Salvaged"))
levels(pnr15$Treatment)
str(pnr15)
colSums(is.na(pnr15))
pnr15 <- pnr15 %>%  mutate(across(6:last_col(), ~ as.numeric(as.character(.))))

# Dicyrtomidae and Katiannidae combined into Symphypleona
pnr15 <- pnr15[, !(names(pnr15) %in% c("Dicyrtomidae", "Katiannidae"))]

# ---- 2022 data ----
pnr2022 <- read.csv("Data/Processed/Powdermill_Community_2022.csv")
pnr22<- na.omit(pnr2022)
pnr22 <- pnr22 %>% select(-c(Site_ID_2,Notes, Sorter))
pnr22 <- pnr22 %>% rename(Quadrat = Site_ID)  # for naming consistency
str(pnr22)
colSums(is.na(pnr22))


pnr22$Treatment <- as.factor(pnr22$Treatment)
levels(pnr22$Treatment)
pnr22 <- pnr22 %>% mutate(Treatment = dplyr::recode(Treatment, 
                                                    "Forest " = "Forest", 
                                                    "windthrow" = "Windthrow"))
levels(pnr22$Treatment)

# Removing Orders
pnr22 <- pnr22 %>% dplyr::select(-Coleoptera, -Collembola, -Protura, -Diplura)

# Combined into Symphypleona
pnr22 <- pnr22 %>% dplyr::select(-Dicyrtomidae, -Katiannidae, -Bourletiellidae, -Sminthuridae)

pnr22 <- pnr22 %>% rename(Ptiliidae = Ptilidae)




# ---- Combining 2015 and 2022 ----
pnr15$Year <- "2015"
pnr22$Year <- "2022"

# species_cols_pnr15 <- colnames(pnr15)[!(colnames(pnr15) %in% c("Quadrat", "Treatment", "Interval", "DateSet" , "DateColl", "Year"))]
# species_cols_pnr22 <- colnames(pnr22)[!(colnames(pnr22) %in% c("Quadrat", "Treatment", "Interval", "DateSet" , "DateColl", "Year"))]
# 
# all_species <- union(species_cols_pnr15, species_cols_pnr22)  # all species in both years combined
# 
# missing_in_pnr15 <- setdiff(all_species, species_cols_pnr15)
# 
# missing_in_pnr22 <- setdiff(all_species, species_cols_pnr22)



pnr22$Odontellidae <- NULL

colnames(pnr15)
colnames(pnr22)

all_data <- rbind(pnr15, pnr22)


# Extracting Transect info

traploc <- read.csv("Data/Processed/PNR_PitfallTrapLocations.csv")

traploc <- traploc %>% rename(Quadrat = Plot)
traploc <- traploc[ , c(1,4)]




all_data <- left_join(all_data, traploc, by = "Quadrat")


all_data <- all_data %>%  relocate(51, .before = 3)
all_data <- all_data %>%  relocate(51, .before = 3)


colSums(is.na(all_data))


all_data <- all_data %>%  mutate(across(8:last_col(), ~ as.numeric(as.character(.))))


# Calculating trap duration 
all_data$DateColl <- as.Date(all_data$DateColl, format = "%m/%d/%Y")
all_data$DateSet  <- as.Date(all_data$DateSet, format = "%m/%d/%Y")

all_data$TrapTime <- as.numeric(all_data$DateColl - all_data$DateSet)

all_data <- all_data %>%  relocate(52, .before = 8)


write.csv(all_data, "Data/Saved/pnr_clean.csv", row.names = FALSE)


# ---- Cleaned and Combined Data (2015 and 2022) ----

all_data <- read.csv("Data/Saved/pnr_clean.csv")

all_data <- all_data %>%  filter(Year != 2022 | Interval %in% 1:6)

str(all_data)


# Checking how Intervals (1-6) in 2015 and 2022 differ in set time 

library(lubridate)

# all_data$StartJulian <- as.POSIXlt(all_data$DateSet)$yday + 1
# all_data$EndJulian   <- as.POSIXlt(all_data$DateColl)$yday + 1
# 
# interval_dates <- all_data %>%
#   group_by(Year, Interval) %>%
#   summarise(
#     StartDate = first(DateSet),
#     EndDate = first(DateColl),
#     StartJulian = first(StartJulian),
#     EndJulian = first(EndJulian),
#     .groups = "drop"
#   )


all_data <- all_data %>% dplyr::select(-c(Interval, DateSet, DateColl))


# ---- Pooled by Interval ----

pooled <- all_data %>%  group_by(Quadrat, Treatment, Year, Transect) %>%
  dplyr::summarise(across(everything(), ~ sum(.x, na.rm = TRUE)),.groups = "drop")

pooled$abundance <- rowSums(pooled[,6:49], na.rm=TRUE)

library(hillR) # just FYI, can also use this package (KIP)
# Hill numbers
# q = 0 (default) for species richness
# q = 1 for get shannon entropy
# q = 2 for inverse Simpson

pooled$hill.richness <- hill_taxa(pooled[,6:49], q = 0, MARGIN = 1)
pooled$hill.shannon <- hill_taxa(pooled[,6:49], q = 1, MARGIN = 1)

str(pooled)

colSums(pooled[,6:49])

# ---- GLMM for pooled data ----
mean(pooled$abundance)
var(pooled$abundance)


summary_table <- pooled %>%
  group_by(Year, Treatment) %>%
  summarise(
    mean   = mean(abundance, na.rm = TRUE),
    median = median(abundance, na.rm = TRUE),
    max   = max(abundance, na.rm = TRUE),
    min   = min(abundance, na.rm = TRUE),
    var    = var(abundance, na.rm = TRUE),
    .groups      = "drop"
  )
summary_table

library(glmmTMB)
library(DHARMa)
library(car)


pooled <- pooled %>% mutate(across(c(Year, Treatment, Transect), as.factor))

lapply(pooled[c("Year", "Treatment", "Transect")], levels)

with(all_data, table(Transect, Treatment, Year))
# for these models, did you run with a poisson distribution first? (KIP)

#### Abundance ----

# With random effect and poisson
model1_RE <- glmmTMB(abundance ~ Treatment * Year + (1|Transect) + offset(log(TrapTime)), family=poisson(), data = pooled)
Anova(model1_RE,  type = "III")
summary(model1_RE)

simulateResiduals(model1_RE, n = 1000, re.form = NULL, plot = TRUE)   # Dispersion test significant - overdispersion 

# with random effect and negative binomial 
model1_RE_nb <- glmmTMB(abundance ~ Treatment * Year + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = pooled)
Anova(model1_RE_nb,  type = "III")
summary(model1_RE_nb)

simulateResiduals(model1_RE_nb, n = 1000, re.form = NULL, plot = TRUE) # Diagnostics ok

AIC(model1_RE, model1_RE_nb) 

VarCorr(model1_RE_nb)

# The among-transect variation seems to be close to 0

# Model without the random effect 

model1 <- glmmTMB(abundance ~ Treatment * Year + offset(log(TrapTime)), family=nbinom2(), data = pooled)

summary(model1)
simulateResiduals(model1, plot = TRUE) # Diagnostics ok

AIC(model1_RE_nb, model1)   # Delta = 2



#### Richness ----

mean(pooled$hill.richness)
var(pooled$hill.richness)

# With random effect and poisson
model2_RE <- glmmTMB(hill.richness ~ Treatment * Year + (1|Transect) + offset(log(TrapTime)), family=poisson(), data = pooled)
Anova(model2_RE,  type = "III")
summary(model2_RE)

simulateResiduals(model2_RE, plot = TRUE)   # Dispersion test significant - underdispersion 

# with random effect and compois
model2_RE_cp <- glmmTMB(hill.richness ~ Treatment * Year + (1|Transect) + offset(log(TrapTime)), family=compois(), data = pooled)
Anova(model2_RE_cp,  type = "III")
summary(model2_RE_cp)

simulateResiduals(model2_RE_cp, plot = TRUE) # Diagnostics ok

AIC(model2_RE, model2_RE_cp) # Delta > 2

VarCorr(model2_RE_cp)


# Model without the random effect 

model2 <- glmmTMB(hill.richness ~ Treatment * Year + offset(log(TrapTime)), family=compois(), data = pooled)

summary(model2)
simulateResiduals(model2, plot = TRUE) # Diagnostics ok

AIC(model2_RE_cp, model2)   # Delta > 2   # lower with random effect


#### Shannon Diversity ----
model3_RE <- glmmTMB(hill.shannon ~ Treatment * Year + (1|Transect), family=gaussian(), data = pooled)
summary(model3_RE)
simulateResiduals(model3_RE, plot = TRUE)   # Diagnostic ok


model3 <- glmmTMB(hill.shannon ~ Treatment * Year, family=gaussian(), data = pooled)
summary(model3)
simulateResiduals(model3, plot = TRUE)   # Diagnostic ok


AIC(model3_RE, model3)   # Delta < 1    # lower with random effect












library(emmeans)
library(multcomp)
library(multcompView)

# model1l <- emmeans(model1, ~ Treatment * Year)
# cld(model1l, Letters = letters)
# emmip(model1, Treatment ~ Year, type = "response")
# 
# 
# emmeans(model1, pairwise ~ Year)
# emmeans(model1, pairwise ~ Treatment)







# ---- Repeated Measures GLMM - Unpooled Data ----

unpooled <- read.csv("Data/Saved/pnr_clean.csv")
str(unpooled)

unpooled <- unpooled %>% mutate(across(c(Quadrat, Treatment, Transect, Year, Interval), as.factor))

lapply(unpooled[c("Year", "Treatment", "Transect", "Interval")], levels)

str(unpooled)

unpooled$abundance <- rowSums(unpooled[,9:52], na.rm=TRUE)

library(hillR) 

unpooled$hill.richness <- hill_taxa(unpooled[,9:52], q = 0, MARGIN = 1)
unpooled$hill.shannon <- hill_taxa(unpooled[,9:52], q = 1, MARGIN = 1)




# Just Quadrat (trap) as a random effect 
mod_abundance <- glmmTMB(abundance ~ Treatment * Year + (1|Quadrat) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)  #poisson does not fit well
summary(mod_abundance)
simulateResiduals(mod_abundance, n = 1000, re.form = NULL,plot = TRUE)


# Quadrat and Interval as random effects
mod_abundance_int <- glmmTMB(abundance ~ Treatment * Year + (1|Quadrat) + (1|Interval) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
summary(mod_abundance_int)
simulateResiduals(mod_abundance_int, n = 1000, re.form = NULL,plot = TRUE)


# Quadrat, Interval, and Transect as random effects
mod_abundance_int_tra <- glmmTMB(abundance ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
summary(mod_abundance_int_tra)
simulateResiduals(mod_abundance_int_tra, n = 1000, re.form = NULL,plot = TRUE)


AIC(mod_abundance, mod_abundance_int,mod_abundance_int_tra)

# AIC(mod_abundance, mod_abundance_int,mod_abundance_int_tra)
# df      AIC
# mod_abundance          8 4201.947
# mod_abundance_int      9 4190.675
# mod_abundance_int_tra 10 4192.675




mean(unpooled$hill.richness)
var(unpooled$hill.richness)
hist(unpooled$hill.richness)


mod_richness <- glmmTMB(hill.richness ~ Treatment * Year + (1|Quadrat) + offset(log(TrapTime)), family=poisson(), data = unpooled)
summary(mod_richness)
simulateResiduals(mod_richness,n = 1000, re.form = NULL,  plot = TRUE)  # Poisson does not fit well , there is underdispersion


# Just Quadrat (trap) as a random effect 
mod_richness_compois <- glmmTMB(hill.richness ~ Treatment * Year + (1|Quadrat) + offset(log(TrapTime)), family=compois(), data = unpooled)
summary(mod_richness_compois)
simulateResiduals(mod_richness_compois,n = 1000, re.form = NULL, plot = TRUE)

# Quadrat and Interval as random effects
mod_richness_compois_int <- glmmTMB(hill.richness ~ Treatment * Year + (1|Quadrat) + (1|Interval) + offset(log(TrapTime)), family=compois(), data = unpooled)
summary(mod_richness_compois_int)
simulateResiduals(mod_richness_compois_int,n = 1000, re.form = NULL, plot = TRUE)

# Quadrat, Interval, and Transect as random effects
mod_richness_compois_int_tra <- glmmTMB(hill.richness ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=compois(), data = unpooled)
summary(mod_richness_compois_int_tra)
simulateResiduals(mod_richness_compois_int_tra,n = 1000, re.form = NULL, plot = TRUE)


AIC(mod_richness_compois, mod_richness_compois_int,mod_richness_compois_int_tra)


# AIC(mod_richness_compois, mod_richness_compois_int,mod_richness_compois_int_tra)
# df      AIC
# mod_richness_compois          8 1573.167
# mod_richness_compois_int      9 1531.978
# mod_richness_compois_int_tra 10 1525.994


model3_RE <- glmmTMB(hill.shannon ~ Treatment * Year + (1|Quadrat), family=gaussian(), data = unpooled)
summary(model3_RE)
simulateResiduals(model3_RE, n = 1000, re.form = NULL, plot = TRUE)   # Diagnostic ok

model3_RE_In <- glmmTMB(hill.shannon ~ Treatment * Year + (1|Quadrat) + (1|Interval), family=gaussian(), data = unpooled)
summary(model3_RE_In)
simulateResiduals(model3_RE_In, n = 1000, re.form = NULL, plot = TRUE)   # Diagnostic ok

model3_RE_In_Tra <- glmmTMB(hill.shannon ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect), family=gaussian(), data = unpooled)
summary(model3_RE_In_Tra)
simulateResiduals(model3_RE_In_Tra, n = 1000, re.form = NULL, plot = TRUE)   # Diagnostic ok

model3 <- glmmTMB(hill.shannon ~ Treatment * Year, family=gaussian(), data = unpooled)
summary(model3)
simulateResiduals(model3, plot = TRUE)   # Diagnostic ok


AIC(model3_RE, model3_RE_In, model3_RE_In_Tra, model3)   
# AIC(model3_RE, model3_RE_In, model3_RE_In_Tra, model3)    # lowest with all random effects
# df      AIC
# model3_RE         8 1238.595
# model3_RE_In      9 1220.315
# model3_RE_In_Tra 10 1218.709
# model3            7 1239.043




Anova(mod_abundance_int, type = "III")
emmeans(mod_abundance_int, pairwise ~ Year)
emmeans(mod_abundance_int, pairwise ~ Treatment)


# ---- NMDS ----

all_data <- read.csv("Data/Saved/pnr_clean.csv")
nmds_df <- pooled

library(vegan)
 community <- nmds_df[, -c(1:5, 50)]
 meta <- nmds_df[, c(2,3)]

 
 nmds1 <- metaMDS(community, 
                distance = "bray",   
                k = 2,               
                trymax = 100)
 
 
 scores_df <- as.data.frame(scores(nmds1, display = "sites"))
 scores_df <- cbind(scores_df, meta)
 
 library(ggplot2)
 
 ggplot(scores_df, aes(NMDS1, NMDS2, color = Treatment, shape = as.factor(Year))) +
   geom_point(size = 3) +
   theme_minimal()

# ---- NMDS 2015 ----
 
 all_2015 <- pooled[which(pooled$Year== "2015"), ]
 
 community15 <- all_2015[, -c(1:5, 50)]
 meta15 <- all_2015[, c(1,2)]
 
 
 nmds15 <- metaMDS(community15, 
                  distance = "bray",   
                  k = 3,               
                  trymax = 100)
 
 
 scores_15 <- as.data.frame(scores(nmds15, display = "sites"))
 scores_15 <- cbind(scores_15, meta15)
 
 library(ggplot2)
 
 ggplot(scores_15, aes(NMDS1, NMDS2, color = Treatment, fill = Treatment)) +
   geom_point(size = 3) +
   stat_ellipse(geom = "polygon", alpha = 0.2, color = NA)+
      theme_minimal() +
   labs(title = "PNR 2015") +
   theme(plot.title = element_text(hjust = 0.5))
 
 ggplot(scores_15, aes(NMDS1, NMDS3, color = Treatment, fill = Treatment)) +
   geom_point(size = 3) +
   stat_ellipse(geom = "polygon", alpha = 0.2, color = NA)+
   theme_minimal() +
   labs(title = "PNR 2015") +
   theme(plot.title = element_text(hjust = 0.5))
 
 adonis15 <- adonis2(community15 ~ Treatment, data = meta15, method = "bray")
 adonis15

# ---- NMDS 2022 ----
 
 all_2022 <- pooled[which(pooled$Year== "2022"), ]
 
 community22 <- all_2022[, -c(1:5, 50)]
 meta22 <- all_2022[, c(1,2)]
 
 
 nmds22 <- metaMDS(community22, 
                   distance = "bray",   
                   k = 3,               
                   trymax = 100)
 
 
 scores_22 <- as.data.frame(scores(nmds22, display = "sites"))
 scores_22 <- cbind(scores_22, meta22)
 

 
 ggplot(scores_22, aes(NMDS1, NMDS2, color = Treatment, fill = Treatment)) +
   geom_point(size = 3) +
   stat_ellipse(geom = "polygon", alpha = 0.2, color = NA)+
   theme_minimal() +
   labs(title = "PNR 2022") +
   theme(plot.title = element_text(hjust = 0.5))
 
 ggplot(scores_22, aes(NMDS1, NMDS3, color = Treatment, fill = Treatment)) +
   geom_point(size = 3) +
   stat_ellipse(geom = "polygon", alpha = 0.2, color = NA)+
   theme_minimal() +
   labs(title = "PNR 2022") +
   theme(plot.title = element_text(hjust = 0.5))
 
 
 
 
 
 adonis22 <- adonis2(community22 ~ Treatment, data = meta22, method = "bray")
 adonis22
 
 
 

# ---- Species summaries ----
 colnames(all_data)
 library(dplyr)

 species_summary <- all_data %>%
   dplyr::select(-Quadrat, -Treatment, -Transect, -TrapTime, -Year) %>%
   dplyr::summarise(across(everything(), ~ mean(. > 0, na.rm = TRUE))) %>%
   tidyr::pivot_longer(cols = everything(), names_to = "Species", values_to = "Proportion") %>%
   arrange(desc(Proportion)) #%>%
   #dplyr::filter(Proportion >= 0.30) 
 # dplyr::pull(Species)
 
 print(species_summary$Species)
 
 species_to_analyze <- species_summary %>%
   filter(Proportion >= 0.3) %>%
   pull(Species)
 
# ---- GLMM for taxa ----

 
 library(emmeans)
 library(multcomp)
 mean(unpooled$Araneae)
 var(unpooled$Araneae)
 
 aran <- glmmTMB(Araneae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(aran)
 simulateResiduals(aran, n = 1000, re.form = NULL, plot = TRUE)
 aranl <- emmeans(aran, ~ Treatment * Year)
 cld(aranl, Letters = letters)
 emmip(aran, Treatment ~ Year)
 
 
 tomo <- glmmTMB(Tomoceridae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(tomo)
 simulateResiduals(tomo, n = 1000, re.form = NULL, plot = TRUE)
 tomol <- emmeans(tomo, ~ Treatment * Year)
 cld(tomol, Letters = letters)
 emmip(tomo, Treatment ~ Year)
 
 form <- glmmTMB(Formicidae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(form)
 simulateResiduals(form, n = 1000, re.form = NULL, plot = TRUE)
 forml <- emmeans(form, ~ Treatment * Year)
 cld(forml, Letters = letters)
 emmip(form, Treatment ~ Year)
 

 symp <- glmmTMB(Symphypleona ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(symp)
 simulateResiduals(symp, n = 1000, re.form = NULL, plot = TRUE)
 sympl <- emmeans(symp, ~ Treatment * Year)
 cld(sympl, Letters = letters)
 emmip(symp, Treatment ~ Year)
 
 
 hypo <- glmmTMB(Hypogastruridae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(hypo)
 simulateResiduals(hypo, n = 1000, re.form = NULL, plot = TRUE)
 hypol <- emmeans(hypo, ~ Treatment * Year)
 cld(hypol, Letters = letters)
 emmip(hypo, Treatment ~ Year)
 
 
 stap <- glmmTMB(Staphylinidae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(stap)
 simulateResiduals(stap, n = 1000, re.form = NULL, plot = TRUE)
 stapl<-emmeans(stap, ~ Treatment * Year)
 cld(stapl, Letters = letters)
 emmip(stap, Treatment ~ Year)
 
 
 isot <- glmmTMB(Isotomidae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(isot)
 simulateResiduals(isot, n = 1000, re.form = NULL, plot = TRUE)
 isotl<-emmeans(isot, ~ Treatment * Year)
 cld(isotl, Letters = letters)
 emmip(isot, Treatment ~ Year)
 
 
 poly <- glmmTMB(Polydesmida ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(poly)
 simulateResiduals(poly, n = 1000, re.form = NULL, plot = TRUE)
 polyl<-emmeans(poly, ~ Treatment * Year)
 cld(polyl, Letters = letters)
 emmip(poly, Treatment ~ Year)
 
 ento <- glmmTMB(Entomobryidae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(ento)
 simulateResiduals(ento, n = 1000, re.form = NULL, plot = TRUE)
 entol<-emmeans(ento, ~ Treatment * Year)
 cld(entol, Letters = letters)
 emmip(ento, Treatment ~ Year)
 
 cara <- glmmTMB(Carabidae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(cara)
 simulateResiduals(cara, n = 1000, re.form = NULL, plot = TRUE)
 caral<-emmeans(cara, ~ Treatment * Year)
 cld(caral, Letters = letters)
 emmip(cara, Treatment ~ Year)
 
 

 snai <- glmmTMB(Snails ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(snai)
 simulateResiduals(snai, n = 1000, re.form = NULL, plot = TRUE)
 snail<-emmeans(snai, ~ Treatment * Year)
 cld(snail, Letters = letters)
 emmip(snai, Treatment ~ Year)
 
 
 nean <- glmmTMB(Neanuridae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(nean)
 simulateResiduals(nean, n = 1000, re.form = NULL, plot = TRUE)
 neanl<-emmeans(nean, ~ Treatment * Year)
 cld(neanl, Letters = letters)
 emmip(nean, Treatment ~ Year)
 
 
 opil <- glmmTMB(Opiliones ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(opil)
 simulateResiduals(opil, n = 1000, re.form = NULL, plot = TRUE)
 opill<-emmeans(opil, ~ Treatment * Year)
 cld(opill, Letters = letters)
 emmip(opil, Treatment ~ Year)
 
 
 niti <- glmmTMB(Nitidulidae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(niti)
 simulateResiduals(niti, n = 1000, re.form = NULL, plot = TRUE)
 nitil<-emmeans(niti, ~ Treatment * Year)
 cld(nitil, Letters = letters)
 emmip(niti, Treatment ~ Year)
 
 
 scol <- glmmTMB(Scolytinae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(scol)
 simulateResiduals(scol, n = 1000, re.form = NULL, plot = TRUE)
 scoll<-emmeans(scol, ~ Treatment * Year)
 cld(scoll, Letters = letters, adjust="tukey")
 emmip(scol, Treatment ~ Year)
 
 
 rhap <- glmmTMB(Rhaphidophoridae ~ Treatment * Year + (1|Quadrat) + (1|Interval) + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = unpooled)
 summary(rhap)
 simulateResiduals(rhap, n = 1000, re.form = NULL, plot = TRUE)
 rhapl<-emmeans(rhap, ~ Treatment * Year)
 cld(rhapl, Letters = letters)
 emmip(rhap, Treatment ~ Year)
 
 
 
 
 
 
 
 
 

 
  
 
 
 
 
 
  