# ---- Powdermill Nature Preserve Community Data ----

library(dplyr)


# ---- 2015 data ----
pnr2015 <- read.csv("Data/Processed/PNR2015_InvertebrateCommunity.csv")
pnr15 <- pnr2015[pnr2015$Quadrat >= 41, ]   # only keeping the treatments (quadrats) that overlap with 2022 data , i.e., 41-64
pnr15 <- na.omit(pnr15)
pnr15 <- pnr15[, -c(3,4,8)]
pnr15$Treatment <- as.factor(pnr15$Treatment)
levels(pnr15$Treatment)


pnr15 <- pnr15 %>% mutate(Treatment = dplyr::recode(Treatment, 
                                                    "B" = "Windthrow", 
                                                    "BS" = "Salvaged"))
levels(pnr15$Treatment)
str(pnr15)
colSums(is.na(pnr15))
pnr15 <- pnr15 %>%  mutate(across(6:last_col(), ~ as.numeric(as.character(.))))


# Combining Dicyrtomidae and Katiannidae as Globular Springtails
pnr15 <- pnr15 %>% 
  mutate(Glob_Springt = Dicyrtomidae + Katiannidae, .keep = "unused")

# ---- 2022 data ----
pnr2022 <- read.csv("Data/Processed/Powdermill_Community_2022.csv")
pnr22<- na.omit(pnr2022)
pnr22 <- pnr22[, -c(2,4,8)]
pnr22 <- pnr22 %>% rename(Quadrat = Site_ID)  # for naming consistency
str(pnr22)
colSums(is.na(pnr22))
# pnr22$Collembola[is.na(pnr22$Collembola)] <- 0
# colSums(is.na(pnr22))

pnr22$Treatment <- as.factor(pnr22$Treatment)
levels(pnr22$Treatment)
pnr22 <- pnr22 %>% mutate(Treatment = dplyr::recode(Treatment, 
                                                    "Forest " = "Forest", 
                                                    "windthrow" = "Windthrow"))
levels(pnr22$Treatment)

pnr22 <- pnr22 %>% select(-Coleoptera, -Collembola, -Protura, -Diplura)

# Combining Dicyrtomidae, Katiannidae, Bourletiellidae and Sminthuridae as Globular Springtails
pnr22 <- pnr22 %>% 
  mutate(Glob_Springt = Dicyrtomidae + Katiannidae + Bourletiellidae + Sminthuridae, .keep = "unused")

pnr22 <- pnr22 %>% rename(Ptiliidae = Ptilidae)

# # Adding Paronellidae into Entomobryidae 
# 
# pnr22 <- pnr22 %>% mutate(Entomobryidae = Entomobryidae + Paronellidae) %>%  select(-Paronellidae)
#  
# 
# # Removing Poduridae
# pnr22 <- pnr22 %>% select(-Poduridae)    
# 



# ---- Combining 2015 and 2022 ----
pnr15$Year <- "2015"
pnr22$Year <- "2022"

species_cols_pnr15 <- colnames(pnr15)[!(colnames(pnr15) %in% c("Quadrat", "Treatment", "Interval", "DateSet" , "DateColl", "Year"))]
species_cols_pnr22 <- colnames(pnr22)[!(colnames(pnr22) %in% c("Quadrat", "Treatment", "Interval", "DateSet" , "DateColl", "Year"))]


all_species <- union(species_cols_pnr15, species_cols_pnr22)  # all species in both years combined



missing_in_pnr15 <- setdiff(all_species, species_cols_pnr15)
pnr15[missing_in_pnr15] <- 0


missing_in_pnr22 <- setdiff(all_species, species_cols_pnr22)
# pnr22[missing_in_pnr22] <- 0



pnr15 <- pnr15[, c("Quadrat", "Treatment", "Interval", "DateSet" , "DateColl", "Year", all_species)]
pnr22 <- pnr22[, c("Quadrat", "Treatment", "Interval", "DateSet" , "DateColl", "Year", all_species)]



all_data <- rbind(pnr15, pnr22)


# Extracting Transect info

traploc <- read.csv("Data/Processed/PNR_PitfallTrapLocations.csv")

traploc <- traploc %>% rename(Quadrat = Plot)
traploc <- traploc[ , c(1,4)]




all_data <- left_join(all_data, traploc, by = "Quadrat")


all_data <- all_data %>%  relocate(53, .before = 3)


colSums(is.na(all_data))
# all_data$Collembola[is.na(all_data$Collembola)] <- 0
# colSums(is.na(all_data))

all_data <- all_data %>%  mutate(across(8:last_col(), ~ as.numeric(as.character(.))))

library(vegan)

all_data$Abundance <- rowSums(all_data[,8:53], na.rm=TRUE)
all_data$Richness <- apply(all_data[,8:53]>0,1,sum)
all_data$Diversity <- all_data %>% select(-c(1:7)) %>%  diversity(index = "shannon")
str(all_data)
all_data$Richness <- as.numeric(all_data$Richness)


# Calculating trap duration 
all_data$DateColl <- as.Date(all_data$DateColl, format = "%m/%d/%Y")
all_data$DateSet  <- as.Date(all_data$DateSet, format = "%m/%d/%Y")

all_data$TrapTime <- as.numeric(all_data$DateColl - all_data$DateSet)

all_data <- all_data %>%  relocate(57, .before = 7)


write.csv(all_data, "Data/Saved/pnr_clean.csv", row.names = FALSE)


# ---- Cleaned and Combined Data (2015 and 2022) ----

all_data <- read.csv("Data/Saved/pnr_clean.csv")





# ---- GLMM Combined----
mean(all_data$Abundance)
var(all_data$Abundance)


summary_table <- all_data %>%
  group_by(Year, Treatment) %>%
  summarise(
    mean   = mean(Abundance, na.rm = TRUE),
    median = median(Abundance, na.rm = TRUE),
    max   = max(Abundance, na.rm = TRUE),
    min   = min(Abundance, na.rm = TRUE),
    var    = var(Abundance, na.rm = TRUE),
    .groups      = "drop"
  )
summary_table

library(glmmTMB)
library(DHARMa)
all_data$Year <- as.factor(all_data$Year)
all_data$Treatment <- as.factor(all_data$Treatment)
all_data$Transect <- as.factor(all_data$Transect)


levels(all_data$Treatment)


model1 <- glmmTMB(Abundance ~ Treatment * Year + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = all_data)

summary(model1)


library(emmeans)


emmeans(model1, ~ Treatment * Year)

mean(all_data$Richness)
var(all_data$Richness)


model2 <- glmmTMB(Richness ~ Treatment * Year + (1|Transect)+ offset(log(TrapTime)), family=poisson(), data = all_data)

summary(model2)

emmeans(model2, ~ Treatment * Year)




mean(all_data$Diversity)
var(all_data$Diversity)
hist(all_data$Diversity)

model3 <- glmmTMB(Diversity ~ Treatment * Year + (1|Transect)+ offset(TrapTime), family=gaussian(), data = all_data)
summary(model3)

emmeans(model3, ~ Treatment * Year)
emmeans(model3, ~ Year)


# ---- NMDS ----

library(vegan)
 community <- all_data[, -c(1:8, 55:57)]
 meta <- all_data[, c(2,8)]

 
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
 
 all_2015 <- all_data[which(all_data$Year== "2015"), ]
 
 community15 <- all_2015[, -c(1:8, 55:57)]
 meta15 <- all_2015[, c(1,2)]
 
 
 nmds15 <- metaMDS(community15, 
                  distance = "bray",   
                  k = 2,               
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
 
 adonis15 <- adonis2(community15 ~ Treatment, data = meta15, method = "bray")
 adonis15

# ---- NMDS 2022 ----
 
 all_2022 <- all_data[which(all_data$Year== "2022"), ]
 
 community22 <- all_2022[, -c(1:8, 55:57)]
 meta22 <- all_2022[, c(1,2)]
 
 
 nmds22 <- metaMDS(community22, 
                   distance = "bray",   
                   k = 2,               
                   trymax = 100)
 
 
 scores_22 <- as.data.frame(scores(nmds22, display = "sites"))
 scores_22 <- cbind(scores_22, meta22)
 

 
 ggplot(scores_22, aes(NMDS1, NMDS2, color = Treatment, fill = Treatment)) +
   geom_point(size = 3) +
   stat_ellipse(geom = "polygon", alpha = 0.2, color = NA)+
   theme_minimal() +
   labs(title = "PNR 2022") +
   theme(plot.title = element_text(hjust = 0.5))
 
 
 adonis22 <- adonis2(community22 ~ Treatment, data = meta22, method = "bray")
 adonis22
 
 
 
# ---- GLMM for 2015 ----
 model1 <- glmmTMB(Abundance ~ Treatment + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = all_2015)
 
 summary(model1)
 
 emmeans(model1, ~Treatment)
 
 model2 <- glmmTMB(Richness ~ Treatment + (1|Transect) + offset(log(TrapTime)), family=poisson(), data = all_2015)
 
 summary(model2)
 emmeans(model2, ~Treatment)
 
 model3 <- glmmTMB(Diversity ~ Treatment + (1|Transect) + offset(TrapTime), family=gaussian(), data = all_2015)
 
 summary(model3)
 emmeans(model3, ~Treatment)
 
# ---- GLMM for 2022 ----
 model1 <- glmmTMB(Abundance ~ Treatment + (1|Transect) + offset(log(TrapTime)), family=nbinom2(), data = all_2022)
 
 summary(model1)
 
 emmeans(model1, ~Treatment)
 
 model2 <- glmmTMB(Richness ~ Treatment + (1|Transect) + offset(log(TrapTime)), family=poisson(), data = all_2022)
 
 summary(model2)
 emmeans(model2, ~Treatment)
 
 model3 <- glmmTMB(Diversity ~ Treatment + (1|Transect) + offset(TrapTime), family=gaussian(), data = all_2022)
 
 summary(model3)
 emmeans(model3, ~Treatment)
 
# ---- Species summaries ----
 colnames(all_data)
 library(dplyr)
 library(multcomp)
 species_summary <- all_data %>%
   dplyr::select(-Quadrat, -Treatment, -Transect, -Interval, -DateSet, -DateColl, -TrapTime, -Year, -Abundance, -Richness, -Diversity) %>%
   dplyr::summarise(across(everything(), ~ mean(. > 0, na.rm = TRUE))) %>%
   tidyr::pivot_longer(cols = everything(), names_to = "Species", values_to = "Proportion") %>%
   arrange(desc(Proportion)) #%>%
   #dplyr::filter(Proportion >= 0.30) 
 # dplyr::pull(Species)
 
 print(species_summary$Species)
 
# ---- GLMM for taxa ----
 mean(all_data$Araneae)
 var(all_data$Araneae)
 
 aran <- glmmTMB(Araneae ~ Treatment * Year + (1|Transect)+ offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(aran)
 emmeans(aran, ~ Treatment * Year)
 
 tomo <- glmmTMB(Tomoceridae ~ Treatment * Year + (1|Transect)+  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(tomo)
 emmeans(tomo, ~ Treatment * Year)

 dicy <- glmmTMB(Dicyrtomidae ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(dicy)
 emmeans(dicy, ~ Treatment * Year)
 
 
 hypo <- glmmTMB(Hypogastruridae ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(hypo)
 emmeans(hypo, ~ Treatment * Year)
 
 
 stap <- glmmTMB(Staphylinidae ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(stap)
 stapl<-emmeans(stap, ~ Treatment * Year)
 cld(stapl, Letters = letters, adjust = "tukey")
 
 poly <- glmmTMB(Polydesmida ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(poly)
 polyl<-emmeans(poly, ~ Treatment * Year)
 cld(polyl, Letters = letters, adjust = "tukey")
 
 ento <- glmmTMB(Entomobryidae ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(ento)
 entol<-emmeans(ento, ~ Treatment * Year)
 cld(entol, Letters = letters, adjust = "tukey")
 
 
 cara <- glmmTMB(Carabidae ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(cara)
 caral<-emmeans(cara, ~ Treatment * Year)
 cld(caral, Letters = letters, adjust = "tukey")

 coll <- glmmTMB(Collembola ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(coll)
 colll<-emmeans(coll, ~ Treatment * Year)
 cld(colll, Letters = letters, adjust = "tukey")
 
 
 coll <- glmmTMB(Collembola ~ Treatment * Year +  offset(log(TrapTime)), family=nbinom2(), data = all_data)
 summary(coll)
 colll<-emmeans(coll, ~ Treatment * Year)
 cld(colll, Letters = letters, adjust = "tukey")
 
 
 
 
 

 
  
 
 
 
 
 
  