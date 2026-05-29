# ---- Powdermill Nature Preserve Community Data ----

library(dplyr)
library(vegan)

# ---- 2015 data ----
pnr2015 <- read.csv("Data/Processed/PNR2015_InvertebrateCommunity.csv")
pnr15 <- pnr2015[pnr2015$Quadrat >= 41, ]   # only keeping the treatments (quadrats) that overlap with 2022 data , i.e., 41-64
pnr15 <- na.omit(pnr15)
pnr15 <- pnr15[, -c(3:5,8)]
pnr15$Treatment <- as.factor(pnr15$Treatment)
levels(pnr15$Treatment)


pnr15 <- pnr15 %>% mutate(Treatment = dplyr::recode(Treatment, 
                                                    "B" = "Windthrow", 
                                                    "BS" = "Salvaged"))
levels(pnr15$Treatment)
str(pnr15)
colSums(is.na(pnr15))
pnr15 <- pnr15 %>%  mutate(across(5:last_col(), ~ as.numeric(as.character(.))))


# ---- 2022 data ----
pnr2022 <- read.csv("Data/Processed/Powdermill_Community_2022.csv")
pnr22<- na.omit(pnr2022)
pnr22 <- pnr22[, -c(2,4,5,8)]
pnr22 <- pnr22 %>% rename(Quadrat = Site_ID)  # for naming consistency
str(pnr22)
colSums(is.na(pnr22))
pnr22$Collembola[is.na(pnr22$Collembola)] <- 0
colSums(is.na(pnr22))

pnr22$Treatment <- as.factor(pnr22$Treatment)
levels(pnr22$Treatment)
pnr22 <- pnr22 %>% mutate(Treatment = dplyr::recode(Treatment, 
                                                    "Forest " = "Forest", 
                                                    "windthrow" = "Windthrow"))
levels(pnr22$Treatment)

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

species_cols_pnr15 <- colnames(pnr15)[!(colnames(pnr15) %in% c("Quadrat", "Treatment", "DateSet" , "DateColl", "Year"))]
species_cols_pnr22 <- colnames(pnr22)[!(colnames(pnr22) %in% c("Quadrat", "Treatment", "DateSet" , "DateColl", "Year"))]


all_species <- union(species_cols_pnr15, species_cols_pnr22)  # all species in both years combined



missing_in_pnr15 <- setdiff(all_species, species_cols_pnr15)
pnr15[missing_in_pnr15] <- 0


missing_in_pnr22 <- setdiff(all_species, species_cols_pnr22)
pnr22[missing_in_pnr22] <- 0



pnr15 <- pnr15[, c("Quadrat", "Treatment", "DateSet" , "DateColl", "Year", all_species)]
pnr22 <- pnr22[, c("Quadrat", "Treatment", "DateSet" , "DateColl", "Year", all_species)]



all_data <- rbind(pnr15, pnr22)


# Extracting Transect info

traploc <- read.csv("Data/Processed/PNR_PitfallTrapLocations.csv")

traploc <- traploc %>% rename(Quadrat = Plot)
traploc <- traploc[ , c(1,4)]




all_data <- left_join(all_data, traploc, by = "Quadrat")


all_data <- all_data %>%  relocate(90, .before = 3)


colSums(is.na(all_data))
all_data$Collembola[is.na(all_data$Collembola)] <- 0
colSums(is.na(all_data))

all_data <- all_data %>%  mutate(across(7:last_col(), ~ as.numeric(as.character(.))))


all_data$Abundance <- rowSums(all_data[,7:90], na.rm=TRUE)
all_data$Richness <- apply(all_data[,7:90]>0,1,sum)
all_data$Diversity <- all_data %>% select(-c(1:6)) %>%  diversity(index = "shannon")
str(all_data)
all_data$Richness <- as.numeric(all_data$Richness)


# Calculating trap duration 
all_data$DateColl <- as.Date(all_data$DateColl, format = "%m/%d/%Y")
all_data$DateSet  <- as.Date(all_data$DateSet, format = "%m/%d/%Y")

all_data$TrapTime <- as.numeric(all_data$DateColl - all_data$DateSet)

all_data <- all_data %>%  relocate(94, .before = 6)


# write.csv(all_data, "Data/Saved/pnr_clean.csv", row.names = FALSE)


# ---- Cleaned and Combined Data (2015 and 2022) ----

all_data <- read.csv("Data/Saved/pnr_clean.csv")





# ---- GLMM ----
mean(all_data$Abundance)
var(all_data$Abundance)



library(glmmTMB)
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










