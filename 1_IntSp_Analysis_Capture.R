#Code to test for differences in length and condition at capture for JSOES
#Interior Spring Chinook yearlings from 2015-2019 and 2021
#Updated on 10/16/2025


# Vector of required packages
pkgs <- c(
  "rstudioapi", "tidyverse", "ggplot2", "dplyr", "ggpubr", "MuMIn", 
  "emmeans", "car", "nlme", "lme4", "ggeffects", 
   "performance", "broom", "multcompView" 
)

# Install any missing packages
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(pkgs[!installed])
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #

#####-------------Import data for 2015-2021 capture comparison------#######

trawl_2015_2021 <- read.csv("1_TrawlIDs_May2015-2021.csv") %>% 
  #filter(Year >2014) %>% 
  #filter(Year < 2022) %>%
  mutate(origin = case_when(Tag == 'YES' ~ PBTHatcheryAssignment, TRUE ~ "Wild")) %>% 
  mutate(origin= ifelse(origin== "Chief_Joseph_sp", "Chief Joseph", origin)) %>% 
  mutate(origin= ifelse(origin== "Little_White_Salmon_sp", "Little White Salmon", origin)) %>% 
  filter(!is.na(origin)) %>% 
  filter(!is.na(FultonCondition)) %>% 
  filter(FINAL_GSI == "Interior_Sp")

morethanfive <-trawl_2015_2021  %>%
  mutate(origin.f = as.factor(origin)) %>% 
  mutate(year.f = as.factor(Year)) %>%
  mutate(origin.f = na_if(origin.f, "")) %>% 
  dplyr::group_by(origin.f, year.f ) %>%  # Group by 'origin' and 'Year'
  filter(n() >= 5) %>%
  ungroup() %>% 
  mutate(year.f = droplevels(year.f)) %>% 
  mutate(origin.f = droplevels(origin.f)) 


summarybyorigin <- morethanfive%>%
  filter(!is.na(origin.f)) %>% 
  dplyr::group_by(origin.f, year.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group


summarybyoriginonly <- morethanfive%>%
  filter(!is.na(origin.f)) %>% 
  mutate(origin.f = as.factor(origin)) %>% 
  dplyr::group_by(origin.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group


unique(morethanfive$origin.f)

morethanfive  <- morethanfive %>% 
  filter(!is.na(origin.f)) %>% 
  filter(!origin.f %in% c("Klickitat", "Methow", "Warm Springs", "SF_Walla_Walla"))
         
unique(morethanfive$origin.f)

summarybyoriginonly <- morethanfive%>%
  filter(!is.na(origin.f)) %>% 
  mutate(origin.f = as.factor(origin)) %>% 
  dplyr::group_by(origin.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group


trawl_2015_2021 <- morethanfive

summary1 <- trawl_2015_2021 %>%
  dplyr::group_by(origin, Year) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), mean_length = mean(FieldLength_mm), mean_mass = mean(fish_mass), na.rm = TRUE, .groups = "drop")  # Count rows for each group


####----Test effect of Origin on size at capture----#####  

sizemodel <- lm(log(FieldLength_mm) ~ origin, data = trawl_2015_2021)
sizemodel.r <- lmer(log(FieldLength_mm) ~ origin + (1|year.f), data = trawl_2015_2021)

AIC(sizemodel.r, sizemodel)
summary(sizemodel)
check_model(sizemodel)
Anova(sizemodel)

mydf1 <- ggpredict(sizemodel)
plot(mydf1)

tidy_output <- tidy(sizemodel)
print(tidy_output)
write.csv(tidy_output, "lmm_results.csv", row.names = FALSE)

a <- TukeyHSD(aov(sizemodel))
a

# Extract the p-values from the "origin" term in Tukey output
pvals <- a$origin[, "p adj"]

# Make sure names are in format "group1-group2"
names(pvals) <- rownames(a$origin)

#pairwise comparison letter assignments 
letters <- multcompLetters(pvals, threshold = 0.05)
print(letters$Letters)

write.csv(mydf1, "1_predsizeatcap.csv")

levels(trawl_2015_2021$year.f)

plot(mydf1)+
  theme_classic(base_size = 19)+
  theme(legend.position = "top")+
  theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "", y = "Predicted size at capture, mm")+
  ggtitle(NULL)

capsize <- read.csv("1_predsizeatcap.csv")

capsize <- capsize %>%
  mutate(origin.group = case_when(
  origin.x == "Carson" ~ "Mid",
origin.x == "Chief Joseph" ~ "Upper",
origin.x == "Clearwater" ~ "Snake",
origin.x == "Dworshak" ~ "Snake",
origin.x == "Lookingglass" ~ "Snake",
origin.x == "Leavenworth" ~ "Upper",
origin.x == "Pahsimeroi" ~ "Snake",
origin.x == "Little White Salmon" ~ "Mid",
origin.x == "McCall" ~ "Snake",
origin.x == "Rapid River" ~ "Snake",
origin.x == "Sawtooth" ~ "Snake",
origin.x == "Walla Walla" ~ "Mid",
origin.x == "Wild" ~ "Wild",
origin.x == "Winthrop" ~ "Upper",
origin.x == "Parkdale" ~ "Mid"
))

hatch_colors <- c("Carson" = "#481D6FFF", "Walla Walla" = "#481D6FFF", "Parkdale" = "#481D6FFF", 
                "Chief Joseph" = "#277F8EFF","Winthrop" = "#277F8EFF","Leavenworth" = "#277F8EFF",  "Clearwater" = "#71CF57FF", "Dworshak" = "#71CF57FF", "Little White Salmon" = "#481D6FFF", "Rapid River" = "#71CF57FF", "Lookingglass" = "#71CF57FF", "Sawtooth" = "#71CF57FF", "Wild" = "#FDE725FF")

capsize$origin.x <- factor(capsize$origin.x,
                         levels = c("Carson", "Little White Salmon", "Parkdale", "Chief Joseph", "Winthrop", "Leavenworth", "Clearwater",
                                    "Dworshak", "Lookingglass" , "McCall", "Pahsimeroi", "Rapid River", "Sawtooth",  "Wild"))

# Create a custom color scale
hatch_colors_with_alpha <- scales::alpha(hatch_colors, alpha = 0.6)

hatch_color_scale <- scale_fill_manual(values = hatch_colors_with_alpha)
# Assign colors using the viridis color palette for 3 basin + wild
basin_colors <- c("Mid" = "#481D6FFF", "Upper" = "#277F8EFF", "Snake" = "#71CF57FF", "Wild" = "#FDE725FF")

# Create a custom color scale
basin_colors_with_alpha <- scales::alpha(basin_colors, alpha = 0.6) 

basin_color_scale <- scale_fill_manual(values=basin_colors_with_alpha)

# Assign colors using the viridis color palette for 3 basin + wild
basin_colors <- c("Mid" = "#481D6FFF", "Upper" = "#277F8EFF", "Snake" = "#71CF57FF", "Wild" = "#FDE725FF")

capsize$origin.predicted <- as.numeric(as.character(capsize$origin.predicted))
capsize$origin.std.error <- as.numeric(as.character(capsize$origin.std.error))

# Predicted size at capture 
  
cap_plot <- ggplot(capsize, aes(x = origin.x, y = origin.predicted, color = origin.group)) +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = origin.conf.low, ymax = origin.conf.high),
                width = 0.4, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = basin_colors) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", 
        text = element_text(size = 24),
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust x-axis text rotation angle
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 30, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold")) +
  labs(x = NULL, y = "Predicted size at capture, mm", color = "Basin")
plot(cap_plot)

ggsave("Predsizeatcap_new.png", 
       dpi = 300, width = 15, height = 10, units = "in", bg = "white")

####----Test effect of Origin on condition at capture----#####  

condmodel.nr <- lm(FultonCondition ~ origin, data = trawl_2015_2021)
condmodel <- lmer(FultonCondition~ origin + (1|year.f), data = trawl_2015_2021)

AIC(condmodel.nr, condmodel)

condmodel <- lmer(FultonCondition~ origin + (1|year.f), REML = TRUE, data = trawl_2015_2021)

summary(condmodel)
r.squaredGLMM(condmodel)
check_model(condmodel)
Anova(condmodel)

mydf2 <- ggpredict(condmodel)
plot(mydf2)

tidy_output <- tidy(condmodel)
print(tidy_output)
write.csv(tidy_output, "lmm_results.csv", row.names = FALSE)

emmeans_cond <- emmeans(condmodel, pairwise ~ origin)
summary(emmeans_cond)

emmeans_cond <- emmeans(condmodel, pairwise ~ origin)

# 3. Extract the summary of contrasts (pairwise differences)
pairwise_summary <- summary(emmeans_cond$contrasts)

# 4. Pull out p-values and assign comparison names
pvals <- pairwise_summary$p.value

# Sort names in each contrast alphabetically
pairwise_summary$contrast <- sapply(strsplit(as.character(pairwise_summary$contrast), " - "), function(x) paste(sort(x), collapse = "-"))

# Assign cleaned names
names(pvals) <- pairwise_summary$contrast


# 5. Use multcompLetters to get groupings
letters <- multcompLetters(pvals, threshold = 0.05)

# 6. Print the group letters
print(letters)

summary(condmodel)

Anova(condmodel, type = 2)
plot(residuals(condmodel))
check_model(condmodel)
r.squaredGLMM(condmodel)

mydf2

plot(mydf2)+
  theme_classic(base_size = 19)+
  theme(legend.position = "top")+
  theme(legend.text = element_text(size = 14))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "", y = "Predicted condition at capture")+
  ggtitle(NULL)

write.csv(mydf2, "1_predcapcondidion.csv")

capcondition <- read.csv("1_predcapcondidion.csv")

capcondition <- capcondition %>%
  mutate(origin.group = case_when(
    origin.x == "Carson" ~ "Mid",
    origin.x == "Chief Joseph" ~ "Upper",
    origin.x == "Clearwater" ~ "Snake",
    origin.x == "Dworshak" ~ "Snake",
    origin.x == "Lookingglass" ~ "Snake",
    origin.x == "Leavenworth" ~ "Upper",
    origin.x == "Pahsimeroi" ~ "Snake",
    origin.x == "Little White Salmon" ~ "Mid",
    origin.x == "McCall" ~ "Snake",
    origin.x == "Rapid River" ~ "Snake",
    origin.x == "Sawtooth" ~ "Snake",
    origin.x == "Walla Walla" ~ "Mid",
    origin.x == "Wild" ~ "Wild",
    origin.x == "Winthrop" ~ "Upper",
    origin.x == "Parkdale" ~ "Mid"
  ))


capcondition$origin.predicted <- as.numeric(as.character(capcondition$origin.predicted))
capcondition$origin.std.error <- as.numeric(as.character(capcondition$origin.std.error))

capcondition$origin.x <- factor(capcondition$origin.x, levels = c("Carson", "Little White Salmon", "Parkdale",
                                                                  "Chief Joseph", "Winthrop", "Leavenworth", "Clearwater",
                                                    "Dworshak", "Lookingglass", "McCall", "Pahsimeroi", "Rapid River", "Sawtooth",  "Wild"))
# Predicted condition at capture

cap_cond <- ggplot(capcondition, aes(x = origin.x, y = origin.predicted, color = origin.group)) +
  geom_point(size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = origin.conf.low, ymax = origin.conf.high),
                width = 0.4, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = basin_colors) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", 
        text = element_text(size = 24),
        legend.text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust x-axis text rotation angle
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 30, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold")) +
  labs(x = NULL, y = "Predicted condition at capture, K", color = "Basin")

plot(cap_cond)

ggsave("Predcondatcap_new.png", 
       dpi = 300, width = 15, height = 10, units = "in", bg = "white")

######