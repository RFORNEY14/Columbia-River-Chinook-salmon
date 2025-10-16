#Code to test for effects of length, mass, and condition at capture and origin on smolt-to-adult
#return ratios for JSOES Interior Spring Chinook yearlings from 2015-2019 and 2021
#Data downloaded from Fish Passage Center (fpc.org).
#MidColumbia = BON to BON dam, Snake = LGD to BON, corrected with hatchery-specific instream survival estimates and percent transported from LGD to BON so effectively
#BON-BON, Upper Columbia = MCN to BON, corrected with annual MCN-BON survival from NOAA
#Updated on 10/15/2025 


# Vector of required packages
pkgs <- c(
  "rstudioapi", "tidyverse", "ggplot2", "dplyr", "ggpubr", "MuMIn", "lmerTest", 
  "emmeans", "datawizard", "sjPlot", "car", "nlme", "lme4", "ggeffects", 
  "LaplacesDemon", "performance", "sjPlot", "boot"
)

# Install any missing packages
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(pkgs[!installed])
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Import SARs

SAR <- read.csv("4_SAR_match.csv") %>% 
  mutate(origin.f = as.factor(unique_origin))


#Import information on catch

trawl_2015_2021 <- read.csv("1_TrawlIDs_May2015-2021.csv") %>% 
  filter(!is.na(PBTHatcheryAssignment)) %>% 
  filter(!is.na(FultonCondition)) %>% 
  mutate(origin = case_when(Tag == 'YES' ~ PBTHatcheryAssignment, TRUE ~ "Wild")) %>% 
  mutate(origin= ifelse(origin== "Chief_Joseph_sp", "Chief Joseph", origin)) %>% 
  mutate(origin= ifelse(origin== "Little_White_Salmon_sp", "Little White Salmon", origin)) %>% 
  filter(FINAL_GSI == "Interior_Sp")

#Select hatcheries

trawl_2015_2021 <- trawl_2015_2021 %>% 
  mutate(PBTHatcheryAssignment = as.factor(PBTHatcheryAssignment)) %>% 
  filter(PBTHatcheryAssignment == "Carson" | 
           PBTHatcheryAssignment == "Dworshak" |
           PBTHatcheryAssignment == "Clearwater"|
           PBTHatcheryAssignment == "Leavenworth" | 
           PBTHatcheryAssignment == "McCall"|
           PBTHatcheryAssignment == "Pahsimeroi" |
           PBTHatcheryAssignment == "Rapid River" |
           PBTHatcheryAssignment == "Sawtooth" |
           PBTHatcheryAssignment == "Winthrop" &
           PBTHatcheryAssignment != "Warm Springs") %>%
  mutate(PBTHatcheryAssignment = droplevels(PBTHatcheryAssignment)) %>% 
  mutate(year.f = as.factor(Year)) %>% 
  mutate(origin.f = as.factor(PBTHatcheryAssignment)) 

unique(trawl_2015_2021$PBTHatcheryAssignment)


####----PBTxSAR----##### 

SAR <- SAR %>% 
  select(oceandist_km, Basin, MigrYr, corSAR_transport, origin.f, dams, Basin) %>% 
  mutate(Year = MigrYr) %>% 
  filter(origin.f != "Warm Springs") 

SAR <- na.omit(SAR) %>% 
  mutate(origin.f = droplevels(origin.f))

unique(SAR $origin.f)

PBTxSAR <- left_join(trawl_2015_2021, SAR, by = c("Year"= "Year","origin.f"= "origin.f")) %>% 
  mutate(year.f = as.factor(Year))

PBTxSAR$SAR_mod_prop <- PBTxSAR$corSAR_transport/100
PBTxSAR <- PBTxSAR[!is.na(PBTxSAR$SAR_mod_prop), ]
PBTxSAR$logitSAR <- logit(PBTxSAR$SAR_mod_prop)

PBTxSAR.1<- PBTxSAR %>% 
  select(Basin, FultonCondition, FieldLength_mm, fish_mass, origin.f, year.f, corSAR_transport, SAR_mod_prop, oceandist_km, logitSAR, dams, origin.f, PBTHatcheryAssignment)

PBTxSAR.1 <- PBTxSAR.1[!is.na(PBTxSAR.1$SAR_mod_prop), ]


PBTxSAR.1  <-PBTxSAR.1   %>%
  dplyr::group_by(origin.f, year.f ) %>%  # Group by 'origin' and 'Year'
  filter(n() >= 5) %>%
  ungroup() %>% 
  mutate(year.f = droplevels(year.f)) %>% 
  mutate(origin.f = droplevels(origin.f)) 


####----USE ALL data to model SARS----####

#Drop 2017 due to only Carson Hatchery represented
PBTxSAR.1 <- PBTxSAR.1 %>% 
  filter(year.f != "2017")

summaryxhatchery <- PBTxSAR.1 %>%
  dplyr::group_by(origin.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group

PBTxSAR.1.n <- PBTxSAR.1 %>% 
  select(year.f, origin.f, FieldLength_mm) %>%
  group_by(origin.f, year.f) %>% 
  dplyr::summarise_all(funs(n = n()))

#Model selection

lm <- lm(logitSAR ~ origin.f , data = PBTxSAR.1)
lm1 <- lmer(logitSAR ~ origin.f + (1|year.f), data = PBTxSAR.1)
lm2 <- lmer(logitSAR ~ origin.f + FultonCondition + (1|year.f), data = PBTxSAR.1, REML = FALSE)
lm3 <- lmer(logitSAR ~ origin.f * FultonCondition + (1|year.f), data = PBTxSAR.1, REML = FALSE)
lm4 <- lmer(logitSAR ~ origin.f * FieldLength_mm + (1|year.f), data = PBTxSAR.1, REML = FALSE)
lm5 <- lmer(logitSAR ~ origin.f + FieldLength_mm + (1|year.f), data = PBTxSAR.1, REML = FALSE)

AIC(lm, lm1, lm2, lm3, lm4, lm5) 

check_collinearity(lm3)
check_collinearity(lm4)
check_collinearity(lm2)

summary(lm2)

vif(lm3) 
vif(lm4) 
vif(lm2) 
 
#Final Model

lm2 <- lmer(logitSAR ~ origin.f + FultonCondition + (1|year.f), data = PBTxSAR.1, REML = TRUE)

r.squaredGLMM(lm2)
summary(lm2)
Anova(lm2)
check_model(lm2)

table <- tab_model(lm2, show.se = TRUE, show.p = TRUE, p.style = "stars")
print(table)


#####-----Compare ocean distance, dams and condition models-----######

lm2 <- lmer(logitSAR ~ origin.f + FultonCondition + (1|year.f), data = PBTxSAR.1, REML = FALSE)
lm2.1 <- lmer(logitSAR ~ oceandist_km +FultonCondition + (1|year.f), data = PBTxSAR.1, REML =FALSE)
lm3.1 <- lmer(logitSAR ~ dams +FultonCondition + (1|year.f), data = PBTxSAR.1, REML =FALSE)
lm4.11 <- lmer(logitSAR ~ oceandist_km + dams + FultonCondition +(1|year.f), data = PBTxSAR.1, REML =FALSE)

AIC(lm2, lm2.1, lm3.1, lm4.11)

lm2 <- lmer(logitSAR ~ origin.f + FultonCondition + (1|year.f), data = PBTxSAR.1, REML = TRUE)
lm2.1 <- lmer(logitSAR ~ oceandist_km +FultonCondition + (1|year.f), data = PBTxSAR.1, REML =TRUE)
lm3.1 <- lmer(logitSAR ~ dams +FultonCondition + (1|year.f), data = PBTxSAR.1, REML =TRUE)
lm4.11 <- lmer(logitSAR ~ oceandist_km + dams + FultonCondition +(1|year.f), data = PBTxSAR.1, REML =TRUE)

r.squaredGLMM(lm2) #baseline
r.squaredGLMM(lm2.1) #oceandist_km +FultonCondition 
r.squaredGLMM(lm3.1) #dams +FultonCondition 
r.squaredGLMM(lm4.11) # oceandist_km + dams 
vif(lm3.1)



####----Extract EM Mariginal Means from Final SARs Model----####

emmeans <- emmeans(lm2 ,  ~  origin.f )
emmeans
summary(emmeans)

emmeansall <- emmeans(lm2 , ~  FultonCondition * origin.f )
summary(emmeansall)


####----Extract marginal means for full model using ggpredict----####

mydf <- predict_response(lm2, terms = c( "FultonCondition", "origin.f"))

mydf$FultonCondition  <- mydf$x
mydf$origin.f <- mydf$group
mydf$SAR_mod<- invlogit(mydf$predicted)*100 #back-transform predictions

sarplot <- ggplot(NULL, mapping = aes()) + 
  geom_point(data =PBTxSAR.1, aes(x = FultonCondition , y = corSAR_transport, color = origin.f), size = 3, alpha = .8)+
  
  geom_line(data = mydf, aes(x = FultonCondition, y = SAR_mod), size =1.5, color = "darkgrey")+
  labs(x = "Condition at capture, mm", y = "Predicted SAR")+
  scale_color_manual(values = c("Carson" = "#481D6FFF", 
                                "Leavenworth" = "#277F8EFF", 
                                "Winthrop" ="#277F8EFF",
                                "Clearwater" = "#71CF57FF", 
                                "Dworshak" = "#71CF57FF", 
                                "McCall" = "#71CF57FF", 
                                "Pahsimeroi" = "#71CF57FF", 
                                "Rapid River" = "#71CF57FF", 
                                "Sawtooth" = "#71CF57FF")) +
  ggtitle("") +
  #ylim(0,1.75)+
  #facet_wrap(~origin.f) +
  facet_wrap(~factor(origin.f, c("Carson", "Leavenworth", "Winthrop", "Clearwater", "Dworshak",  "McCall", "Pahsimeroi",  "Rapid River", "Sawtooth")))+
  theme_classic(base_size =18, base_family = "")+
  theme(legend.position = "none")

sarplot

ggsave("FinalSARWoJACKmodplotINORDER_CORSAR.png", 
       dpi = 300, width = 12, height = 8, units = "in", bg = "white")


####-----Annual means by hatchery plot----#####

PBTxSAR.1$predicted <- predict(lm2, PBTxSAR.1, re.form = NULL)  # Includes random effects


# Define a function to generate predictions for bootstrapping
boot_predict <- function(model, data) {
  predict(model, newdata = data, re.form = NULL)
}

# Perform bootstrapping
set.seed(123)  # For reproducibility
boot_results <- boot(data = PBTxSAR.1, statistic = function(data, i) {
  boot_predict(lm2, data[i, ])
}, R = 1000)  # R is the number of bootstrap samples

# Calculate confidence intervals
ci <- apply(boot_results$t, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
PBTxSAR.1$conf.low <- ci[1, ]
PBTxSAR.1$conf.high <- ci[2, ]


######-------Plot predicted means per year and origin------######

hatch_colors <- c("Carson" = "#645394","Winthrop" = "lightblue", "Leavenworth" = "#2171B5",
                  "Clearwater" = "forestgreen", "Dworshak" = "yellowgreen","McCall" = "mediumseagreen","Pahsimeroi" ="palegreen" , "Rapid River" ="#808000" , "Sawtooth" ="#71CF57FF" )

PBTxSAR.1$origin.f <- factor(PBTxSAR.1$origin.f,
                             levels = c("Carson", "Sawtooth", "Winthrop", "Clearwater", "Leavenworth",
                                        "Pahsimeroi", "Dworshak",
                                        "McCall", "Rapid River"))
summary_df <- PBTxSAR.1 %>%
  group_by(year.f, origin.f) %>%
  summarise(
    mean_predicted = mean(predicted, na.rm = TRUE),
    sd_predicted = sd(predicted, na.rm = TRUE),
    n = n(),
    se_predicted = sd_predicted / sqrt(n),
    ci_lower = mean_predicted - qt(0.975, n - 1) * se_predicted,
    ci_upper = mean_predicted + qt(0.975, n - 1) * se_predicted
  )

SARxORIGINeachYear <- ggplot(summary_df, aes(x = year.f, y = mean_predicted, color = origin.f, group = origin.f)) +
  geom_errorbar(aes(ymin = mean_predicted - sd_predicted, ymax = mean_predicted + sd_predicted), color = "black", width = 0.2) + # For SD
  geom_point(size = 4) +
  scale_color_manual(values = hatch_colors) +
  labs(x = "", y = "Predicted SAR", title = "") +
  guides(col = guide_legend(title = "")) +
  theme_classic(base_size = 24) +
  theme(legend.position = "top")


SARxORIGINeachYear


ggsave("SARWOJACKCORSARxORIGINeachYear.topleg.png", 
       dpi = 300, width = 14, height = 8, units = "in", bg = "white")

#####