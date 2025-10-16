#Code to test for differences in length, mass, and condition at the hatchery and at capture for JSOES
#Interior Spring Chinook yearlings from 2015-2019 and 2021
#Updated on 10/16/2025

# Vector of required packages
pkgs <- c(
  "rstudioapi", "tidyverse", "ggplot2", "dplyr", "ggpubr", "car", "nlme", "MuMIn", "lme4", 
  "emmeans", "ggeffects", "performance"
)

# Install any missing packages
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(pkgs[!installed])
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #


#####-------------Import data for JSOES - HATCHERY SIZE AND MASS 2015-2021 capture comparison------#######


compare <- read.csv("3_SIZEMASSCOMP.csv") %>%
  mutate(year = as.factor(year),origin = as.factor(origin)) %>%
  filter(year != "2017")

compare <- compare %>% 
  filter(!(origin == "Chief Joseph" & year == "2021")) %>% 
  mutate(year = droplevels(year)) %>% 
  mutate(origin = droplevels(origin)) 

comparelength <- lmer(cap.length ~ h.length + (1|year), data = compare, REML= TRUE)

summary(comparelength)
Anova(comparelength)
r.squaredGLMM(comparelength)

check_model(comparelength)

comparemass <- lmer(cap.mass~ h.mass  + (1|year), data = compare, REML= TRUE)
comparemass
r.squaredGLMM(comparemass)
Anova(comparemass)

ggplot(compare, aes(x = h.length, y = cap.length, color = year)) +
  geom_point(size = 6) +
  geom_smooth(method = lm, se= FALSE, size =2)+
  theme_classic()+
  theme_classic(base_size = 11) +
  theme(legend.position = "right", 
        text = element_text(size = 24),
        legend.text = element_text(size = 30),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 30, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold")) +
  labs(x = "Hatchery fork length, mm", y = "Fork length at capture, mm")

ggsave("HatchvsCaplength.png", 
       dpi = 300, width = 15, height = 10, units = "in", bg = "white")

oe <- compare %>% 
  filter(!is.na(oe.ow))
  
cor.test(oe$h.length, oe$oe.ow, data = oe)  

######