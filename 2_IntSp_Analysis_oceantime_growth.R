# Code to examine differences in size at ocean entry, timing of ocean entry, and ocean
# growth rate for Interior Spring Chinook salmon 
#Updated on 10/16/2025


# Vector of required packages
pkgs <- c(
  "rstudioapi", "tidyverse", "ggplot2", "dplyr", "ggpubr", "MuMIn", 
  "emmeans", "car", "nlme", "lme4", "ggeffects", "hrbrthemes",
  "LaplacesDemon", "performance", "plyr", "broom.mixed", "report", "ggridges", "viridis"
)

# Install any missing packages
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(pkgs[!installed])
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#####-------------Import data for 2021 oe timing-----#######

oe <- read.csv("2_OceanMetrics.csv") %>% 
  mutate(origin = as.factor(PBTHatcheryAssignment)) %>% 
  mutate(year.f = as.factor(Year)) %>% 
  mutate(origin.f = as.factor(origin))

oe <-oe %>%
  group_by(PBTHatcheryAssignment, Year) %>%
  filter(n() >= 5) %>%
  ungroup()

summary <- oe%>%
  dplyr::group_by(origin, year.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group

summary
sum(summary$sample_size)

smolt_windows <- data.frame(
  year.f = c("2016", "2016", 
             "2018", 
             "2021", "2021", "2021", "2021"),
  Basin = c("Wild", "Snake", 
            "Snake", 
            "Wild", "Mid", "Upper", "Snake"),
  xmin = c(105, 116,
           117,
           123, 107, 125, 128),
  xmax = c(153, 146,
           156,
           166, 145, 156, 156),
  ymin = -Inf,
  ymax = Inf
)


####----Test effect of Origin on timing of entry----#####  

timing <-oe %>%
  filter(!is.na(oe.julian.day))%>%
  dplyr::select(origin.f, year.f, oe.julian.day)%>%
  mutate(oe.julian.day = as.numeric(oe.julian.day)) %>%
  filter(year.f != "2019") %>%
  mutate(year.f = droplevels(year.f)) %>%
  group_by(origin.f, year.f) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  mutate(origin.f = droplevels(origin.f))

summarytiming <- timing%>%
  dplyr::group_by(origin.f, year.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group

summarytimingxorigin <- timing%>%
  dplyr::group_by(origin.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group

oeo_model <- lmer(oe.julian.day ~ origin.f + (1|year.f), data = timing) 
summary(oeo_model)
Anova(oeo_model)

r.squaredGLMM(oeo_model )
check_model(oeo_model)

tidy_output <- tidy(oeo_model)
print(tidy_output)
write.csv(tidy_output, "lmm_results.csv", row.names = FALSE)

report_parameters(oeo_model)
report_statistics(oeo_model)

emms_oeo <- emmeans(oeo_model, "origin.f")
pairwise_oeo  <- pairs(emms_oeo , adjust = 'bonferroni')
summary(pairwise_oeo )


####----Entry Timing Plot----####

timing2 <-oe %>%
  filter(!is.na(oe.julian.day))%>% 
  dplyr::select(origin.f, year.f, oe.julian.day, Basin)%>%
  filter(year.f != "2015") %>% 
  filter(year.f != "2017") %>%
  filter(year.f != "2019") %>%
  mutate(year.f = droplevels(year.f)) %>% 
  group_by(origin.f, year.f) %>%
  filter(n() >= 5) %>%
  ungroup()

timing2$origin.f <- factor(timing2$origin.f, levels = rev(c(
  "Carson", "Little White Salmon", "Walla Walla", "Chief Joseph", "Sawtooth", "Dworshak", "Lookingglass", 
  "Rapid River", "Pahsimeroi", "McCall", "Clearwater", "Wild"
)))

timing2$Basin <- factor(timing2$Basin, levels = c("Mid", "Upper", "Snake", "Wild"))

#PLOT

y_positions <- data.frame(
  origin.f = c("Mid", "Upper", "Snake", "Wild"),
  y = c(4, 3, 2, 1)  
)

# Join to smolt_windows and calculate ymin/ymax offsets
smolt_windows <- smolt_windows %>%
  left_join(y_positions, by = c("Basin" = "origin.f")) %>%
  mutate(
    ymin = y - 0.45,
    ymax = y - 0.25
  )

ggplot(timing2, aes(x = oe.julian.day, y = origin.f, fill = Basin)) +
  
  geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = 0.6) +
  
  # smolt over bon shading
  geom_rect(data = smolt_windows,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Basin),
            inherit.aes = FALSE,
            alpha = 0.8, color = NA) +
  # Cruise window for 2016
  geom_rect(data = data.frame(year.f = "2016", xmin = 145, xmax = 151,
                              ymin = -Inf, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "white", alpha = 0.25) +
  geom_vline(data = data.frame(year.f = "2016", day = c(145, 151)),
             aes(xintercept = day),
             inherit.aes = FALSE, linetype = "dashed", color = "grey30", size = 0.9) +
  
  # Cruise window for 2018
  geom_rect(data = data.frame(year.f = "2018", xmin = 144, xmax = 148,
                              ymin = -Inf, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "white", alpha = 0.25) +
  geom_vline(data = data.frame(year.f = "2018", day = c(144, 148)),
             aes(xintercept = day),
             inherit.aes = FALSE, linetype = "dashed", color = "grey30", size = 0.9) +
  
  # Cruise window for 2021
  geom_rect(data = data.frame(year.f = "2021", xmin = 142, xmax = 146,
                              ymin = -Inf, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "white", alpha = 0.25) +
  geom_vline(data = data.frame(year.f = "2021", day = c(142, 146)),
             aes(xintercept = day),
             inherit.aes = FALSE, linetype = "dashed", color = "grey30", size = 0.9) +
  
  scale_fill_manual(values = c(
    "Mid" = "#481D6FFF",
    "Upper" = "#277F8EFF",
    "Snake" = "#71CF57FF",
    "Wild" = "#FDE725FF"
  )) +
  xlab("Marine entry day of year") +
  ylab("") +
  theme_ipsum() +
  facet_wrap(~ year.f) +
  theme(
    legend.position = "right",
    panel.spacing = unit(1.2, "lines"),
    strip.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24, face = "bold", hjust = 0.5) ,
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold")
    
  )


ggsave("PIT_TimingofEntry.png", 
       dpi = 300, width = 14, height = 8, units = "in", bg = "white")


####----Test effect of Origin on entry size----#####  

oesize <- oe %>%
  filter(!is.na(oe_ow))%>% 
  dplyr::select(origin.f, year.f, oe_ow)%>%
  filter(year.f != "2019") %>%  
  mutate(year.f = droplevels(year.f)) %>% 
  group_by(origin.f, year.f) %>%
  filter(n() >= 5) %>%
  ungroup()%>% 
  mutate(origin.f = droplevels(origin.f)) 

summaryoesizexorigin <- oesize%>%
  dplyr::group_by(origin.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group

oesize_model <- lmer(oe_ow ~ origin.f + (1|year.f), data = oesize) # just do this as question in origin
summary(oesize_model)
Anova(oesize_model)
r.squaredGLMM(oesize_model)

tidy_output <- tidy(oesize_model)
print(tidy_output)
write.csv(tidy_output, "lmm_results.csv", row.names = FALSE)

report_parameters(oesize_model)
report_statistics(oesize_model)

emms_oesize <- emmeans(oesize_model, "origin.f")
pairwise_oesize  <- pairs(emms_oeo , adjust = 'bonferroni')
summary(pairwise_oesize)


####----Size at Entry Plot----####


oesize2 <- oe %>%
  filter(!is.na(oe.julian.day))%>% 
  dplyr::select(origin.f, year.f, oe_ow, Basin)%>%
  filter(year.f != "2015") %>% 
  filter(year.f != "2017") %>%
  filter(year.f != "2019") %>%  
  mutate(year.f = droplevels(year.f)) %>% 
  group_by(origin.f, year.f) %>%
  filter(n() >= 5) %>%
  ungroup()

oesize2$year <- as.numeric(oesize2$year.f)

oesize2$origin.f <- factor(oesize2$origin.f, levels = rev(c(
  "Carson", "Little White Salmon", "Walla Walla", "Chief Joseph", "Sawtooth", "Dworshak", "Lookingglass", 
  "Rapid River", "Pahsimeroi", "McCall", "Clearwater", "Wild"
)))

oesize2$Basin <- factor(oesize2$Basin, levels = c("Mid", "Upper", "Snake", "Wild"))

# Plot

ggplot(oesize2, aes(x = oe_ow, y = origin.f, fill = Basin)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = 0.6) +
  scale_fill_manual(values = c(
    "Mid" = "#481D6FFF",
    "Upper" = "#277F8EFF",
    "Snake" = "#71CF57FF",
    "Wild" = "#FDE725FF"
  )) +
  theme_ipsum() +
  facet_wrap(~ year.f)+
  xlab("Otolith size at marine entry (µm)")+
  ylab("") +
  theme(
    legend.position = "right",
    panel.spacing = unit(1.2, "lines"),
    strip.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24, face = "bold", hjust = 0.5) ,
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"))

ggsave("SizeatEntry.png", 
       dpi = 300, width = 14, height = 8, units = "in", bg = "white")



####----Test effect of Origin on marine growth----#####  


oegrow <- oe %>%
  filter(!is.na(ocean.growth))%>% 
  dplyr::select(origin.f, year.f, ocean.growth)%>%
  filter(year.f != "2019") %>%  
  mutate(year.f = droplevels(year.f)) %>% 
  group_by(origin.f, year.f) %>%
  filter(n() >= 5) %>%
  ungroup()%>% 
  mutate(origin.f = droplevels(origin.f)) 

oegrow$ocean.growth <- as.numeric(oegrow$ocean.growth)

summaryoegrowxorigin <- oegrow%>%
  dplyr::group_by(origin.f) %>%  # Group by 'origin' and 'Year'
  dplyr::summarise(sample_size = n(), .groups = "drop")  # Count rows for each group


oegrow_model <- lmer(ocean.growth~ origin.f + (1|year.f), data = oegrow) # just do this as question in origin
summary(oegrow_model)
Anova(oegrow_model)
r.squaredGLMM(oegrow_model)

check_model(oegrow_model)

tidy_output <- tidy(oegrow_model)
print(tidy_output)
write.csv(tidy_output, "lmm_results.csv", row.names = FALSE)

report_parameters(oegrow_model)
report_statistics(oegrow_model)

emms_oegrow  <- emmeans(oegrow_model , "origin.f")
pairwise_oegrow   <- pairs(emms_oegrow, adjust = 'bonferroni')
summary(pairwise_oegrow )


####----Marine Growth Plot----####


oegrow2 <- oe %>%
  filter(!is.na(ocean.growth))%>% 
  dplyr::select(origin.f, year.f, ocean.growth, Basin)%>%
  filter(year.f != "2015") %>% 
  filter(year.f != "2017") %>%
  filter(year.f != "2019") %>%  
  mutate(year.f = droplevels(year.f)) %>% 
  group_by(origin.f, year.f) %>%
  filter(n() >= 5) %>%
  ungroup()

oegrow2$ocean.growth <- as.numeric(oegrow2$ocean.growth)

oegrow2$origin.f <- factor(oegrow2$origin.f, levels = rev(c(
  "Carson", "Little White Salmon", "Walla Walla", "Chief Joseph", "Sawtooth", "Dworshak", "Lookingglass", 
  "Rapid River", "Pahsimeroi", "McCall", "Clearwater", "Wild"
)))

oegrow2$Basin <- factor(oegrow2$Basin, levels = c("Mid", "Upper", "Snake", "Wild"))

# Plot
ggplot(oegrow2, aes(x = ocean.growth, y = origin.f, fill = Basin)) +
  geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = 0.6) +
  scale_fill_manual(values = c(
    "Mid" = "#481D6FFF",
    "Upper" = "#277F8EFF",
    "Snake" = "#71CF57FF",
    "Wild" = "#FDE725FF"
  )) +
  xlab("Otolith marine growth (µm)")+
  ylab("") +
  theme_ipsum() +
  facet_wrap(~ year.f)+
  theme(
    legend.position = "right",
    panel.spacing = unit(1.2, "lines"),
    strip.text.x = element_text(size = 24),
    axis.title.x = element_text(size = 24, face = "bold", hjust = 0.5) ,
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, face = "bold"))

ggsave("Marinegrow.png", 
       dpi = 300, width = 14, height = 8, units = "in", bg = "white")

######
