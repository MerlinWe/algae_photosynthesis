###############################################################################################################
#################################### Temperature and Light Stress on Algae #################################### 
###############################################################################################################

# Hesse et al. (2025)

rm(list = ls()) # clean global environment
options(scipen = 999, digits = 4) 

library(cowplot)     
library(MuMIn)      
library(lme4)        
library(lmerTest)   
library(emmeans)    
library(car)        
library(report)     
library(MASS)       
library(knitr)      
library(flextable)   
library(officer)    
library(broom)      
library(ggthemes)   
library(tidyverse)  

path <- "/Users/serpent/Documents/Projects/Giglio/Code/" # to WD 

algae <- read_csv(paste0(path, "giglio_algae.csv"))
glimpse(algae)

dat <- algae %>%
	
	# filtering
	dplyr::filter(Depth == 20) %>%
	dplyr::filter(Temperature != 34) %>%
	dplyr::select(-SampleiD, -Depth, -Date) %>%
	
	# categorize temperature
	dplyr::mutate(
		Temp_cat = case_when(
			Temperature == 21 ~ "control",
			Temperature ==  26 ~ "warm",
			Temperature == 30 ~ "hot",
			TRUE ~ "NA")) %>%
	
	# set formats 
	dplyr::mutate(
		Respiration = as.double(Respiration),
		netPhotosynthesis = as.double(netPhotosynthesis),
		grossPhotosynthesis = as.double(grossPhotosynthesis),
		PR_ratio = as.double(PR_ratio),
		
		Algae = as.factor(Algae),
		Species = as.factor(Species),
		Light = factor(Light, levels = c("control", "medium", "high")),
		Temp_cat = factor(Temp_cat, levels = c("control", "warm", "hot"))) %>%
	
	mutate(PR_ratio = grossPhotosynthesis / abs(Respiration)) # fix pr ratio


# Linearity? Can we treat temperature & light as continues variable? 

dat %>% ggplot(aes(x = Temperature, y = netPhotosynthesis, color = Species)) +
	geom_point(alpha = .5) +
	geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +  # Linear fit
	geom_smooth(method = "loess", linetype = "dashed", se = FALSE) +  # Non linear smoothed fit
	facet_wrap(~Species) +
	theme_bw()

dat %>% ggplot(aes(x = Temperature, y = PR_ratio, color = Species)) +
	geom_point(alpha = .5) +
	geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +  # Linear fit
	geom_smooth(method = "loess", linetype = "dashed", se = FALSE) +  # Non linear smoothed fit
	facet_wrap(~Species) +
	theme_bw()

dat %>% 
	mutate(Light = case_when(
			Light == "control" ~ "1",
			Light ==  "medium" ~ "2",
			Light == "high"    ~ "3",
			TRUE ~ "NA"),
			Light = as.numeric(Light)) %>%
	ggplot(aes(x = Light, y = netPhotosynthesis, color = Species)) +
	geom_point(alpha = .5) +
	geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +  # Linear fit
	geom_smooth(method = "loess", linetype = "dashed", se = FALSE) +  # Non linear smoothed fit
	facet_wrap(~Species) +
	theme_bw()

dat %>% 
	mutate(Light = case_when(
		Light == "control" ~ "1",
		Light ==  "medium" ~ "2",
		Light == "high"    ~ "3",
		TRUE ~ "NA"),
		Light = as.numeric(Light)) %>%
	ggplot(aes(x = Light, y = PR_ratio, color = Species)) +
	geom_point(alpha = .5) +
	geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +  # Linear fit
	geom_smooth(method = "loess", linetype = "dashed", se = FALSE) +  # Non linear smoothed fit
	facet_wrap(~Species) +
	theme_bw()

# not really... 

##### Fig 1: descriptive plot

temperature_labels <- c("control" = "21 °C", "warm" = "26 °C", "hot" = "30 °C")
p1 <- ggplot(dat, aes(x = Temp_cat, y = netPhotosynthesis, color = Species)) +
	geom_point(alpha = 0.4, position = position_jitter(width = 0.2, height = 0)) + 
	geom_smooth(aes(group = Species), method = "loess", se = FALSE) + 
	scale_x_discrete(labels = temperature_labels) +
	scale_colour_manual(values = c("sienna", "darkgreen", "red")) +
	labs(x = "Temperature (°C)", y = "Net Photosynthesis (µmol O2 cm² h⁻¹)", title = "Main effect of temperature") +
	theme_few() +
	theme(
		text = element_text(size = 12, family = "Sans"),
		strip.text = element_text(size = 12, face = "bold"),
		plot.title =  element_text(size = 12, face = "bold"),
		legend.position = "none")

p2 <- ggplot(dat, aes(x = Light, y = netPhotosynthesis, color = Species)) +
	geom_point(alpha = 0.4, position = position_jitter(width = 0.2, height = 0)) + 
	geom_smooth(aes(group = Species), method = "loess", se = FALSE) +   
	scale_colour_manual(values = c("sienna", "darkgreen", "red")) +
	labs(x = "Light intensity", y = NULL, title = "Main effect of light") +
	theme_few() +
	theme(
		text = element_text(size = 12, family = "Sans"),
		strip.text = element_text(size = 12, face = "bold"),
		plot.title =  element_text(size = 12, face = "bold"),
		legend.position = "none")

p3 <- ggplot(dat, aes(x = Temperature, y = netPhotosynthesis, color = Species, group = Species)) +
	geom_point(alpha = 0.4, position = position_jitter(width = 0.2, height = 0)) + 
	geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +  
	facet_wrap(~ Light) +  
	scale_colour_manual(values = c("sienna", "darkgreen", "red")) +	
	labs(
		x = "Temperature (°C)",
		y = "Net Photosynthesis (µmol O2 cm² h⁻¹)",
		title = "Temperature x Light Interaction",
		color = NULL) +
	theme_few() + 
	theme(text = element_text(size = 12, family = "Sans"), 
				plot.title =  element_text(size = 12, face = "bold"),
		strip.text = element_text(size = 12), 
		legend.position = "bottom")

plot_desc <- plot_grid(plot_grid(p1, p2, rel_widths = c(.54, .5)), 
											 p3, ncol = 1, nrow = 2,
											 rel_heights = c(.4, .6))
ggsave(filename = "/Users/serpent/Desktop/desc_plot.png",
			 plot = plot_desc, 
			 bg = "transparent",
			 width = 290, 
			 height = 220, 
			 units = "mm", 
			 dpi = 800)

##### Initial model (full factorial LMM) for netPhotosynthesis

# with random intercept
mod1 <- lmer(netPhotosynthesis ~ Temp_cat * Light * Species + (1 | Species), data = dat)

# with random slope 
mod2 <- lmer(netPhotosynthesis ~ Temp_cat * Light * Species + (Temp_cat + Light | Species), data = dat)

anova(mod1, mod2) # no sig improvement for random slopes but overcomplicated... 

# continue with mod1
summary(mod1)

# species variation is very small... prob no random effect needed 

mod3 <- lm(netPhotosynthesis ~ Temp_cat * Light * Species, data = dat)

anova(mod3)
summary(mod3)
par(mfrow=c(2,2))
plot(mod3) # not too bad 
hist(mod3$residuals)

# do a backwards stepwise selection

mod_final <- stepAIC(mod3, direction = "backward")
summary(mod_final)

summary(mod_final)
AIC(mod_final)

plot(fitted(mod_final), residuals(mod_final),
		 xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red", lty = 2)

qqnorm(residuals(mod_final))
qqline(residuals(mod_final), col = "red")

vif(mod_final, type = "predictor") # ok 

summary(mod_final)
rm(mod1, mod2, mod3)


# Get a clean table for the paper
model_results <- tidy(mod_final)
model_results <- model_results[, c("term", "estimate", "std.error", "statistic", "p.value")]
colnames(model_results) <- c("Parameter", "Estimate", "SE", "T Ratio", "p-value")
model_table <- flextable(model_results) %>%
	set_caption("Table 1. Parameter estimates from the final linear mixed-effects model.") %>%
	colformat_double(j = c("Estimate", "SE", "T Ratio", "p-value"), digits = 3) %>%
	autofit()
doc <- read_docx() %>%
	body_add_flextable(model_table) %>%
	body_add_par("")
# print(doc, target = "model_estimates.docx") 

# Extract estimated marginal means
emm_results <- emmeans(mod_final, ~ Species | Temp_cat * Light)

# Plot estimated marginal means
emm_df <- as.data.frame(emm_results)
temperature_labels <- c("21 °C", "26 °C", "30 °C")
names(temperature_labels) <- c("control", "warm", "hot")

emm_plot <- ggplot(emm_df, aes(x = Species, y = emmean, fill = Species)) +
	
	geom_bar(stat = "identity", position = "dodge", color = "black", alpha = .5, width = .8) +
	geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) + 
	
	scale_fill_manual(values = c("sienna", "darkgreen", "red")) +
	
	facet_grid(Light ~ Temp_cat, scales = "free_y", labeller = labeller(Temp_cat = temperature_labels)) +
	
	labs(
		x = NULL,
		y = "Net Photosynthesis (µmol O2 cm² h⁻¹)") +
		
	theme_few() +
		
	theme(
		text = element_text(size = 12, family = "Sans"),
		strip.text = element_text(size = 12, face = "bold"),
		legend.position = "none")

ggsave(filename = "/Users/serpent/Desktop/emm_factorial_plot.png",
			 plot = emm_plot, 
			 bg = "transparent",
			 width = 290, 
			 height = 190, 
			 units = "mm", 
			 dpi = 800)

# Get a table of pairwise comparisons
contrast <- pairs(emm_results, adjust = "tukey")
contrast <- as.data.frame(summary(contrast))

# Create a formatted flextable
ft <- flextable(contrast) %>%
	set_caption("Pairwise Comparisons of Species") %>%
	colformat_double(j = c("estimate", "SE", "t.ratio", "p.value"), digits = 3) %>%
	autofit()

# Export to Word
doc <- read_docx() %>%   # Create a new Word document
	body_add_flextable(ft) %>%  # Add the flextable
	body_add_par("")  # Add a space after the table
# print(doc, target = "contrast_results.docx")  # Save as Word file



#######################################################################################################
######################################## Analyzing at PR ratio ######################################## 
#######################################################################################################


# with random intercept
mod1 <- lmer(PR_ratio ~ Temp_cat * Light * Species + (1 | Species), data = dat)

# with random slope 
mod2 <- lmer(PR_ratio ~ Temp_cat * Light * Species + (Temp_cat + Light | Species), data = dat)

anova(mod1, mod2) # random slopes do no good 

# continue with mod1
summary(mod1)
AIC(mod1)

anova(mod1)

plot(mod1) # problematic
hist(residuals(mod1, method = "pearson"))

qqnorm(residuals(mod1))
qqline(residuals(mod1), col = "red") # not great


# species variation is again very small... prob no random effect needed 
mod3 <- lm(PR_ratio ~ Temp_cat * Light * Species, data = dat)

anova(mod3)
summary(mod3)
par(mfrow=c(2,2))
plot(mod3) 
dev.off()
hist(mod3$residuals)

# maybe try gamma 

mod_4 <- glm(PR_ratio ~ Temp_cat * Light * Species, 
							 family = Gamma(link = "log"), data = dat)

anova(mod_4)
summary(mod_4)
par(mfrow=c(2,2))
plot(mod_4) # better

# still some normality problems, but we're not looking at p-values so it's fine. 
# Estimates and emmeans should be reliable enough... 

dev.off()
hist(mod_4$residuals) # better than lm 

mod_final <- stepAIC(mod_4, direction = "backward")
summary(mod_final) # same as before 

# this is the top model 
mod_pr <- glm(PR_ratio ~ Temp_cat * Light * Species, 
						 family = Gamma(link = "log"), data = dat)


# estimate marginal means 
emmeans(mod_pr, ~ Species | Temp_cat * Light, type = "response")
emm_pr <- as.data.frame(emmeans(mod_pr, ~ Species | Temp_cat * Light, type = "response"))

# colours and breaks
pr_colors <- c("darkred", "orange", "darkgreen", "blue")
pr_breaks <- c(0.5, 1.0, 3.0, 5.0, max(emm_pr$response))

# heatmap
pr_plot <- emm_pr %>%
	mutate(label = paste0(round(response, 2), "\n± ", round(SE, 2))) %>%  # Add line break before SE
	ggplot(aes(x = Temp_cat, y = Species, fill = response)) +
	geom_tile(color = "white", size = 0.3) + 
	facet_wrap(~ Light, labeller = labeller(Light = c("control" = "Control Light", 
																										"medium" = "Medium Light", 
																										"high" = "High Light"))) + 
	scale_fill_gradientn(
		colors = pr_colors, 
		values = scales::rescale(pr_breaks), # scale breaks to fit gradient
		limits = c(0, max(emm_pr$response)), 
		name = "P:R Ratio") +
	scale_x_discrete(labels = temperature_labels) +
	geom_text(aes(label = label), color = "white", fontface = "bold", size = 4.5) + 
	labs(
		x = NULL,
		y = NULL) +
	theme_few(base_size = 14, base_family = "Sans") + 
	theme(
		text = element_text(size = 12),
		strip.text = element_text(size = 12, face = "bold"), 
		legend.position = "right",
		axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
		axis.text.y = element_text(size = 12),
		plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(filename = "/Users/serpent/Desktop/pr_plot.png",
			 plot = pr_plot, 
			 bg = "transparent",
			 width = 290, 
			 height = 190, 
			 units = "mm", 
			 dpi = 800)
