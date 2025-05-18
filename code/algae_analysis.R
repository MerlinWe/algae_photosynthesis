###############################################################################################################
######################### Temperature and Light Stress on Mediterranean Macro Algae ###########################
###############################################################################################################

rm(list = ls()) # make sure environment is clean 
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
library(multcomp)
library(ggthemes)   
library(tidyverse)  

setwd("/Users/merlin/Documents/Projects/Giglio") # set WD 
export = TRUE # export?  

algae <- read_csv(file.path("code", "giglio_algae.csv"))
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
	
	# Also; transform photosynthesis and respiration units from µmol O2 cm² h⁻¹ to nmol O2 m² s⁻¹ 
	
	# Convert from µmol O₂ cm⁻² h⁻¹ to nmol O₂ m⁻² s⁻¹:
	# Multiply by 1,000 to convert µmol to nmol,
	# Multiply by 10,000 to convert cm² to m²,
	# Divide by 3,600 to convert hours to seconds.
	# Total conversion factor = 1,000 × 10,000 ÷ 3,600 ≈ 2,777.78
	
	dplyr::mutate(
		netPhotosynthesis = netPhotosynthesis * 1e3 * 1e4 / 3600,
		grossPhotosynthesis = grossPhotosynthesis * 1e3 * 1e4 / 3600,
		Respiration = Respiration * 1e3 * 1e4 / 3600) %>%
	
	# fix PR ratio (there was a mistake in the data)
	mutate(PR_ratio = grossPhotosynthesis / abs(Respiration))



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

# not really... keep them categorical for analysis 

## =========================== Analyzing net photosynthesis ===========================

# Start with full-factorial mixed effects model

# with random intercept
mod1 <- lmer(netPhotosynthesis ~ Temp_cat * Light * Species + (1 | Species), data = dat) # unstable

# with random slope 
mod2 <- lmer(netPhotosynthesis ~ Temp_cat * Light * Species + (Temp_cat + Light | Species), data = dat) # unstable

anova(mod1, mod2) # no significant improvement for random slopes but overcomplicated... 

# continue with mod1
summary(mod1)

# species variation is very small... probably no random effect needed 

mod3 <- lm(netPhotosynthesis ~ Temp_cat * Light * Species, data = dat)

anova(mod3)
summary(mod3)

par(mfrow=c(2,2))
plot(mod3) # not too bad 

par(mfrow=c(1,1))
hist(mod3$residuals) # good 

# do a backwards stepwise selection

mod_final <- stepAIC(mod3, direction = "backward")
summary(mod_final)

summary(mod_final)
AIC(mod_final)


par(mfrow=c(1,2))
plot(fitted(mod_final), residuals(mod_final), xlab = "Fitted Values", ylab = "Residuals")
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

if (export) { print(doc, target = file.path("manuscript", "tables", "np_model_summary.docx")) }

# Extract estimated marginal means
emm_results <- emmeans(mod_final, ~ Species | Temp_cat * Light)

# Plot estimated marginal means with a raw data overlay
emm_df <- as.data.frame(emm_results)
temperature_labels <- c("21 °C", "26 °C", "30 °C")
names(temperature_labels) <- c("control", "warm", "hot")

# Extract letters for compact letter display
cld_df <- cld(emm_results, Letters = letters, adjust = "sidak") %>%
	as.data.frame() %>%
	dplyr::rename(label = .group) %>%
	mutate(y.pos = emmean + upper.CL * 0.25) %>%
	mutate(label = gsub(" ", "", label))

# Build Plot (Fig. 2)
emm_plot <- ggplot(emm_df, aes(x = Species, y = emmean, fill = Species)) +
	
	# EMM bars + CI 
	geom_bar(stat = "identity", position = "dodge", color = "black", alpha = .5, width = .8) +
	geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) + 
	
	# Add raw data as a jitter overlay to add transparency 
	geom_jitter(data = dat, aes(x = Species, y = netPhotosynthesis, color = Species),
							position = position_jitter(width = 0.15), alpha = 0.4, size = 1.2, inherit.aes = FALSE) +
	
	scale_fill_manual(values = c("sienna", "olivedrab3", "red1")) +
	scale_colour_manual(values = c("sienna", "olivedrab3", "red1")) +
	
	geom_text(data = cld_df, aes(x = Species, y = emmean * 0.05, label = label), 
		vjust = 0, size = 4, fontface = "bold") +

	facet_grid(Light ~ Temp_cat, scales = "free_y", labeller = labeller(Temp_cat = temperature_labels)) +
	
	labs(
		x = NULL,
		y = "Net Photosynthesis (nmol O₂ m⁻² s⁻¹)") +
	
	theme_few() +
	
	theme(
		text = element_text(size = 12, family = "sans"),
		strip.text = element_text(size = 12, face = "bold"),
		legend.position = "none")

if (export) {
ggsave(filename = file.path("manuscript", "figures", "fig_2.png"),
			 plot = emm_plot, 
			 bg = "transparent",
			 width = 290, 
			 height = 200, 
			 units = "mm", 
			 dpi = 600)
}

# Caption: Figure 2. Estimated marginal means (EMMs) of net photosynthesis across all combinations of temperature (21 °C, 26 °C, 30 °C) 
# and light (control, medium, high) for three Mediterranean macroalgal taxa (Cystoseira, Flabellia, Phyllophora). Bars show model-estimated 
# means from the full-factorial linear model with 95% confidence intervals (black whiskers), overlaid with raw data points. Letters within 
# each panel denote results of post-hoc pairwise comparisons (sidak-adjusted) between taxa: taxa that do not share a letter differ 
# significantly at α = 0.05. Note that identical predicted means across model and raw data are expected due to the balanced full-factorial design, 
# but confidence intervals reflect inferential uncertainty based on model structure and residual variance.

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
if (export) { print(doc, target = file.path("manuscript", "tables", "pairwise_comparisons.docx")) }

rm(cld_df, contrast, curWarnings, doc, emm_df, emm_results, ft, model_results, model_table) # clean environment

## =========================== Analyzing P:R ratio ===========================

# Again, start with full-factorial mixed effects model

# with random intercept
mod1 <- lmer(PR_ratio ~ Temp_cat * Light * Species + (1 | Species), data = dat) # unstable

# with random slope 
mod2 <- lmer(PR_ratio ~ Temp_cat * Light * Species + (Temp_cat + Light | Species), data = dat) # unstable

anova(mod1, mod2) # random slopes do no good 

# continue with mod1
summary(mod1)
AIC(mod1)
anova(mod1)

plot(mod1) # problematic , heteroscedasticity

par(mfrow=c(1,2))
hist(residuals(mod1, method = "pearson"))
qqnorm(residuals(mod1))
qqline(residuals(mod1), col = "red") # not great

# species variation is again very small... prob no random effect needed 
mod3 <- lm(PR_ratio ~ Temp_cat * Light * Species, data = dat)

anova(mod3)
summary(mod3)

par(mfrow=c(2,2))
plot(mod3) 
par(mfrow=c(1,1))
hist(mod3$residuals)

# still not great, maybe try gamma 

mod4 <- glm(PR_ratio ~ Temp_cat * Light * Species, 
						 family = Gamma(link = "log"), data = dat)

anova(mod4)
summary(mod4)

par(mfrow=c(2,2))
plot(mod4) # better

# still some deviance from gamma distrubtions problems, especially at high values.
# But we're not looking at p-values so it's fine. Estimates and emmeans should be reliable enough... 

par(mfrow=c(1,1))
hist(mod4$residuals) # better than the linear one 

mod_final <- stepAIC(mod4, direction = "backward")
summary(mod_final) # same as before , full-factorial model is best 

# this is the top model 
mod_pr <- glm(PR_ratio ~ Temp_cat * Light * Species, 
							family = Gamma(link = "log"), data = dat)

rm(mod1, mod2, mod3, mod4, mod_final) # clean environment

# estimate marginal means 
emmeans(mod_pr, ~ Species | Temp_cat * Light, type = "response")
emm_pr <- as.data.frame(emmeans(mod_pr, ~ Species | Temp_cat * Light, type = "response"))

# colours and breaks
pr_colors <- c("darkred", "orange", "darkgreen", "blue")
pr_breaks <- c(0.5, 1.0, 3.0, 5.0, max(emm_pr$response))

# Heatmap
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
	theme_few(base_size = 14, base_family = "sans") + 
	theme(
		text = element_text(size = 12),
		strip.text = element_text(size = 12, face = "bold"), 
		legend.position = "right",
		axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
		axis.text.y = element_text(size = 12),
		plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

if (export) {
	ggsave(filename = file.path("manuscript", "figures", "fig_3.png"),
				 plot = pr_plot, 
				 bg = "transparent",
				 width = 290, 
				 height = 190, 
				 units = "mm", 
				 dpi = 600)
}

## =========================== Supplementary ===========================

summary_data <- dat %>%
	group_by(Algae, Temperature, Light, Species) %>%
	summarize(
		netPhotosynthesis_mean = mean(netPhotosynthesis),
		grossPhotosynthesis_mean = mean(grossPhotosynthesis),
		Respiration_mean = mean(Respiration),
		PR_ratio_mean = mean(PR_ratio),
		PR_ratio_sd = sd(PR_ratio),
		netPhotosynthesis_se = sd(netPhotosynthesis) / sqrt(n()),
		grossPhotosynthesis_se = sd(grossPhotosynthesis) / sqrt(n()),
		Respiration_se = sd(Respiration) / sqrt(n())
	) %>%
	ungroup()

plot_data <- summary_data %>%
	pivot_longer(cols = c(netPhotosynthesis_mean, grossPhotosynthesis_mean, Respiration_mean),
							 names_to = "Variable", values_to = "Value") %>%
	mutate(Variable = factor(Variable, levels = c("netPhotosynthesis_mean", "grossPhotosynthesis_mean", "Respiration_mean")))

temperature_labels <- c("21 °C", "26 °C", "30 °C")
names(temperature_labels) <- c("21", "26", "30")

summary_data <- summary_data %>%
	group_by(Temperature, Light, Species) %>%
	mutate(max_value = max(c(netPhotosynthesis_mean, grossPhotosynthesis_mean, abs(Respiration_mean))) + 0.001)

summary_plot <- ggplot(plot_data, aes(x = Species, y = Value, fill = Variable)) +
	
	geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8, color = "black", alpha = .8) +
	
	geom_errorbar(aes(
		ymin = Value - case_when(
			Variable == "netPhotosynthesis_mean" ~ netPhotosynthesis_se,
			Variable == "grossPhotosynthesis_mean" ~ grossPhotosynthesis_se,
			Variable == "Respiration_mean" ~ Respiration_se),
		ymax = Value + case_when(
			Variable == "netPhotosynthesis_mean" ~ netPhotosynthesis_se,
			Variable == "grossPhotosynthesis_mean" ~ grossPhotosynthesis_se,
			Variable == "Respiration_mean" ~ Respiration_se)), 
	position = position_dodge(width = 0.9), width = 0.2) +
	
	facet_grid(Light ~ Temperature, scales = "free_y", labeller = labeller(Temperature = temperature_labels)) +
	
	geom_text(data = summary_data,
						aes(x = Species, y = max_value + 2.5,
							label = sprintf("%.2f (± %.2f)", PR_ratio_mean, PR_ratio_sd)),
						position = position_dodge(width = 0.9), inherit.aes = FALSE, vjust = -0.5, size = 3) +
	
	labs(x = NULL, y = "Oxygen Flux (nmol O₂ m⁻² s⁻¹)", fill = "Variable") +
	
	scale_fill_manual(
		values = c("netPhotosynthesis_mean" = "palegreen4", "grossPhotosynthesis_mean" = "dodgerblue2", "Respiration_mean" = "hotpink3"),
		labels = c("Net Photosynthesis", "Gross Photosynthesis", "Respiration")) +
	
	theme_bw() +
	theme(
		text = element_text(size = 12, family = "sans"),
		strip.text = element_text(size = 14, face = "bold"),
		legend.position = "top") +
	
	ylim(-15.5, 37)

if (export) {
	ggsave(filename = file.path("manuscript", "figures", "supplementary_1.png"),
				 plot = summary_plot, 
				 bg = "white",
				 width = 290, 
				 height = 230, 
				 units = "mm", 
				 dpi = 600)
}

## ========== END ==========

