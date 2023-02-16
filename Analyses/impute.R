######################
## setup
######################

# setup environment
library(haven)
library(readxl)
library(tidyverse)
library(mice)
library(miceadds)
library(ggmice)
library(lme4)
`%nin%` <- Negate(`%in%`)
modus <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

######################
## data
######################

# load spss dataset
dat_raw <- read_sav("data.zsav")

# rename and clean labels
dat_rename <- janitor::clean_names(dat_raw)
vars <- names(dat_rename)

# add empty observation where necessary
dat_meting <- expand.grid(id = unique(dat_rename$id), meting = 1:5)
dat_long <- full_join(dat_rename, dat_meting)

# recode isco to oesch
oesch <- read_excel("oesch.xlsx")
dat_oesch <- dat_long %>% 
  mutate(isco08 = as.character(isco08))
dat_oesch$oesch <- plyr::mapvalues(dat_oesch$isco08, from = oesch$ISCO, to = oesch$OESCH)

# clean oesch codes
dat_clean <- dat_oesch %>% 
  select(-isco08, -weegfactor) %>% 
  mutate(across(everything(), as.numeric),
         poly = meting * meting,
         oesch = na_if(oesch, 99),
         oesch = na_if(oesch, 110),
         oesch = na_if(oesch, 310),
         oesch = na_if(oesch, 4111),
         oesch = na_if(oesch, 7536),
         .after = meting) %>% 
  rowwise() %>% 
  mutate(q43 = mean(c(q43_wave1, q43_wave2tot5), na.rm = TRUE),
         q43 = na_if(q43, "NaN"),
         .keep = "unused")

# # plot outcome
# ggplot(dat_clean, aes(meting, q60_1)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = TRUE)
# 
# ggplot(dat_clean, aes(meting, q60_1)) +
#   geom_point() +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 4), se = TRUE) +
#   theme_classic()
# ggplot(dat_clean[1:200, ], aes(meting, q60_1, color = as.factor(id))) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   theme(legend.position = "none")

# create variable vectors for imputation models
vars_outcome <- c("q60_1", "q49_1", "q49_2", "q49_3", "q49_4", "q49_5", "q49_6")
vars_fixed <- c("edu", "oesch", "geboortejaar", "q29", "q52", "q53", "q54", "sekse", "etn_herkomst")
vars_random <- c("meting", "poly", "q43", "q48_1", "q48_9", "q48_10", "q48_11", "q48_12", "q47_1", "q47_2", "q47_3", "q33_1", "q33_2", "q33_3", "q33_4", "q33_5", "q33_6")
vars_covariates <- c("q34_1", "q34_2", "q34_3", "q34_4", "q34_5", "q38_1", "q38_2", "q38_3", "q38_4", "q31_1", "q31_2", "q31_3", "q31_4", "q31_5", "q31_6", "q31_7", "q59", "q60_2", "q60_3", "q60_4", "q60_8", "q60_10", "q60_12")
vars_analyse <- c("id", vars_outcome, vars_fixed, vars_random)
vars <- names(dat_clean)

######################
## imputation
######################

# impute fixed variables
modes <- dat_clean %>% 
  group_by(id) %>% 
  summarise(across(vars_fixed, ~round(mean(.x, na.rm = TRUE), 0))) 
dat_fixed <- dat_clean[, vars %nin% vars_fixed] %>% 
  full_join(modes, by = "id")
# dat_fixed <- dat_clean %>% 
#   group_by(id) %>% 
#   mutate(across(vars_fixed, ~round(mean(.x, na.rm = TRUE), 0)))
# sum(is.na(dat_clean)) - sum(is.na(dat_fixed))

# select relevant variables
dat_imp <- dat_fixed[c(vars_analyse, vars_covariates)]
dat_imp[sapply(dat_imp, is.nan)] <- NA
vars <- names(dat_imp)

# inspect missingness
md.pattern(dat_imp[vars_fixed], rotate.names = TRUE)
md.pattern(dat_imp[vars_outcome], rotate.names = TRUE)
md.pattern(dat_imp[vars_random], rotate.names = TRUE, plot = FALSE)
ggmice::plot_corr(dat_imp[vars_analyse], label = TRUE, square = FALSE, rotate = TRUE)

# create methods
meth_imp <- mice::make.method(dat_imp)
vars_complete <- vars[colSums(is.na(dat_imp)) == 0]
meth_imp[vars %in% vars_fixed & vars %nin% vars_complete] <- "2lonly.pmm" #constante variablen
meth_imp[vars %in% vars_outcome] <- "2l.pmm"
meth_imp[vars %in% vars_random] <- "pmm"

# create pred matrix
pred_imp <- mice::quickpred(dat_imp, mincor = 0.1) #correlaties r >= 0.06
plot_pred(pred_imp)
# pred_imp[pred_imp == 1] <- 2 #vervang 1en door 2en vanwege conventies multilevelmodellen mice
pred_imp[, vars_outcome] <- 2 #uitkomst analysemodel als predictor voor alle variabelen
pred_imp[vars_outcome, ] <- 1 #alle variabelen als predictor voor uitkomst analysemodel
pred_imp[vars_random, "meting"] <- 2 #meetmoment als voorspeller voor varierende variabelen
pred_imp[vars_random, "poly"] <- 2 #meetmoment als voorspeller voor varierende variabelen
pred_imp[vars_complete, ] <- 0 #volledig geobserveerd, dus verwijderen
pred_imp[, "id"] <- -2 #clustervariabele gespecificeerd in elk model
diag(pred_imp) <- 0 #variabelen niet als predictor voor zichzelf
plot_pred(pred_imp)

# impute the data
imp_ruw <- mice(dat_imp, m = 5, maxit = 2, method = meth_imp, predictorMatrix = pred_imp, printFlag= FALSE, seed = 104, visitSequence = "monotone")
imp_it <- mice.mids(imp_ruw, maxit = 10)
saveRDS(imp_ruw, "imps.RDS")
saveRDS(imp_it, "imps_12it.RDS")
imp_it20 <- mice.mids(imp_it, maxit = 8)
saveRDS(imp_it20, "imps_20it.RDS")

######################
## factors
######################

# post-process imputaties zodat factor levels e.d. kloppen
long <- mice::complete(imp_20it, "long", include = TRUE)

# voeg kleine huishoudenscategorieen samen
long$q29[long$q29 == 6]<- 5
long$q29[long$q29 == 7]<- 5
# voeg opleidingsniveaus samen
long$edu[long$edu < 4] <- 1
long$edu[long$edu == 4| long$edu == 5] <- 2
long$edu[long$edu > 5] <- 3

# corrigeer factor levels
vars_factor <- c("edu", "oesch", "q29", "sekse", "etn_herkomst")
long_factor <- long %>% 
  mutate(
    across(vars_factor, round),
    across(vars_factor, as.factor))
levels(long_factor$sekse) <- c("Man", "Vrouw")
levels(long_factor$edu)<- c("Hoog", "Middel", "Laag")
# levels(long_factor$BG_AGE4_MIN)<- c("18-34","35-49","50-64","65+")
levels(long_factor$etn_herkomst)<- c("Nederlandse achtergrond", "Westerse migratieachtergrond", "Niet-westerse migratieachtergrond")
levels(long_factor$q29)<- c("Alleenstaand zonder kinderen", "Alleenstaand met kinderen", "Samenwonend met partner met thuiswonende kinderen",
                            "Samenwonend met partner zonder thuiswonende kinderen", "Anders")

# post-processing of scale variables
# ment = q49_1/2/3/4/5/6 #mentale gevolgen
# econ = q52/q53/q54 #economische
# onvr = q48_1/9/10/11/12 #onvrede met beleid
# drei = q47_1/2/3 #dreiging ervaren
# cont = q33_1/2/3/4/5/6 #sociaal contact
# eenz = q34_ (eenzaamheid)
# hulp = q38_ (hulp ontvangen)
# medi = q31_ (mediagebruik)
# rond = q59 (rondkomen)
# - q60_2/3/4/8/10/12 (vertrouwen)

long_scales <- long_factor %>% 
  rowwise() %>% 
  mutate(ment = mean(c(q49_1, q49_2, q49_3,  q49_4, q49_5, q49_6)),
         econ = mean(c(q52, q53, q54)),
         onvr = mean(c(q48_1, q48_9, q48_10,  q48_11, q48_12)),
         drei = mean(c(q47_1, q47_2, q47_3)),
         cont = mean(c(q33_1, q33_2, q33_3,  q33_4, q33_5, q33_6)),
         eenz = mean(c(q34_1, q34_2, q34_3,  q34_4, q34_5)),
         hulp = mean(c(q38_1, q38_2, q38_3,  q38_4)),
         medi = mean(c(q31_1, q31_2, q31_3,  q31_4, q31_5, q31_6, q31_7)), 
         .keep = "unused")

# TODO: convert age to categories
# TODO: rename q43 to gezondheid, q59 to rondkomen
# TODO: recode meting into months

# converteer terug naar mids object
imp_post <- as.mids(long_scales)
saveRDS(imp_post, "./imputaties.RDS")

######################
## analyses
######################

# fit_all <- imp_post %>% 
#   complete("all") %>% 
#   purrr::map(~lme4::lmer(q60_1 ~ . * meting - id + edu:etn_herkomst:econ + (1|id), data = .x)) 
# fit_all %>% 
#   pool() %>% 
#   broom.mixed::tidy()
# results_vert <- mitml::testEstimates(as.mitml.result(fit_all), extra.pars = TRUE)
# write.csv2(results_vert$estimates, "vertrouwen2.csv")
# saveRDS(fit_all, "vertrouwen.RDS")

fit_all <- imp %>% 
  complete("all") %>% 
  purrr::map(~lme4::lmer(q60_1 ~ . * meting - id + edu:etn_herkomst:econ + (1 + meting|id), data = .x)) 
fit_all %>% 
  pool() %>% 
  broom.mixed::tidy()
results_vert <- mitml::testEstimates(as.mitml.result(fit_all), extra.pars = TRUE)
write.csv2(results_vert$estimates, "vertrouwen_random_time_effect.csv")
saveRDS(fit_all, "vertrouwen_random_time_effect.RDS")

# fit_all <- imp_post %>% 
#   mice::complete("all") %>% 
#   purrr::map(~lme4::lmer(ment ~ . * meting - id + edu:etn_herkomst:econ + (1|id), data = .x)) 
# fit_all %>% 
#   mice::pool() %>% 
#   broom.mixed::tidy()
# results_ment <- mitml::testEstimates(as.mitml.result(fit_all), extra.pars = TRUE)
# write.csv2(results_ment$estimates, "mentaal2.csv")
# saveRDS(fit_all, "mentaal.RDS")

fit_all <- imp_post %>% 
  mice::complete("all") %>% 
  purrr::map(~lme4::lmer(ment ~ . * meting - id + edu:etn_herkomst:econ + (1 + meting |id), data = .x)) 
fit_all %>% 
  mice::pool() %>% 
  broom.mixed::tidy()
results_ment <- mitml::testEstimates(as.mitml.result(fit_all), extra.pars = TRUE)
write.csv2(results_ment$estimates, "mentaal_random_time_effects.csv")
saveRDS(fit_all, "mentaal_random_time_effect.RDS")

# # analyse the incomplete data
# mod_fixed <- lme4::lmer(q60_1 ~ . + (1|id), data = dat_imp) %>%
#   broom.mixed::tidy()
# mod_meting <- lme4::lmer(q60_1 ~ meting + (meting|id), data = dat_imp) %>%
#   broom.mixed::tidy()
# 
# # analyse the imputed data
# fit_fixed <- imp_ruw %>%
#   with(lme4::lmer(q60_1 ~ 1 + (1|id)))
# fit_fixed %>%
#   pool() %>%
#   broom.mixed::tidy()
# mitml::testEstimates(as.mitml.result(fit_fixed), extra.pars = TRUE)

# fit_meting <- imp_ruw %>% 
#   with(lme4::lmer(q60_1 ~ meting + (1|id) + (0 + meting | id))) 
# fit_meting %>% 
#   pool() %>% 
#   broom.mixed::tidy()
# mitml::testEstimates(as.mitml.result(fit_meting), extra.pars = TRUE)
# 
# fit_all <- imp_ruw %>% 
#   complete("all") %>% 
#   purrr::map(~lme4::lmer(q60_1 ~ . + (1|id), data = .x)) 
# fit_all %>% 
#   pool() %>% 
#   broom.mixed::tidy()
# mitml::testEstimates(as.mitml.result(fit_all), extra.pars = TRUE)
# 
# fit_all <- imp_ruw %>% 
#   complete("all") %>% 
#   purrr::map(~lme4::lmer(q60_1 ~ . * meting + (1|id), data = .x)) 
# fit_all %>% 
#   pool() %>% 
#   broom.mixed::tidy()
# mitml::testEstimates(as.mitml.result(fit_all), extra.pars = TRUE)