#Code for Thesis Project
#Subject: Estimating the Observed Power of Moderator tests in Meta-Analyses
#Name: Ernani Hazbolatow

path <- "C://Thesis/Code"
addTaskCallback(function(...) {set.seed(322);TRUE})

# Packages ----------------------------------------------------------------
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library(tidyverse) # Easily Install and Load the 'Tidyverse'
library(metafor) # Meta-Analysis Package for R
library(metapower) # Power Analysis for Meta-Analysis
options(pillar.width = Inf) 

setwd(path)


# Functions ---------------------------------------------------------------
d_splitter <- function(d, cat, steps, var_equal){
  #Splits d as per E-mail
  #Steps: Defines if we split d among the categories
  #Var_equal: Nulmat gives all combinations, when vars equal, no difference
  if (steps == TRUE){
    null_vec <- c()
    steps_vec <- d/(cat - 1)
    
    for (i in 1:cat){
      null_vec[i] <- d - (d/(cat-1) * (i - 1))
    } 
    
    return(null_vec)
    
  } else if (steps == FALSE){
    null_mat <- matrix(0, nrow = cat, ncol = cat)
    diag(null_mat) <- d
    if (var_equal == TRUE){
      null_vec <- null_mat[1,]
      return(null_vec)
    } else {
      return(null_mat)
    }
  } else (
    print("Invalid Input: Needs either TRUE or FALSE for step in function(d, cat, steps)")
  )
}

pwr_bounds <- function(es_i, v_i, mod_col, alpha, tau2, d){
  max <- meta_mod(es_i, v_i, mod_col, alpha, tau2, 
                  d_splitter(d, nlevels(mod_col), FALSE, TRUE))[[8]]
  
  min <- meta_mod(es_i, v_i, mod_col, alpha, tau2, 
                  d_splitter(d, nlevels(mod_col), TRUE, TRUE))[[8]]
  pwr_bounds <- c(min, max)
  return(pwr_bounds)
}

post_hoc_contrasts <- function(yi, vi, mod_col, tau2, contrast_g, alpha, contrast_parameter){
  df <- data.frame(yi, vi, mod_col)
  
  df$'wi' <- 1/(df$vi + tau2)
  df$mod_col <- as.factor(df$mod_col)
  
  df <- df %>% filter(vi != 0)
  
  eps <- c()
  order_ <- c()
  vps <- c()
  
  print(df)
  for (i in levels(df$mod_col)){
    tempdf <- df %>% filter(df$mod_col == i)
    es_p <- sum(tempdf$yi * tempdf$wi)/sum(tempdf$wi) 
    v_p <- 1/sum(tempdf$wi)
    
    eps <- c(eps, assign(paste0("es_", i), es_p))
    vps <- c(vps, assign(paste0("v_", i), v_p))
    order_ <- c(order_ , i)
  }
  print(length(eps))
  if (sum(contrast_g) != 0 ){
    print("Contrast is not equal to 0!")
    contrast_g <- readline("Give an appropriate contrast set: ")
  } 
  
  #  if (length(eps) != length(contrast_g)){
  #    print("Contrast length is not equal to p!")
  #    contrast_g <- readline("Give an appropriate contrast set: ")
  ##  }
  
  G <- sum(eps * contrast_g)
  v_G <- sum(vps * (contrast_g^2))
  Z_G <- G/sqrt(v_G)
  
  cv_G <- qnorm(1-(alpha/2), 0, 1)
  
  if (abs(Z_G) > cv_G){
    print("Null Hypothesis (y=0) rejected")
  } else {
    print("Null Hypothesis (y=0) not rejected")
  }
  
  p_list <- c()
  for (i in contrast_parameter){
    cv_h1_scheffe <- qchisq(1-alpha, length(eps)-1) 
    power <- 1 - pnorm(sqrt(cv_h1_scheffe) - (i/sqrt(v_G))) 
    + pnorm(-sqrt(cv_h1_scheffe) - (i/sqrt(v_G)))
    p_list <- c(p_list ,c(i, power))
  }
  print(eps)
  return(p_list)
}

meta_mod <- function(yi, vi, mod_col, alpha, tau2, null){
  #a revamp of earlier functions
  #provides a series of results together 
  #creates computational ease
  
  df <- data.frame(yi, vi, mod_col)
  
  df$'wi' <- 1/(df$vi + tau2)
  df$mod_col <- as.factor(df$mod_col)
  
  df <- df %>% filter(vi != 0)
  
  eps <- c()
  order_ <- c()
  vps <- c()
  wps <- c()
  
  for (i in levels(df$mod_col)){
    tempdf <- df %>% filter(df$mod_col == i)
    es_p <- sum(tempdf$yi * tempdf$wi)/sum(tempdf$wi) 
    v_p <- 1/sum(tempdf$wi)
    #es_p <- round(es_p, 2)
    #v_p <- round(v_p, 4)
    
    eps <- c(eps, assign(paste0("es_", i), es_p))
    vps <- c(vps, assign(paste0("v_", i), v_p))
    order_ <- c(order_ , i)
  }
  
  global_es <- sum(eps * (1/vps)) / sum((1/vps))
  
  if (typeof(null) == "double"){
    w_scalar <- (1/mean(vps))
    #Assume equal variances = 1, however we must still define a scalar w (Summation identity)
    theta_ii <- mean(null)
    
    #best_est_v <- ((length(eps)-1)*sum(1/vps))/(sum(1/vps)^2 - sum((1/vps)^2))
    #best_est_w <- 1/best_est_v
    lambda_b <- w_scalar*(sum(((null - theta_ii)^2)))
    
    #print(best_est_w)
    #theta_ii <- sum(null * (1/vps)) / sum((1/vps))
    #lambda_b <- sum((1/vps) * (null - theta_ii)^2)
    crit_value <- qchisq(p = 1-alpha, 
                         df = length(eps) - 1)
    
    obs_power <- pchisq(crit_value,
                        df = length(eps) - 1,
                        ncp = lambda_b,
                        lower.tail = FALSE)
    
    return(list(eps, order_, vps, global_es, theta_ii, lambda_b, crit_value, obs_power))
    
  } else {
    QB <- sum((1/vps) * (eps - global_es)^2)
    
    crit_value <- qchisq(p = 1-alpha, 
                         df = length(eps) - 1)
    
    obs_power <- pchisq(crit_value,
                        df = length(eps) - 1,
                        ncp = QB,
                        lower.tail = FALSE)
    
    return(list(eps, order_, vps, global_es, QB, crit_value, obs_power))
  }
  
}

compute_t2 <- function(vi, i2){
  #Function that converts given I2 to T2
  #Uses Best estimator for "typical" within-study variation by Higgins (2002)
  wi <- 1/vi
  estimated_variance <- (sum(wi)*(length(wi) - 1)) / ((sum(wi))^2 - sum(wi^2))
  t2 <- (estimated_variance * i2) / (1 - i2)
  return(t2)
}


# Data Wrangling ----------------------------------------------------------------

#Data Cleaning Lindberg Table 2 2010
Lindberg2010A <- readxl::read_xlsx("Lindberg2010A.xlsx", col_names = TRUE) %>% 
  select("ES (unbiased effect size)", "SE (standard error of the effect size)", "age", "nationality", "ethnicity", "ability") %>%
  mutate(Var = `SE (standard error of the effect size)`^2)

Lindberg2010A <- Lindberg2010A %>% rename(ES = `ES (unbiased effect size)`) 

Lindberg2010A <- Lindberg2010A %>% filter(!is.na(Var))

Lindberg2010A_mod1 <- Lindberg2010A
Lindberg2010A_mod1$ability <- as.factor(Lindberg2010A_mod1$ability)
Lindberg2010A_mod1$ability <- recode_factor(Lindberg2010A_mod1$ability, 
                                            `1` = "Low", 
                                            `2` = "General",
                                            `3` = "Moderate",
                                            `4` = "Highly")

Lindberg2010A_mod2 <- Lindberg2010A
Lindberg2010A_mod2 <- Lindberg2010A_mod2 %>% filter(as.numeric(nationality) !=10, as.numeric(nationality) !=8)
Lindberg2010A_mod2$nationality  <- as.factor(Lindberg2010A_mod2$nationality)
Lindberg2010A_mod2$nationality  <- recode_factor(Lindberg2010A_mod2$nationality, 
                                                 `1` = "US", 
                                                 `2` = "Canada",
                                                 `3` = "Mexico/Caribbean & Central/South America",
                                                 `4` = "Europe",
                                                 `5` = "Australia/New Zealand",
                                                 `6` = "Asia", 
                                                 `7` = "Africa",
                                                 `9` = "Middle East")
Lindberg2010A_mod2 <- Lindberg2010A_mod2 %>% Filter(as.numeric(nationality) !=10)

Lindberg2010A_mod2 <- Lindberg2010A_mod2 %>% dplyr::filter(!(nationality %in% c("Unknown")))

Lindberg2010A_mod3 <- Lindberg2010A
Lindberg2010A_mod3 <- Lindberg2010A_mod3 %>% filter(ethnicity != 9)
Lindberg2010A_mod3 <- Lindberg2010A_mod3 %>% filter(!is.na(ethnicity))
Lindberg2010A_mod3$ethnicity  <- as.factor(Lindberg2010A_mod3$ethnicity)
Lindberg2010A_mod3$ethnicity  <- recode_factor(Lindberg2010A_mod3$ethnicity, 
                                               `1` = "Primarily Euro-American", 
                                               `2` = "Primarily minority (any group)")

Lindberg2010A_mod4 <- Lindberg2010A
Lindberg2010A_mod4 <- Lindberg2010A_mod4 %>% filter(age != 9)
Lindberg2010A_mod4$age  <- as.factor(Lindberg2010A_mod4$age)
Lindberg2010A_mod4$age  <- recode_factor(Lindberg2010A_mod4$age, 
                                         `1` = "Pre-School", 
                                         `2` = "Elementary School",
                                         `3` = "Middle School",
                                         `4` = "High School",
                                         `5` = "College",
                                         `6` = "Adult/General Population")



moderators_LA <- c("age", "ethnicity", "nationality", "ability")
meta_mod(Lindberg2010A_mod1$ES, Lindberg2010A_mod1$Var, Lindberg2010A_mod1$ability, c(0,0,0,0), 0.05, .070)
meta_mod(Lindberg2010A_mod2$ES, Lindberg2010A_mod2$Var, Lindberg2010A_mod2$nationality, c(0,0,0,0), 0.05, .070)
meta_mod(Lindberg2010A_mod3$ES, Lindberg2010A_mod3$Var, Lindberg2010A_mod3$ethnicity, c(0,0,0,0), 0.05, .070)
meta_mod(Lindberg2010A_mod4$ES, Lindberg2010A_mod4$Var, Lindberg2010A_mod4$age, c(0,0,0,0), 0.05, .070)


Lindberg2010B <- readxl::read_xlsx("Lindberg2010B.xlsx", col_names = TRUE) 

#Flynn Effect
Trahan2014 <- readxl::read_xls("Trahan2014.xls", col_names = TRUE, col_types = "text") %>%
  rename("ES" = `MeanDiffPerYear`, "var" = `v3`) %>% 
  select("ES", "Testgroup", "var", "Order")

Trahan2014$ES <- as.numeric(Trahan2014$ES)
Trahan2014$var <- as.numeric(Trahan2014$var)

Trahan2014$Testgroup <- as.factor(Trahan2014$Testgroup)

compute_t2(vi = Trahan2014$var, 0.052)


for (i in 1:length(Trahan2014$Order)){
  if (Trahan2014$Order[i] <= 25) {
    Trahan2014$Order[i] = 0
  } else if ((25 <= Trahan2014$Order[i]) & (Trahan2014$Order[i] < 75)) {
    Trahan2014$Order[i] = 50
  } else if (Trahan2014$Order[i] >= 75) {
    Trahan2014$Order[i] = 100
  }
}

Trahan2014_mod1 <- Trahan2014 %>% filter(Order == 50)

#Meta Trahan
meta_mod(Trahan2014_mod1$ES, Trahan2014_mod1$var, Trahan2014_mod1$Testgroup, 0.05, 0.8, FALSE)

boxplot(ES ~ Testgroup, data = Trahan2014_mod1)


#Retrieval paper
Murayama2014 <- readxl::read_xlsx("RIFMA_UVT.xlsx", col_names = TRUE)


Murayama2014_mod1 <- Murayama2014

## Recoding Murayama2014_mod1$design
Murayama2014_mod1$design <- Murayama2014_mod1$design %>%
  fct_recode(
    "Between" = "B",
    "Within" = "W"
  )

Murayama2014_mod1 <- Murayama2014_mod1 %>%
  filter(population == "H")



Uttal2012 <- readxl::read_xlsx("Uttal2012.xlsx", col_names = TRUE) %>%
  filter(outlier == 0) 

Uttal2012$"design" <- case_when(Uttal2012$between == 1 ~ "Between",
                                Uttal2012$within == 1 ~ "Within",
                                Uttal2012$mixed == 1 ~ "Mixed")
Uttal2012$"design" <- as.factor(Uttal2012$"design")

Uttal2012$sex <- case_when(Uttal2012$female == 1 ~ "Female",
                           Uttal2012$male == 1 ~ "Male",
                           Uttal2012$sexns == 1 ~ "No sig. improvement")
Uttal2012$sex <- as.factor(Uttal2012$sex)
Uttal2012$control <- case_when(Uttal2012$spatialfillerns == 1 ~ "No Control Group",
                               Uttal2012$spatialfiller == 1 ~ "Control Group w/ Spatial filler",
                               Uttal2012$nospatialfiller == 1 ~ "Control Group w/ different filler")
Uttal2012$control <- as.factor(Uttal2012$control)

## Recoding Uttal2012$transfer into Uttal2012$Transferable
Uttal2012$transferable <- Uttal2012$transfer %>%
  as.character() %>%
  fct_recode(
    "No transfer" = "1",
    "Within cell" = "2",
    "Across cell" = "3"
  )

Uttal2012$transferable <- as.factor(Uttal2012$transferable)

Uttal2012$training <- case_when(Uttal2012$course == 1 ~ "Course",
                                Uttal2012$vg == 1 ~ "Video games",
                                Uttal2012$spatialtask == 1 ~ "Spatial task")
Uttal2012$training <- as.factor(Uttal2012$training)
Uttal2012$durable <- case_when(Uttal2012$immediate == 1 ~ "Immediate",
                               Uttal2012$delayed == 1 ~ "Delayed")
Uttal2012$durable <- as.factor(Uttal2012$durable)

var_to_keep <- c("effect.size", "var.eff.size", "training", "durable", "control", "sex", "design", "transferable")
Uttal2012_slim <- Uttal2012 %>% select(all_of(var_to_keep))

Uttal2012_mod1 <- Uttal2012_slim %>%
  filter(!(is.na(Uttal2012_slim$durable)))

Uttal2012_mod2 <- Uttal2012_slim %>%
  filter(!(is.na(Uttal2012_slim$training)))

Uttal2012_mod3 <- Uttal2012_slim %>%
  filter(!(is.na(Uttal2012_slim$sex)))

Uttal2012_mod4 <- Uttal2012_slim %>%
  filter(!(is.na(Uttal2012_slim$design)))

Uttal2012_mod5 <- Uttal2012_slim %>%
  filter(!(is.na(Uttal2012_slim$transferable)))


for (i in 1:4){
  df <- readxl::read_xlsx("Rhodes2012.xlsx", col_names = TRUE, sheet = i)
  
  df$`Subject Age Group` <- as.factor(df$`Subject Age Group`)
  df$`Type of Encoding` <- as.factor(df$`Type of Encoding`)
  df$`Type of Encoding` <- recode_factor(df$`Type of Encoding`, 
                                         `1` = "Self-directed", 
                                         `2` = "Elaborative",
                                         `4` = "Item-specific",)
  
  assign(paste0("Rhodes2012_df", "_", i), df)
}
  

# Data Analysis -----------------------------------------------------------


## Rhodes Meta-Analysis 1 -------------

### Moderator 1: Age Group----  
table(Rhodes2012_df_1$`Subject Age Group`)
compute_t2(Rhodes2012_df_1$`Std Error (g)` ^ 2, 0.6166)

#Power with Obs. ES
obs_Rhodes2012_1_mod1 <- meta_mod(Rhodes2012_df_1$g,  (Rhodes2012_df_1$`Std Error (g)` ^ 2), 
         Rhodes2012_df_1$`Subject Age Group`, 0.05, 0.07409887, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_1_mod1_0.2 <- pwr_bounds(Rhodes2012_df_1$g, Rhodes2012_df_1$`Std Error (g)` ^ 2, Rhodes2012_df_1$`Subject Age Group`, 0.05, 0.07409887, 0.2)
bounds_Rhodes2012_1_mod1_0.5 <- pwr_bounds(Rhodes2012_df_1$g, Rhodes2012_df_1$`Std Error (g)` ^ 2,Rhodes2012_df_1$`Subject Age Group`, 0.05, 0.07409887, 0.5)

### Moderator 2: Encoding---- 
table(Rhodes2012_df_1$`Type of Encoding`)

#Power with Obs. ES
obs_Rhodes2012_1_mod2 <- meta_mod(Rhodes2012_df_1$g, (Rhodes2012_df_1$`Std Error (g)` ^ 2), 
         Rhodes2012_df_1$`Type of Encoding`, 0.05, 0.07409887, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_1_mod2_0.2 <- pwr_bounds(Rhodes2012_df_1$g, Rhodes2012_df_1$`Std Error (g)` ^ 2, Rhodes2012_df_1$`Type of Encoding`, 0.05, 0.07409887, 0.2)
bounds_Rhodes2012_1_mod2_0.5 <- pwr_bounds(Rhodes2012_df_1$g, Rhodes2012_df_1$`Std Error (g)` ^ 2, Rhodes2012_df_1$`Type of Encoding`, 0.05, 0.07409887, 0.5)



## Rhodes Meta-Analysis 2 --------------

### Moderator 1: Age Group----
table(Rhodes2012_df_2$`Subject Age Group`)
compute_t2(Rhodes2012_df_2$`Std Error (g)` ^ 2, 0.7282)

#Power with Obs. ES
obs_Rhodes2012_2_mod1 <- meta_mod(Rhodes2012_df_2$g,  (Rhodes2012_df_2$`Std Error (g)` ^ 2), 
         Rhodes2012_df_2$`Subject Age Group`, 0.05, 0.1291136, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_2_mod1_0.2 <- pwr_bounds(Rhodes2012_df_2$g, (Rhodes2012_df_2$`Std Error (g)` ^ 2), Rhodes2012_df_2$`Subject Age Group`, 0.05, 0.1291136, 0.2)
bounds_Rhodes2012_2_mod1_0.5 <- pwr_bounds(Rhodes2012_df_2$g, (Rhodes2012_df_2$`Std Error (g)` ^ 2), Rhodes2012_df_2$`Subject Age Group`, 0.05, 0.1291136, 0.5)

### Moderator 2: Encoding----
table(Rhodes2012_df_2$`Type of Encoding`)

#Power with Obs. ES
obs_Rhodes2012_2_mod2 <- meta_mod(Rhodes2012_df_2$g, (Rhodes2012_df_2$`Std Error (g)` ^ 2), 
         Rhodes2012_df_2$`Type of Encoding`, 0.05, 0.1291136, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_2_mod2_0.2 <- pwr_bounds(Rhodes2012_df_2$g, (Rhodes2012_df_2$`Std Error (g)` ^ 2),Rhodes2012_df_2$`Type of Encoding`, 0.05, 0.1291136, 0.2)
bounds_Rhodes2012_2_mod2_0.5 <- pwr_bounds(Rhodes2012_df_2$g, (Rhodes2012_df_2$`Std Error (g)` ^ 2),Rhodes2012_df_2$`Type of Encoding`, 0.05, 0.1291136, 0.5)


## Rhodes Meta-Analysis 3 ----------------

### Moderator 1: Age Group----
table(Rhodes2012_df_3$`Subject Age Group`)
compute_t2(Rhodes2012_df_3$`Std Error (g)` ^ 2, 0.7113)

#Power with Obs. ES
obs_Rhodes2012_3_mod1 <- meta_mod(Rhodes2012_df_3$g,  (Rhodes2012_df_3$`Std Error (g)` ^ 2), 
         Rhodes2012_df_3$`Subject Age Group`, 0.05, 0.1086997, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_3_mod1_0.2 <- pwr_bounds(Rhodes2012_df_3$g, (Rhodes2012_df_3$`Std Error (g)` ^ 2), Rhodes2012_df_3$`Subject Age Group`, 0.05, 0.1086997, 0.2)
bounds_Rhodes2012_3_mod1_0.5 <- pwr_bounds(Rhodes2012_df_3$g, (Rhodes2012_df_3$`Std Error (g)` ^ 2), Rhodes2012_df_3$`Subject Age Group`, 0.05, 0.1086997, 0.5)
### Moderator 2: Encoding----
table(Rhodes2012_df_3$`Type of Encoding`)
compute_t2(Rhodes2012_df_3$`Std Error (g)` ^ 2, 0.7282)

#Power with Obs. ES
obs_Rhodes2012_3_mod2 <-  meta_mod(Rhodes2012_df_3$g, (Rhodes2012_df_3$`Std Error (g)` ^ 2), 
         Rhodes2012_df_3$`Type of Encoding`, 0.05, 0.1086997, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_3_mod2_0.2 <- pwr_bounds(Rhodes2012_df_3$g, (Rhodes2012_df_3$`Std Error (g)` ^ 2), Rhodes2012_df_3$`Type of Encoding`, 0.05, 0.1086997, 0.2)
bounds_Rhodes2012_3_mod2_0.5 <- pwr_bounds(Rhodes2012_df_3$g, (Rhodes2012_df_3$`Std Error (g)` ^ 2), Rhodes2012_df_3$`Type of Encoding`, 0.05, 0.1086997, 0.5)


## Rhodes Meta-Analysis 4 ----------------

### Moderator 1: Age Group----  
table(Rhodes2012_df_4$`Subject Age Group`)
compute_t2(Rhodes2012_df_4$`Std Error (g)` ^ 2, 0.7176)

#Power with Obs. ES
obs_Rhodes2012_4_mod1 <- meta_mod(Rhodes2012_df_4$g,  (Rhodes2012_df_4$`Std Error (g)` ^ 2), 
         Rhodes2012_df_4$`Subject Age Group`, 0.05, 0.09090928, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_4_mod1_0.2 <- pwr_bounds(Rhodes2012_df_4$g, (Rhodes2012_df_4$`Std Error (g)` ^ 2), Rhodes2012_df_4$`Subject Age Group`, 0.05, 0.09090928, 0.2)
bounds_Rhodes2012_4_mod1_0.5 <- pwr_bounds(Rhodes2012_df_4$g, (Rhodes2012_df_4$`Std Error (g)` ^ 2), Rhodes2012_df_4$`Subject Age Group`, 0.05, 0.09090928, 0.5)
### Moderator 2: Encoding----
table(Rhodes2012_df_4$`Type of Encoding`)
compute_t2(Rhodes2012_df_4$`Std Error (g)` ^ 2, 0.7282)

#Power with Obs. ES
obs_Rhodes2012_4_mod2 <- meta_mod(Rhodes2012_df_4$g, (Rhodes2012_df_4$`Std Error (g)` ^ 2), 
         Rhodes2012_df_4$`Type of Encoding`, 0.05, 0.09090928, FALSE)

#Power bounds (Min-max)
bounds_Rhodes2012_4_mod2_0.2 <- pwr_bounds(Rhodes2012_df_4$g, (Rhodes2012_df_4$`Std Error (g)` ^ 2), Rhodes2012_df_4$`Type of Encoding`, 0.05, 0.09090928, 0.2)
bounds_Rhodes2012_4_mod2_0.5 <- pwr_bounds(Rhodes2012_df_4$g, (Rhodes2012_df_4$`Std Error (g)` ^ 2), Rhodes2012_df_4$`Type of Encoding`, 0.05, 0.09090928, 0.5)

## Trahan Meta-Analysis ----------------

### Moderator 1: Test Group----
table(Trahan2014_mod1$Testgroup)
compute_t2(Trahan2014_mod1$var, 0.94)

#Power with Obs. ES
obs_trahan_mod1 <- meta_mod(Trahan2014_mod1$ES, Trahan2014_mod1$var, Trahan2014_mod1$Testgroup, 0.05, 0.04764577, FALSE)

#Power bounds (Min-max)
bounds_trahan_mod1_0.2 <- pwr_bounds(Trahan2014_mod1$ES, Trahan2014_mod1$var, Trahan2014_mod1$Testgroup, 0.05, 0.04764577, 0.2)
bounds_trahan_mod1_0.5 <- pwr_bounds(Trahan2014_mod1$ES, Trahan2014_mod1$var, Trahan2014_mod1$Testgroup, 0.05, 0.04764577, 0.5)


## Uttal Meta-analysis ----------------
### Moderator 1: Durable ----
table(Uttal2012_mod1$durable)

#Power with Obs. ES
obs_uttal_mod1 <- meta_mod(Uttal2012_mod1$effect.size, Uttal2012_mod1$var.eff.size, Uttal2012_mod1$durable, 0.05, 0.185, FALSE)

#Power bounds (Min-max)
bounds_uttal_mod1_0.2 <- pwr_bounds(Uttal2012_mod1$effect.size, Uttal2012_mod1$var.eff.size, Uttal2012_mod1$durable, 0.05, 0.185, 0.2)
bounds_uttal_mod1_0.5 <- pwr_bounds(Uttal2012_mod1$effect.size, Uttal2012_mod1$var.eff.size, Uttal2012_mod1$durable, 0.05, 0.185, 0.5)

#Contrast Pwr
post_hoc_contrasts(Uttal2012_slim$effect.size, Uttal2012_slim$var.eff.size, Uttal2012_slim$durable, 0.185, c(1,-1), 0.05, c(0.2, 0.5))

### Moderator 2: Training ----
table(Uttal2012_mod2$training)

#Power with Obs. ES
obs_uttal_mod2 <- meta_mod(Uttal2012_mod2$effect.size, Uttal2012_mod2$var.eff.size, Uttal2012_mod2$training, 0.05, 0.185, FALSE)

#Power bounds (Min-max)
bounds_uttal_mod2_0.2 <- pwr_bounds(Uttal2012_mod2$effect.size, Uttal2012_mod2$var.eff.size, Uttal2012_mod2$training, 0.05, 0.185, 0.2)
bounds_uttal_mod2_0.5 <- pwr_bounds(Uttal2012_mod2$effect.size, Uttal2012_mod2$var.eff.size, Uttal2012_mod2$training, 0.05, 0.185, 0.5)

### Moderator 3: Gender ----
table(Uttal2012_mod3$sex)

#Power with Obs. ES
obs_uttal_mod3 <- meta_mod(Uttal2012_mod3$effect.size, Uttal2012_mod3$var.eff.size, Uttal2012_mod3$sex, 0.05, 0.185, FALSE)

#Power bounds (Min-max)
bounds_uttal_mod3_0.2 <- pwr_bounds(Uttal2012_mod3$effect.size, Uttal2012_mod3$var.eff.size, Uttal2012_mod3$sex, 0.05, 0.185, 0.2)
bounds_uttal_mod3_0.5 <- pwr_bounds(Uttal2012_mod3$effect.size, Uttal2012_mod3$var.eff.size, Uttal2012_mod3$sex, 0.05, 0.185, 0.5)
### Moderator 4: Design ----
table(Uttal2012_mod4$design)

#Power with Obs. ES
obs_uttal_mod4 <- meta_mod(Uttal2012_mod4$effect.size, Uttal2012_mod4$var.eff.size, Uttal2012_mod4$design, 0.05, 0.185, FALSE)

#Power bounds (Min-max)
bounds_uttal_mod4_0.2 <- pwr_bounds(Uttal2012_mod4$effect.size, Uttal2012_mod4$var.eff.size, Uttal2012_mod4$design, 0.05, 0.185, 0.2)
bounds_uttal_mod4_0.5 <- pwr_bounds(Uttal2012_mod4$effect.size, Uttal2012_mod4$var.eff.size, Uttal2012_mod4$design, 0.05, 0.185, 0.5)

### Moderator 5: Transferability ----
table(Uttal2012_mod5$transferable)

#Power with Obs. ES
obs_uttal_mod5 <-  meta_mod(Uttal2012_mod5$effect.size, Uttal2012_mod5$var.eff.size, Uttal2012_mod5$transferable, 0.05, 0.185, FALSE)
#Power bounds (Min-max)
bounds_uttal_mod5_0.2 <- pwr_bounds(Uttal2012_mod5$effect.size,  Uttal2012_mod5$var.eff.size, Uttal2012_mod5$transferable, 0.05, 0.185, 0.2)
bounds_uttal_mod5_0.5 <- pwr_bounds(Uttal2012_mod5$effect.size,  Uttal2012_mod5$var.eff.size, Uttal2012_mod5$transferable, 0.05, 0.185, 0.5)

## Lindberg 2010A Table 2 Analysis --------

### Moderator 1: Ability ----
table(Lindberg2010A_mod1$ability) 

#Power with Obs. ES
obs_Lindberg_mod1 <- meta_mod(Lindberg2010A_mod1$ES, Lindberg2010A_mod1$Var, Lindberg2010A_mod1$ability, 0.05, 0.070, null = FALSE)

#Power bounds (Min-max)
bounds_Lindberg_mod1_0.2 <- pwr_bounds(Lindberg2010A_mod1$ES, Lindberg2010A_mod1$Var,Lindberg2010A_mod1$ability, 0.05, 0.070, 0.2)
bounds_Lindberg_mod1_0.5 <- pwr_bounds(Lindberg2010A_mod1$ES, Lindberg2010A_mod1$Var,Lindberg2010A_mod1$ability, 0.05, 0.070, 0.5)
### Moderator 2: Nationality ----
table(Lindberg2010A_mod2$nationality)

#Power with Obs. ES
obs_Lindberg_mod2 <- meta_mod(Lindberg2010A_mod2$ES, Lindberg2010A_mod2$Var, Lindberg2010A_mod2$nationality, 0.05, 0.070, null = FALSE)

#Power bounds (Min-max)
bounds_Lindberg_mod2_0.2 <-  pwr_bounds(Lindberg2010A_mod2$ES,  Lindberg2010A_mod2$Var, Lindberg2010A_mod2$nationality, 0.05, 0.070, 0.2)
bounds_Lindberg_mod2_0.5 <- pwr_bounds(Lindberg2010A_mod2$ES,  Lindberg2010A_mod2$Var, Lindberg2010A_mod2$nationality, 0.05, 0.070, 0.5)

### Moderator 3: Ethnicity ----
table(Lindberg2010A_mod3$ethnicity)

#Power with Obs. ES
obs_Lindberg_mod3 <- meta_mod(Lindberg2010A_mod3$ES, Lindberg2010A_mod3$Var, Lindberg2010A_mod3$ethnicity, 0.05, 0.070, null = FALSE)

#Power bounds (Min-max)
bounds_Lindberg_mod3_0.2 <- pwr_bounds(Lindberg2010A_mod3$ES, Lindberg2010A_mod3$Var, Lindberg2010A_mod3$ethnicity, 0.05, 0.070, 0.2)
bounds_Lindberg_mod3_0.5 <- pwr_bounds(Lindberg2010A_mod3$ES, Lindberg2010A_mod3$Var, Lindberg2010A_mod3$ethnicity, 0.05, 0.070, 0.5)

#Power Contrast
post_hoc_contrasts(Lindberg2010A_mod3$ES, Lindberg2010A_mod3$Var, Lindberg2010A_mod3$ethnicity, 0.070, c(1,-1), 0.05, c(0.2, 0.5))

### Moderator 4: Age ----
table(Lindberg2010A_mod4$age)

#Power with Obs. ES
obs_Lindberg_mod4 <- meta_mod(Lindberg2010A_mod4$ES, Lindberg2010A_mod4$Var, Lindberg2010A_mod4$age, 0.05, 0.070, null = FALSE)

#Power bounds (Min-max)
bounds_Lindberg_mod4_0.2 <- pwr_bounds(Lindberg2010A_mod4$ES, Lindberg2010A_mod4$Var, Lindberg2010A_mod4$age, 0.05, 0.070, 0.2)
bounds_Lindberg_mod4_0.5 <- pwr_bounds(Lindberg2010A_mod4$ES, Lindberg2010A_mod4$Var, Lindberg2010A_mod4$age, 0.05, 0.070, 0.5)


# Plotting ----------------------------------------------------------------
obs_power <- c(obs_Rhodes2012_1_mod1[[7]], obs_Rhodes2012_1_mod2[[7]], obs_Rhodes2012_2_mod1[[7]], obs_Rhodes2012_2_mod2[[7]],
               obs_Rhodes2012_3_mod1[[7]], obs_Rhodes2012_3_mod2[[7]], obs_Rhodes2012_4_mod1[[7]], obs_Rhodes2012_4_mod2[[7]],
               obs_trahan_mod1[[7]], obs_uttal_mod1[[7]], obs_uttal_mod2[[7]], obs_uttal_mod3[[7]],
               obs_uttal_mod4[[7]], obs_uttal_mod5[[7]], obs_Lindberg_mod1[[7]], obs_Lindberg_mod2[[7]], 
               obs_Lindberg_mod3[[7]], obs_Lindberg_mod4[[7]])

bounds_0.2_min <- c(bounds_Rhodes2012_1_mod1_0.2[[1]], bounds_Rhodes2012_1_mod2_0.2[[1]], bounds_Rhodes2012_2_mod1_0.2[[1]], bounds_Rhodes2012_2_mod2_0.2[[1]],
                   bounds_Rhodes2012_3_mod1_0.2[[1]], bounds_Rhodes2012_3_mod2_0.2[[1]], bounds_Rhodes2012_4_mod1_0.2[[1]], bounds_Rhodes2012_4_mod2_0.2[[1]],
                   bounds_trahan_mod1_0.2[[1]], bounds_uttal_mod1_0.2[[1]], bounds_uttal_mod2_0.2[[1]], bounds_uttal_mod3_0.2[[1]],
                   bounds_uttal_mod4_0.2[[1]], bounds_uttal_mod5_0.2[[1]], bounds_Lindberg_mod1_0.2[[1]], bounds_Lindberg_mod2_0.2[[1]], 
                   bounds_Lindberg_mod3_0.2[[1]], bounds_Lindberg_mod4_0.2[[1]])

bounds_0.2_max <- c(bounds_Rhodes2012_1_mod1_0.2[[2]], bounds_Rhodes2012_1_mod2_0.2[[2]], bounds_Rhodes2012_2_mod1_0.2[[2]], bounds_Rhodes2012_2_mod2_0.2[[2]],
                    bounds_Rhodes2012_3_mod1_0.2[[2]], bounds_Rhodes2012_3_mod2_0.2[[2]], bounds_Rhodes2012_4_mod1_0.2[[2]], bounds_Rhodes2012_4_mod2_0.2[[2]],
                    bounds_trahan_mod1_0.2[[2]], bounds_uttal_mod1_0.2[[2]], bounds_uttal_mod2_0.2[[2]], bounds_uttal_mod3_0.2[[2]],
                    bounds_uttal_mod4_0.2[[2]], bounds_uttal_mod5_0.2[[2]], bounds_Lindberg_mod1_0.2[[2]], bounds_Lindberg_mod2_0.2[[2]], 
                    bounds_Lindberg_mod3_0.2[[2]], bounds_Lindberg_mod4_0.2[[2]])

bounds_0.5_min <- c(bounds_Rhodes2012_1_mod1_0.5[[1]], bounds_Rhodes2012_1_mod2_0.5[[1]], bounds_Rhodes2012_2_mod1_0.5[[1]], bounds_Rhodes2012_2_mod2_0.5[[1]],
                    bounds_Rhodes2012_3_mod1_0.5[[1]], bounds_Rhodes2012_3_mod2_0.5[[1]], bounds_Rhodes2012_4_mod1_0.5[[1]], bounds_Rhodes2012_4_mod2_0.5[[1]],
                    bounds_trahan_mod1_0.5[[1]], bounds_uttal_mod1_0.5[[1]], bounds_uttal_mod2_0.5[[1]], bounds_uttal_mod3_0.5[[1]],
                    bounds_uttal_mod4_0.5[[1]], bounds_uttal_mod5_0.5[[1]], bounds_Lindberg_mod1_0.5[[1]], bounds_Lindberg_mod2_0.5[[1]], 
                    bounds_Lindberg_mod3_0.5[[1]], bounds_Lindberg_mod4_0.5[[1]])
  
  
  
bounds_0.5_max <- c(bounds_Rhodes2012_1_mod1_0.5[[2]], bounds_Rhodes2012_1_mod2_0.5[[2]], bounds_Rhodes2012_2_mod1_0.5[[2]], bounds_Rhodes2012_2_mod2_0.5[[2]],
                   bounds_Rhodes2012_3_mod1_0.5[[2]], bounds_Rhodes2012_3_mod2_0.5[[2]], bounds_Rhodes2012_4_mod1_0.5[[2]], bounds_Rhodes2012_4_mod2_0.5[[2]],
                   bounds_trahan_mod1_0.5[[2]], bounds_uttal_mod1_0.5[[2]], bounds_uttal_mod2_0.5[[2]], bounds_uttal_mod3_0.5[[2]],
                   bounds_uttal_mod4_0.5[[2]], bounds_uttal_mod5_0.5[[2]], bounds_Lindberg_mod1_0.5[[2]], bounds_Lindberg_mod2_0.5[[2]], 
                   bounds_Lindberg_mod3_0.5[[2]], bounds_Lindberg_mod4_0.5[[2]])

meta_name <- c("Rhodes 1 Mod 1", "Rhodes 1 Mod 2", "Rhodes 2 Mod 1", "Rhodes 2 Mod 2",
               "Rhodes 3 Mod 1", "Rhodes 3 Mod 2", "Rhodes 4 Mod 1", "Rhodes 4 Mod 2",
               "Trahan Mod 1", "Uttal Mod 1", "Uttal Mod 2", "Uttal Mod 3", "Uttal Mod 4",
               "Uttal Mod 5", "Lindberg Mod 1", "Lindberg Mod 2", "Lindberg Mod 3", "Lindberg Mod 4")

results_df <- data.frame(obs_power, bounds_0.2_min, bounds_0.2_max, bounds_0.5_min, bounds_0.5_max, meta_name)

df.long <- pivot_longer(results_df, cols=1:5, names_to = "Type", values_to = "Estimate")

ggplot(df.long, aes(x=Estimate,y=meta_name)) +
  geom_point(aes(colour = Type, size))
