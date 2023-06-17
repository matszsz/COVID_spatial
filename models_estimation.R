library('rgdal')
library('sp')
library('RColorBrewer')
library('geosphere')
library('spdep')
library('rgdal')
library('BMS')
library('car')

# eu map
map_eu <- readOGR(".", "NUTS_RG_01M_2021_3035")
map_eu <- spTransform(map_eu, "+proj=longlat")
map_eu <- map_eu[substr(map_eu@data$NUTS_ID, 1, 3) != 'FRY', ]
map_eu <- map_eu[((map_eu@data$LEVL_CODE == 2) & 
                    (substr(map_eu@data$NUTS_ID, 1, 2) != 'DE')) | 
                   ((substr(map_eu@data$NUTS_ID, 1, 2) == 'DE') & map_eu@data$LEVL_CODE == 1),]

# loading data
data_country_specific <- read.csv('country_specific_variables.csv')
data_region_specific <- read.csv('region_specific_variables.csv')
data_region_specific$code_country <- substr(data_region_specific$code, 1, 2)
data <- base::merge(data_region_specific, data_country_specific, by = 'code_country')

# Bayesian model averaging
iters <- 1000000
set.seed(1)

data_bma_cases <- data[c('cases', 
                         'pollution', 
                         'doctors', 
                         'health_expenditure', 
                         'day', 
                         'r1b', 
                         'economic_damage_3_2', 
                         'NACE_B_E', 
                         'population_aged_60', 
                         'GDP', 
                         'density', 
                         'mean_stringency', 
                         'diabetes', 
                         'chronic_lower_respiratory_disease', 
                         'high_blood_pressure',
                         'trust_in_government',
                         'trust_in_health_authorities')]

# Table2
bma_cases <- bms(data_bma_cases, 
                 g = "UIP", 
                 mprior = "random", 
                 mprior.size = 6, 
                 mcmc = "bd", 
                 iter = iters,
                 nmodel = 0,
                 burn = 10000,
                 randomizeTimer = FALSE)

data_bma_deaths <- data[c('deaths_increase',
                          'pollution', 
                          'doctors', 
                          'health_expenditure', 
                          'day', 
                          'cases', 
                          'r1b', 
                          'economic_damage_3_2', 
                          'NACE_B_E', 
                          'population_aged_60', 
                          'GDP', 
                          'density', 
                          'mean_stringency', 
                          'diabetes', 
                          'chronic_lower_respiratory_disease',
                          'trust_in_government',
                          'trust_in_health_authorities')]

# Table 3
bma_deaths <- bms(data_bma_deaths, 
                  g = "UIP", 
                  mprior = "random", 
                  mprior.size = 6, 
                  mcmc = "bd", 
                  iter = iters,
                  nmodel = 0,
                  burn = 10000,
                  randomizeTimer = FALSE)

# drawing maps
spatial_data_plotting <- merge(y = data, x = map_eu[(map_eu@data$NUTS_ID != 'NO0B') &
                                                      (substr(map_eu@data$NUTS_ID,1,3) != 'ES7') &
                                                      (substr(map_eu@data$NUTS_ID,1,3) != 'PT3') &
                                                      (substr(map_eu@data$NUTS_ID,1,3) != 'PT2'), ], by.y = "code", by.x = "NUTS_ID")
pal <- colorRampPalette(brewer.pal(9, "YlOrRd"), bias = 1)

# Figure 1
pdf(file = "plots/number_of_cases.pdf", height = 5)
spplot(spatial_data_plotting, zcol = "cases", colorkey = list(labels = list(cex = 0.4)), col.regions = pal(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')), lwd = 0.4,
       main = FALSE)
dev.off()

# Figure 2
pdf(file = "plots/deaths_increase.pdf", height = 5)
spplot(spatial_data_plotting, zcol = "deaths_increase", colorkey = list(labels = list(cex = 0.4)), col.regions = pal(100), cuts = 99,
       par.settings = list(axis.line = list(col =  'transparent')), lwd = 0.4,
       main = FALSE)
dev.off()

# creating SpatialPolygonsDataFrame
map_eu_filtered <- map_eu[map_eu@data$NUTS_ID %in% data$code, ]
spatial_data <- merge(y = data, x = map_eu_filtered, by.y = "code", by.x = "NUTS_ID")

### creating spatial adjacency matrices
# based on shared borders
cont1 <- poly2nb(spatial_data, queen = T)
W_borders <- nb2mat(cont1, style = 'W', zero.policy = T)
W_borders_list <- mat2listw(W_borders, style = 'W')

# based on distance
distance <- distm(coordinates(spatial_data), fun = distCosine) / 1000
W_distance <- 1/distance
diag(W_distance) <- 0
W_distance <- W_distance/rowSums(W_distance)
W_distance_list <- mat2listw(W_distance, style="W")

# based on belonging to the same country
rownames(spatial_data@data) <- NULL

countries <- unique(substr(spatial_data@data$NUTS_ID, 1, 2))
neighbors <- list()
for(i in 1:length(countries)){
  neighbors[i] <- list(rownames(spatial_data@data[substr(spatial_data@data$NUTS_ID, 1, 2) == countries[i], ]))
}

names(neighbors) <- countries

W_countries <- matrix(rep(0, nrow(spatial_data)^2), nrow = nrow(spatial_data))
for(i in 1:nrow(spatial_data)){
  c <- substr(spatial_data@data$NUTS_ID[i], 1, 2)
  for(j in as.numeric(neighbors[[c]])){
    W_countries[i, j] <- 1
    W_countries[j, i] <- 1
  }
}
diag(W_countries) <- 0
for(i in 1:nrow(spatial_data)){
  W_countries[i,] <- W_countries[i,]/ifelse(sum(W_countries[i,]) == 0, 1, sum(W_countries[i,]))
}

W_countries_list <- mat2listw(W_countries, style = 'W')

### estimating models
source('spatial_models.R')

# Table 4
x_cases <- spatial_data@data[c('economic_damage_3_2', 'NACE_B_E', 'high_blood_pressure', 
                         'chronic_lower_respiratory_disease', 'population_aged_60', 
                         'day', 'GDP', 'r1b','trust_in_health_authorities', 'diabetes')]
y_cases <- spatial_data@data$cases

cases_lm <- lm(cases~economic_damage_3_2+NACE_B_E+high_blood_pressure+chronic_lower_respiratory_disease+
                 population_aged_60+day+GDP+r1b+day+trust_in_health_authorities+diabetes, data = spatial_data@data)
summary(cases_lm)
moran.test(cases_lm$residuals, W_distance_list, zero.policy = T)
moran.test(cases_lm$residuals, W_borders_list, zero.policy = T)
moran.test(cases_lm$residuals, W_countries_list, zero.policy = T)
vif(cases_lm)

# Table A6
slx_mod <- function(x, y){
  country_specific_vars <- c()
  regions_specific_vars <- c()
  
  data_country_cols <- colnames(data_country_specific)
  data_region_cols <- colnames(data_region_specific)
  
  for(col in colnames(x)){
    if(col %in% data_country_cols){
      country_specific_vars <- c(country_specific_vars, col)
    }
    if(col %in% data_region_cols){
      regions_specific_vars <- c(regions_specific_vars, col)
    }
  }
  
  m0 <- as.matrix(x)
  m1 <- W_borders%*%as.matrix(x[regions_specific_vars])
  colnames(m1) <- paste0('lag.borders.',colnames(m1))
  m2 <- W_countries%*%as.matrix(x[regions_specific_vars])
  colnames(m2) <- paste0('lag.countries.',colnames(m2))
  m3 <- W_borders%*%as.matrix(x[country_specific_vars])
  colnames(m3) <- paste0('lag.borders.',colnames(m3))
  
  dat <- cbind(y, m0, m1, m2, m3)
  model <- lm(y~., data = as.data.frame(dat))
  return(model)
}
cases_slx <- slx_mod(x_cases, y_cases)
summary(cases_slx)
moran.test(cases_slx$residuals, W_distance_list, zero.policy = T)
moran.test(cases_slx$residuals, W_borders_list, zero.policy = T)
moran.test(cases_slx$residuals, W_countries_list, zero.policy = T)

cases_mixed_W_SARAR <- spatial.mixed.W(y_cases,
                              x_cases,
                              list(W_borders_list),
                              list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SARAR)
moran.test(cases_mixed_W_SARAR$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR$residuals, W_countries_list, zero.policy = T)

cases_mixed_W_SEM <- spatial.mixed.W(y_cases,
                                   x_cases,
                                   list(),
                                   list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SEM)
moran.test(cases_mixed_W_SEM$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SEM$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SEM$residuals, W_countries_list, zero.policy = T)

# Table 5 - standardized regression coefficients
cases_mixed_W_SEM$Coef[2:length(cases_mixed_W_SEM$Coef)]*(apply(x_cases, 2, sd)/sd(y_cases))
cases_mixed_W_SARAR$Coef[2:length(cases_mixed_W_SARAR$Coef)]*(apply(x_cases, 2, sd)/sd(y_cases))

round(AIC(cases_slx), 2)
round(AIC(cases_lm), 2)
round(2*cases_mixed_W_SARAR$loglik + 2*(ncol(x_cases)+1), 2)
round(2*cases_mixed_W_SEM$loglik + 2*(ncol(x_cases)+1), 2)

cases_SARAR <- spatial.mixed.W(y_cases,
                               x_cases,
                               list(W_borders_list),
                               list(W_borders_list))
summarize(cases_SARAR)
moran.test(cases_SARAR$residuals, W_distance_list, zero.policy = T)
moran.test(cases_SARAR$residuals, W_borders_list, zero.policy = T)
moran.test(cases_SARAR$residuals, W_countries_list, zero.policy = T)

cases_SEM <- spatial.mixed.W(y_cases,
                             x_cases,
                             list(),
                             list(W_borders_list))
summarize(cases_SEM)
moran.test(cases_SEM$residuals, W_distance_list, zero.policy = T)
moran.test(cases_SEM$residuals, W_borders_list, zero.policy = T)
moran.test(cases_SEM$residuals, W_countries_list, zero.policy = T)

# Table 6
x_deaths <- spatial_data@data[c('cases', 'pollution', 'diabetes', 'mean_stringency', 
                                'day', 'doctors', 'health_expenditure')]
y_deaths <- spatial_data@data$deaths_increase

deaths_lm <- lm(deaths_increase~cases+pollution+diabetes+mean_stringency+
                  day+doctors+health_expenditure, data = spatial_data@data)
summary(deaths_lm)
moran.test(deaths_lm$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_lm$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_lm$residuals, W_countries_list, zero.policy = T)
vif(deaths_lm)

# Table A7
deaths_slx <- slx_mod(x_deaths, y_deaths)
summary(deaths_slx)
moran.test(deaths_slx$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_slx$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_slx$residuals, W_countries_list, zero.policy = T)

deaths_SEM <- spatial.mixed.W(y_deaths,
                              x_deaths,
                              list(),
                              list(W_borders_list))

summarize(deaths_SEM)
moran.test(deaths_SEM$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_SEM$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_SEM$residuals, W_countries_list, zero.policy = T)

round(qnorm(c(0.025, 0.975), deaths_SEM$Coef[2], deaths_SEM$coef_se[2]), 3)

# Table 7 - standardized regression coefficients
deaths_SEM$Coef[2:length(deaths_SEM$Coef)]*(apply(x_deaths, 2, sd)/sd(y_deaths))

round(AIC(deaths_slx), 2)
round(AIC(deaths_lm), 2)
round(2*deaths_SEM$loglik + 2*(ncol(x_deaths)+1), 2)

deaths_mixed_W_SEM <- spatial.mixed.W(y_deaths,
                                    x_deaths,
                                    list(),
                                    list(W_borders_list, W_countries_list))

summarize(deaths_mixed_W_SEM)
moran.test(deaths_mixed_W_SEM$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_mixed_W_SEM$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_mixed_W_SEM$residuals, W_countries_list, zero.policy = T)

# Table 8
x_all_cases <- spatial_data@data[c('economic_damage_3_2', 
                                   'NACE_B_E', 
                                   'high_blood_pressure', 
                                   'chronic_lower_respiratory_disease', 
                                   'population_aged_60', 
                                   'day', 
                                   'GDP', 
                                   'r1b',
                                   'trust_in_health_authorities', 
                                   'diabetes',
                                   'pollution',
                                   'mean_stringency',
                                   'doctors',
                                   'health_expenditure',
                                   'trust_in_government',
                                   'density')]

cases_mixed_W_SEM_all_vars <- spatial.mixed.W(y_cases,
                                              x_all_cases,
                                              list(),
                                              list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SEM_all_vars)
moran.test(cases_mixed_W_SEM_all_vars$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SEM_all_vars$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SEM_all_vars$residuals, W_countries_list, zero.policy = T)

cases_mixed_W_SARAR_all_vars <- spatial.mixed.W(y_cases,
                                       x_all_cases,
                                       list(W_borders_list),
                                       list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SARAR_all_vars)
moran.test(cases_mixed_W_SARAR_all_vars$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_all_vars$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_all_vars$residuals, W_countries_list, zero.policy = T)

x_all_deaths <- spatial_data@data[c('economic_damage_3_2', 
                                    'NACE_B_E', 
                                    'high_blood_pressure', 
                                    'chronic_lower_respiratory_disease', 
                                    'population_aged_60', 
                                    'day', 
                                    'GDP', 
                                    'r1b',
                                    'trust_in_health_authorities', 
                                    'diabetes',
                                    'cases',
                                    'pollution',
                                    'mean_stringency',
                                    'doctors',
                                    'health_expenditure',
                                    'trust_in_government',
                                    'density')]

deaths_SEM_all_vars <- spatial.mixed.W(y_deaths,
                              x_all_deaths,
                              list(),
                              list(W_borders_list))

summarize(deaths_SEM_all_vars)
moran.test(deaths_SEM_all_vars$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_SEM_all_vars$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_SEM_all_vars$residuals, W_countries_list, zero.policy = T)

# Table A1
x2_cases <- spatial_data@data[c('economic_damage_2_2', 'NACE_B_E', 'high_blood_pressure', 
                                'chronic_lower_respiratory_disease', 'population_aged_60', 
                                'day', 'GDP', 'r1b','trust_in_health_authorities', 'diabetes')]
x3_cases <- spatial_data@data[c('economic_damage_3_3', 'NACE_B_E', 'high_blood_pressure', 
                                'chronic_lower_respiratory_disease', 'population_aged_60', 
                                'day', 'GDP', 'r1b','trust_in_health_authorities', 'diabetes')]

cases_mixed_W_SARAR_x2 <- spatial.mixed.W(y_cases,
                                       x2_cases,
                                       list(W_borders_list),
                                       list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SARAR_x2)
moran.test(cases_mixed_W_SARAR_x2$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_x2$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_x2$residuals, W_countries_list, zero.policy = T)

cases_mixed_W_SEM_x2 <- spatial.mixed.W(y_cases,
                                     x2_cases,
                                     list(),
                                     list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SEM_x2)
moran.test(cases_mixed_W_SEM_x2$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SEM_x2$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SEM_x2$residuals, W_countries_list, zero.policy = T)

cases_mixed_W_SARAR_x3 <- spatial.mixed.W(y_cases,
                                          x3_cases,
                                          list(W_borders_list),
                                          list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SARAR_x3)
moran.test(cases_mixed_W_SARAR_x3$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_x3$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_x3$residuals, W_countries_list, zero.policy = T)

cases_mixed_W_SEM_x3 <- spatial.mixed.W(y_cases,
                                        x3_cases,
                                        list(),
                                        list(W_borders_list, W_countries_list))
summarize(cases_mixed_W_SEM_x3)
moran.test(cases_mixed_W_SEM_x3$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SEM_x3$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SEM_x3$residuals, W_countries_list, zero.policy = T)

# Table A2
cases_mixed_W_SARAR_distance <- spatial.mixed.W(y_cases,
                                       x_cases,
                                       list(W_distance_list),
                                       list(W_distance_list, W_countries_list))
summarize(cases_mixed_W_SARAR_distance)
moran.test(cases_mixed_W_SARAR_distance$residuals, W_distance_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_distance$residuals, W_borders_list, zero.policy = T)
moran.test(cases_mixed_W_SARAR_distance$residuals, W_countries_list, zero.policy = T)

deaths_mixed_W_SARAR_distance <- spatial.mixed.W(y_deaths,
                                                x_deaths,
                                                list(W_distance_list),
                                                list(W_distance_list, W_countries_list))
summarize(deaths_mixed_W_SARAR_distance)
moran.test(deaths_mixed_W_SARAR_distance$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_mixed_W_SARAR_distance$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_mixed_W_SARAR_distance$residuals, W_countries_list, zero.policy = T)

### Instrumental variable approach
set.seed(1)

# Table A3
data_bma_cases_IV <- data[c('cases', 
                         'pollution', 
                         'doctors', 
                         'health_expenditure', 
                         'day', 
                         'r1b', 
                         'economic_damage_3_2', 
                         'NACE_B_E', 
                         'population_aged_60', 
                         'GDP', 
                         'density', 
                         'mean_stringency_IV', 
                         'diabetes', 
                         'chronic_lower_respiratory_disease', 
                         'high_blood_pressure',
                         'trust_in_government',
                         'trust_in_health_authorities')]

bma_cases_IV <- bms(data_bma_cases_IV, 
                 g = "UIP", 
                 mprior = "random", 
                 mprior.size = 6, 
                 mcmc = "bd", 
                 iter = iters,
                 nmodel = 0,
                 burn = 10000,
                 randomizeTimer = FALSE)

# Table A4
data_bma_deaths_IV <- data[c('deaths_increase',
                          'pollution', 
                          'doctors', 
                          'health_expenditure', 
                          'day', 
                          'cases', 
                          'r1b', 
                          'economic_damage_3_2', 
                          'NACE_B_E', 
                          'population_aged_60', 
                          'GDP', 
                          'density', 
                          'mean_stringency_IV', 
                          'diabetes', 
                          'chronic_lower_respiratory_disease',
                          'trust_in_government',
                          'trust_in_health_authorities')]

bma_deaths_IV <- bms(data_bma_deaths_IV, 
                  g = "UIP", 
                  mprior = "random", 
                  mprior.size = 6, 
                  mcmc = "bd", 
                  iter = iters,
                  nmodel = 0,
                  burn = 10000,
                  randomizeTimer = FALSE)

x_deaths_IV <- spatial_data@data[c('cases', 'pollution', 'diabetes', 'mean_stringency_IV', 
                                   'day', 'doctors', 'health_expenditure')]
y_deaths_IV <- spatial_data@data$deaths_increase

# Table A5
deaths_SEM_IV <- spatial.mixed.W(y_deaths_IV,
                              x_deaths_IV,
                              list(),
                              list(W_borders_list))

summarize(deaths_SEM_IV)
moran.test(deaths_SEM_IV$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_SEM_IV$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_SEM_IV$residuals, W_countries_list, zero.policy = T)

# Table A8
set.seed(1)
data_bma_deaths_net <- data[c('deaths_increase',
                          'pollution', 
                          'doctors', 
                          'health_expenditure', 
                          'day',
                          'r1b', 
                          'economic_damage_3_2', 
                          'NACE_B_E', 
                          'population_aged_60', 
                          'GDP', 
                          'density', 
                          'mean_stringency', 
                          'diabetes', 
                          'chronic_lower_respiratory_disease',
                          'trust_in_government',
                          'trust_in_health_authorities')]

bma_deaths_net <- bms(data_bma_deaths_net, 
                  g = "UIP", 
                  mprior = "random", 
                  mprior.size = 6, 
                  mcmc = "bd", 
                  iter = iters,
                  nmodel = 0,
                  burn = 10000,
                  randomizeTimer = FALSE)

# Table A9
x_deaths_net <- spatial_data@data[c('day', 'pollution', 'population_aged_60', 'trust_in_government', 
                                'r1b', 'NACE_B_E', 'doctors')]
y_deaths_net <- spatial_data@data$deaths_increase

deaths_net_lm <- lm(deaths_increase~day+pollution+population_aged_60+trust_in_government
                    +r1b+NACE_B_E+doctors, data = spatial_data@data)
summary(deaths_net_lm)
moran.test(deaths_net_lm$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_net_lm$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_net_lm$residuals, W_countries_list, zero.policy = T)

deaths_net_mixed_W_SARAR <- spatial.mixed.W(y_deaths_net,
                                       x_deaths_net,
                                       list(W_borders_list),
                                       list(W_borders_list, W_countries_list))
summarize(deaths_net_mixed_W_SARAR)
moran.test(deaths_net_mixed_W_SARAR$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_net_mixed_W_SARAR$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_net_mixed_W_SARAR$residuals, W_countries_list, zero.policy = T)

deaths_net_mixed_W_SEM <- spatial.mixed.W(y_deaths_net,
                                            x_deaths_net,
                                            list(W_borders_list),
                                            list(W_borders_list))
summarize(deaths_net_mixed_W_SEM)
moran.test(deaths_net_mixed_W_SEM$residuals, W_distance_list, zero.policy = T)
moran.test(deaths_net_mixed_W_SEM$residuals, W_borders_list, zero.policy = T)
moran.test(deaths_net_mixed_W_SEM$residuals, W_countries_list, zero.policy = T)