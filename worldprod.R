library(ggplot2)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(FactoMineR) #for PCA
library(factoextra)
library(GGally)
library(cluster)
library(MASS)
library(broom)
library(nortest) #for Anderson-Darling, allows to compute normality when above 5000 observations
library(emmeans)
library(jmv) #for MANCOVA


### 1. Data importation
data = read.table("/Users/orlando/Desktop/DENS/Double_diplôme/Cours_S1/S1_D3/stats/projet_stats/Yield.csv",
                  header = TRUE, 
                  sep = "",
                  dec = ".", 
                  na.strings = c("NA", ""))


by_country = read.table("/Users/orlando/Desktop/DENS/Double_diplôme/Cours_S1/S1_D3/stats/projet_stats/Yieldbycountry.csv",
                        header = TRUE, 
                        sep = "",
                        dec = ".", 
                        na.strings = c("NA", ""))

# Add more accurate data about the rainfall, even though they are also incomplete
real_rain = read.table("/Users/orlando/Desktop/DENS/Double_diplôme/Cours_S1/S1_D3/stats/projet_stats/real_rain.csv",
                        header = TRUE, 
                        skip = 3,
                        sep = ",",
                        dec = ".", 
                        na.strings = c("NA", ""))

year_vector <- character()
year_vector = sprintf("X%d", 1990:2013)
real_rain = real_rain[, c("Country.Name", year_vector)]
names(real_rain) = c("Area", sprintf("X%d", 1990:2013))

unique_areas = unique(data$Area)
real_rain = real_rain[real_rain$Area %in% unique_areas, ]


# Reshape the dataframe using pivot_longer
real_rain <- real_rain %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "Year", 
               values_to = "r_rain", 
               values_drop_na = TRUE)
real_rain$Year = sub("X", "", real_rain$Year)

names(by_country)[names(by_country) == "average_rain_fall_mm_per_year"] = "rain"
names(by_country)[names(by_country) == "pesticides_tonnes"] = "pest"

by_country = merge(by_country, real_rain, by = c("Area", "Year"), all.x = TRUE)

by_country = by_country %>%
  mutate(rain = ifelse(is.na(r_rain), rain, ifelse(rain != r_rain, r_rain, rain))) %>%
  dplyr::select(-r_rain)


data = merge(data, by_country, by = c("Area", "Year"), all.x = TRUE)
data = data[, -c(9:20)]
data = data %>%
  mutate(rain.x = ifelse(is.na(rain.y), rain.x, ifelse(rain.x != rain.y, rain.y, rain.x))) %>%
  dplyr::select(-rain.y)

names(data)[names(data) == "rain.x"] = "rain"
names(data)[names(data) == "pest.x"] = "pest"


# Checking that there are various rain values at least for one country
unique_areas_by_country = unique(by_country$Area)
areas_with_multiple_rain_values = character(0)

for (area in unique_areas_by_country) {
  rain_values = unique(by_country$rain[by_country$Area == area])
  
  if (length(rain_values) > 1) {
    areas_with_multiple_rain_values = c(areas_with_multiple_rain_values, area)
  }
}

print(areas_with_multiple_rain_values) # --> true for 6 countries


# Factorize

data$Item = as.factor(data$Item)
data$Year = as.factor(data$Year)
data$Area = as.factor(data$Area)

by_country$Year = as.factor(by_country$Year)
by_country$Area = as.factor(by_country$Area)


##############################

### 2. Clustering on countries to reduce the modalities

numeric_data = scale(data[, c("temp", "rain")])

numeric_bycountry = scale(by_country[, c("avg_temp", "rain")])

set.seed(123)
rdmset = 40

# Function to perform k-means and return inertia
calculate_inertia = function(k, data) {
  kmeans_result = kmeans(data, centers = k, nstart = rdmset) # nstart chosen randomly
  return(kmeans_result$tot.withinss)
}

# Elbow Method
inertia_values = sapply(1:15, function(k) calculate_inertia(k, numeric_data))

elbow_plot = ggplot(data.frame(k = 1:15, inertia = inertia_values), aes(x = k, y = inertia)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Method", x = "Number of Clusters (k)", y = "Inertia")

print(elbow_plot)

inertia_values_bycountry = sapply(1:15, function(k) calculate_inertia(k, numeric_bycountry))

elbow_plot_bycountry = ggplot(data.frame(k = 1:15, inertia = inertia_values_bycountry), aes(x = k, y = inertia)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Method for by_country", x = "Number of Clusters (k)", y = "Inertia")

print(elbow_plot_bycountry)


# 6 would be the ideal amount of cluster 
k_ideal = 6

best_kmeans_result = kmeans(numeric_data, centers = k_ideal, nstart = rdmset)
bkmeans_res_bycountry = kmeans(numeric_bycountry, centers = k_ideal, nstart = rdmset)

# Access cluster assignments and centroids
best_cluster_assignments = best_kmeans_result$cluster
best_centroids = best_kmeans_result$centers

bca_bycountry = bkmeans_res_bycountry$cluster

# Access the inertia (sum of squared distances within clusters)
best_inertia = best_kmeans_result$tot.withinss
best_iw_bycountry = bkmeans_res_bycountry$tot.withinss

# Print the optimal k and other information
cat("Optimal Number of Clusters (k):", k_ideal, "\n")
cat("Centroids:", best_centroids, "\n")
cat("Inertia k-means:", best_inertia, "\n")
cat("Inertia k-means by_country:", best_iw_bycountry, "\n")

data = cbind(data, Cluster_group = best_cluster_assignments)
by_country = cbind(by_country, Cluster_group = bca_bycountry)


# Check if there are countries that are assigned to different clusters
areas_with_different_clusters = character(0)
awdc_bycountry = character(0)

for (area in unique_areas) {
  cluster_assignments = unique(data$Cluster_group[data$Area == area])
  ca_bycountry = unique(by_country$Cluster_group[by_country$Area == area])
  
  if (length(cluster_assignments) > 1) {
    areas_with_different_clusters = c(areas_with_different_clusters, area)
  }
  if (length(ca_bycountry) > 1) {
    awdc_bycountry = c(awdc_bycountry, area)
  }
}

print(areas_with_different_clusters) 
print(awdc_bycountry)


# DOESN'T WORK (SILHOUETTE_SCORES ALWAYS OVER 1)
# # Silhouette Method
# silhouette_scores = sapply(2:15, function(k) {
#   kmeans_result = kmeans(numeric_data, centers = k, nstart = 200)
#   dis = dist(numeric_data)^2
#   return(mean(silhouette(kmeans_result$cluster, dis)))
# })
# 
# best_k_silhouette <- which.max(silhouette_scores)


# # Create a data frame for ggplot
# silhouette_df <- data.frame(k = 1:15, silhouette = silhouette_scores)
# 
# # Plotting
# silhouette_plot <- ggplot(silhouette_df, aes(x = k, y = silhouette)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Silhouette Method", x = "Number of Clusters (k)", y = "Silhouette Score")
# 
# print(silhouette_plot)



# HAC
hierarchical_cluster = hclust(dist(numeric_data), method = "ward.D2")
hc_bycountry =  hclust(dist(numeric_bycountry), method = "ward.D2")

cut_dendrogram = cutree(hierarchical_cluster, k = k_ideal)  # according to the dendrogram and what we have found before, k = 6
cut_dendro_bycountry = cutree(hc_bycountry, k = k_ideal)

data = cbind(data, Cluster_HAC = cut_dendrogram)
by_country = cbind(by_country, Cluster_HAC = cut_dendro_bycountry)

# Plot the dendrogram
plot(hierarchical_cluster, main = "Hierarchical Clustering Dendrogram", xlab = "Observations", hang = -1)
rect.hclust(hierarchical_cluster, k=k_ideal)

plot(hc_bycountry, main = "Hierarchical Clustering Dendrogram by_country", xlab = "Observations", hang = -1)
rect.hclust(hc_bycountry, k=k_ideal)

# Compute the within-class inertia manually
dw = 0.5*hierarchical_cluster$height^2
iw = rev(c(0, cumsum(dw)))
cat("Inertia HAC:", iw[k_ideal], "\n")

dwbc = 0.5*hc_bycountry$height^2
iwbc = rev(c(0, cumsum(dwbc)))
cat("Inertia HAC by_country:", iwbc[k_ideal], "\n")

# check for miscategorization
areas_diff_clust_HAC = character(0)
adc_HAC_bc = character(0)

for (area in unique_areas) {
  cluster_assignments_HAC = unique(data$Cluster_HAC[data$Area == area])
  ca_HAC_bc = unique(by_country$Cluster_HAC[by_country$Area == area])
  
  if (length(cluster_assignments_HAC) > 1) {
    areas_diff_clust_HAC = c(areas_diff_clust_HAC, area)
  }
  if (length(ca_HAC_bc) > 1) {
    adc_HAC_bc = c(adc_HAC_bc, area)
  }
}

print(areas_diff_clust_HAC) # 2 errors over 101 countries
print(adc_HAC_bc)


if (best_inertia < iw[k_ideal]) {
  data = data[, -which(names(data) == "Cluster_HAC")]
  names(data)[names(data) == "Cluster_group"] = "Cluster"
} else {
  data = data[, -which(names(data) == "Cluster_group")]
  names(data)[names(data) == "Cluster_HAC"] = "Cluster"
}

if (best_iw_bycountry < iwbc[k_ideal]) {
  by_country = by_country[, -which(names(by_country) == "Cluster_HAC")]
  names(by_country)[names(by_country) == "Cluster_group"] = "Cluster"
} else {
  by_country = by_country[, -which(names(by_country) == "Cluster_group")]
  names(by_country)[names(by_country) == "Cluster_HAC"] = "Cluster"
}


# To fix factorization issues (a country associated with two clusters for instance)
data = data %>%
  group_by(Area) %>%
  mutate(Cluster = as.integer(names(which.max(table(Cluster)))))

by_country = by_country %>%
  group_by(Area) %>%
  mutate(Cluster = as.integer(names(which.max(table(Cluster)))))


#Check 
areas_with_different_clusters = character(0)
awdc_bycountry = character(0)

for (area in unique_areas) {
  cluster_assignments = unique(data$Cluster[data$Area == area])
  ca_bycountry = unique(by_country$Cluster[by_country$Area == area])
  
  if (length(cluster_assignments) > 1) {
    areas_with_different_clusters = c(areas_with_different_clusters, area)
  }
  if (length(ca_bycountry) > 1) {
    awdc_bycountry = c(awdc_bycountry, area)
  }
}

print(areas_with_different_clusters) 
print(awdc_bycountry)

# Factorize
data$Cluster = as.factor(data$Cluster)
by_country$Cluster = as.factor(by_country$Cluster)

##############################


### 3. Creating alternative dataframes

## a) Excluding crops

#Those crops are present in less than 55% of the countries
columns_to_exclude = c("Soybeans", "Cassava", "Sweet.potatoes", "Plantains.and.others", "Yams")

# Dataframe by_country with only the dominant crops
somecult = by_country[, !names(by_country) %in% columns_to_exclude]

# somecultnona represents all the countries where Maize, Wheat, Potatoes, Rice and Sorghum are cultivated
somecultnona = na.omit(somecult)


## b) Excluding countries where some years are missing

# Extract the countries with less than 23 observations
country_counts = table(by_country$Area)
selected_countries = names(country_counts[country_counts == 23])
excluded_countries = names(country_counts[country_counts != 23])



# Dataframes with countries that have 23 years of data
fullyeardata = data[data$Area %in% selected_countries, ]
fullyearbycount = by_country[by_country$Area %in% selected_countries, ]


## c) Combining

nona_country_counts = table(somecultnona$Area)
nona_selected_countries = names(nona_country_counts[nona_country_counts == 23])

# by_country for countries that have all the dominant crops for 23 years of data
fullscnona = somecultnona[somecultnona$Area %in% nona_selected_countries, ]

fullscdata = fullyeardata[!fullyeardata$Item %in% c("Soybeans", "Cassava", "Sweet potatoes", "Plantains and others", "Yams"), ]

fullscdata$Item = factor(fullscdata$Item, 
                         levels = c("Maize", "Potatoes", "Rice, paddy", "Sorghum", "Wheat"))

##############################

### 4. Descriptive stats


## a) numerical infos
str(data)
str(by_country)
str(fullscnona)
str(fullscdata)

summary(data)
summary(by_country)
summary(fullscnona)
summary(fullscdata)

print(paste("Nombre de pays :", length(unique(data$Area))))
print(paste("Nombre d'années : ", length(unique(data$Year))))

cat("Nombre de pays dans fullscnona :", length(unique(fullscnona$Area)))

sort(table(data$Area), decreasing = TRUE)
sort(table(data$Year), decreasing = TRUE)
sort(table(data$Item), decreasing = TRUE)

sort(table(by_country$Area), decreasing = TRUE)
sort(table(by_country$Year), decreasing = TRUE)

with(data, table(Item, Year))
with(data, table(Cluster, Item))


#by_country %>% group_by(Area) %>% summarise(Mean = mean(Maize), Sd = sd(Maize))

## b) Plot allowing to observe distribution, correlations and clustering

data$Year = as.numeric(as.character(data$Year))
pairwiseplot = ggpairs(data, columns = c("Year", "rain", "temp", "pest", "yield"), 
                       mapping = aes(color = Cluster))

pairwiseplot

data$Year = as.factor(data$Year)



## c) Boxplots

# Boxplot for yield(crop)
mean_yield = by_country %>%
  dplyr::select(-c(Year, Area)) %>%
  summarise_all(mean, na.rm = TRUE)

mean_yield = mean_yield %>%
  dplyr::select(-c(Area, rain, pest, avg_temp, Cluster))

# Gather the data for box plotting
gathered_data <- mean_yield %>%
  gather(key = "Crop", value = "MeanYield")

ggplot(gathered_data, aes(x = Crop, y = MeanYield, fill = Crop)) +
  geom_boxplot() +
  labs(x = "Crops", y = "Mean Yield", 
       title = "Mean Yield (over all years) Comparison for Different Crops") +
  theme_minimal()


# Boxplot for yield(cluster)
cluster_avg_yield = data %>%
  group_by(Cluster, Area) %>%
  summarise(mean_yield = mean(yield, na.rm = TRUE))

ggplot(cluster_avg_yield, aes(x = as.factor(Cluster), y = mean_yield, fill = as.factor(Cluster))) +
  geom_boxplot() +
  labs(x = "Cluster", y = "Average Yield", title = "Average Yield Comparison for Different Clusters") +
  theme_minimal() 



## d) yield(years) by crop

data$Year = as.numeric(as.character(data$Year))

# Calculate average yield and standard deviation for each crop
summary_stats <- data %>%
  group_by(Item, Year) %>%
  summarize(mean_yield = mean(yield, na.rm = TRUE),
            sd_yield = sd(yield, na.rm = TRUE))

# Create a line plot with shaded area for standard deviation
ggplot(summary_stats, aes(x = Year, y = mean_yield, color = Item, group = Item)) +
  #geom_ribbon(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield, fill = Item), alpha = 0.02) +
  geom_line((aes(y = mean_yield))) +
  labs(y = 'Average Yield') +
  scale_color_brewer(palette = "Spectral") +
  theme_minimal() 


##############################


### 5. PCA

# a) PCA on crops

crop_data = by_country[, c("Soybeans", "Cassava", "Sweet.potatoes", "Plantains.and.others", 
                            "Yams", "Maize", "Wheat", "Rice..paddy", "Sorghum", "Potatoes")]

cropca=PCA(crop_data,scale.unit=T, graph = F)

# Kaiser criterion 
cropcaround = round(cropca$eig,2)
cat("Nb de valeurs propres supérieures à 1 : ",  length(which(cropcaround[,1] > 1)))

#elbow + empirical mean criteria

fviz_eig(cropca, addlabels = TRUE, ylim = c(0, 50)) +
  geom_hline(yintercept = mean(cropcaround[,2]), linetype = "dashed", color = "red") +
  theme_minimal()

# Visualization
p1=fviz_pca_var(cropca, geom = c("text", "arrow"), col.var = "cos2", axes = 1:2)
p2=fviz_pca_var(cropca, geom = c("text", "arrow"), col.var = "cos2", axes = 3:4)
grid.arrange(p1,p2,nrow=1)

#for individuals
# ind_1=fviz_pca_ind(cropca, axes = 1:2,col.ind="cos2")
# ind_2=fviz_pca_ind(cropca, axes = 3:4,col.ind="cos2")
# grid.arrange(ind_1,ind_2,nrow=1)


# b) PCA on some crops

reduced_crop_data = somecult[, c("Maize", "Wheat", "Rice..paddy", "Sorghum", "Potatoes")]

redcropca=PCA(reduced_crop_data,scale.unit=T, graph = F)

# Kaiser criterion 
cat("Nb de valeurs propres supérieures à 1 : ",  length(which(round(redcropca$eig,2)[,1] > 1)))

#elbow + empirical mean criteria

fviz_eig(redcropca, addlabels = TRUE, ylim = c(0, 60)) +
  geom_hline(yintercept = mean(round(redcropca$eig,2)), linetype = "dashed", color = "red") +
  theme_minimal()

# Visualization
red_p1 = fviz_pca_var(redcropca, axes = 1:2, geom = c("text", "arrow"), col.var = "cos2")
red_p1

#red_ind_1= fviz_pca_ind(redcropca, axes = 1:2,col.ind="cos2")


## c) PCA on everything

data$Year = as.numeric(as.character(data$Year))
data$Cluster = as.numeric(as.character(data$Cluster))

quantindex = which(colnames(data) %in% c("Year", "Cluster"))
qualindex = which(colnames(data) %in% c("Area", "Item"))
everypca = PCA(data, scale.unit = T, ncp = 10, quanti.sup = quantindex, 
               quali.sup = qualindex, graph = F)

# Kaiser criterion 
cat("Nb de valeurs propres supérieures à 1 : ",  length(which(round(everypca$eig,2)[,1] > 1)))

#elbow + empirical mean criteria

fviz_eig(everypca, addlabels = TRUE, ylim = c(0, 60)) +
  geom_hline(yintercept = mean(round(everypca$eig,2)), linetype = "dashed", color = "red") +
  theme_minimal()

# Visualization
var = get_pca_var(everypca)
ev_p1 = fviz_pca_var(everypca, geom = c("text", "arrow"), col.var = "cos2", axes = 1:2, repel = T)
ev_p1

ev_ind_1= fviz_pca_ind(everypca, axes = 1:2, col.ind = data$Cluster)
ev_ind_1




##############################

### 4. Inferences

# Convert 'Year' and 'Item' columns to factors using as.factor
data$Year = as.factor(data$Year)
data$Cluster = as.factor(data$Cluster)

## a)	Au sein de pays ayant les mêmes conditions climatiques, la quantité de pesticide utilisée 
#     influence-t-elle le rendement des cultures, pour une culture donnée ?

# on fullscnona

# Loop for every model 
# After trial and error, the predictor variable pest and
# the dependent variable crop should be put under a log

crop_variables = c("Maize", "Wheat", "Rice..paddy", "Sorghum", "Potatoes")

models = list()

for (crop in crop_variables) {
  models[[crop]] <- list(models = list())
  
  for (N in 1:6) {
    regmod = lm(paste("log(",crop, ") ~ log(pest)"), data = fullscnona[fullscnona$Cluster == N, ])
    
    # Anderson-Darling test
    ad_stat = ad.test(residuals(regmod))$p.value # p-value < 0.05 : the sample does not come from a normal distribution
    
    # Shapiro test : 
    shapiro = shapiro.test(residuals(regmod))$p.value
    
    # Summary information
    summary_info = summary(regmod)
      
    # Store models with the best and worst Multiple R-squared

    models[[crop]]$models = c(models[[crop]]$models, list(list(N = N, ad = ad_stat, shap = shapiro,
                                                               model = regmod, summary = summary_info)))        

    }
}

# models[["X"]]$models[[N]]$model to get the reg_lin of the crop X and the cluster N 
# models[["X"]]$models[[N]]$summary for the summary, and summary$r.squared to get the R squared
# plot(models[["X"]]$models[[N]]$model) for a diagnostic graph


# Create an empty data frame to store results
result_model_df = data.frame(crop = character(), N = numeric(), r_squared = numeric(), p_value = numeric())

# Loop through crop variables and clusters
for (crop in crop_variables) {
  for (N in 1:6) {
    # Extract R-squared value
    r_squared_value = models[[crop]]$models[[N]]$summary$r.squared
    p_val = glance(models[[crop]]$models[[N]]$model)$p.value[[1]]
    
    # Append the result to the data frame
    result_model_df = rbind(result_model_df, data.frame(crop = crop, N = N, 
                                                        r_squared = r_squared_value,
                                                        p_value = p_val))
  }
}

# Find the top 5 R-squared values
top_5_models <- result_model_df[result_model_df$p_value < 0.05, ]
top_5_models <- top_5_models[order(-top_5_models$r_squared), ][1:5, ]

# Count the number of pairs with R-squared < 0.25
count_low_r_squared_or_high_p_value <- sum(result_model_df$r_squared < 0.25 | result_model_df$p_value > 0.05)

# Print results
print("Top 5 models:")
print(top_5_models)
print(paste("Number of pairs with R-squared < 0.25 or p-value > 0.05:", count_low_r_squared_or_high_p_value, "over", nrow(result_model_df)))

result_model_df$crop = as.factor(result_model_df$crop)
result_model_df$N = as.factor(result_model_df$N)


# Great graph 
slr_graph_r2 = ggplot(result_model_df, aes(x = factor(crop), y = r_squared, color = N)) +
  geom_point(size = 4) +
  labs(title = "Variation of R-squared for simple regression log(crop) ~ log(pest)",
       x = "Crops",
       y = "R-squared") +
  theme_minimal()

slr_graph_r2


# Diagnostic graph associated to the best value
par(mfrow = c(2,2))
plot(models[[top_5_models[1, "crop"]]]$models[[top_5_models[1, "N"]]]$model)
models[[top_5_models[1, "crop"]]]$models[[top_5_models[1, "N"]]]$ad
models[[top_5_models[1, "crop"]]]$models[[top_5_models[1, "N"]]]$shap

cat("AIC SLR :", extractAIC(models[[top_5_models[1, "crop"]]]$models[[top_5_models[1, "N"]]]$model))

################

## b) Comment expliquer les variations de rendement selon les variables disponibles ?
# On fullscdata - 

#### i. Multiple linear regression

mlr = lm(log(yield) ~ (rain) + log(temp) + log(pest), data = fullscdata)
ad_mlr = ad.test(residuals(mlr))

cat("Test de Anderson-Darling MLR :", ad_mlr$p.value)

par(mfrow = c(2,2)) #partitionner la fenêtre graphique en matrice carrée de dimension 2
plot(mlr)

summary(mlr)

mlr_select = step(lm(log(yield)~1, data = fullscdata), scope = ~rain + log(temp) + log(pest),
                  direction = "both")
summary(mlr_select)

cat("AIC MLR. :", extractAIC(mlr))

#### ii. General linear regression with year, Item, Cluster and pest

# MLR model
anc0 = lm(log(yield) ~ 1, data = fullscdata)

formula = log(yield) ~ Year*Item*Cluster*log(pest)
anc = lm(formula, data = fullscdata)

summary(anc)
par(mfrow = c(2,2)) 
plot(anc)
cat("Test de Anderson-Darling ANC:", ad.test(residuals(anc))$p.value)
#anova(anc0, anc)
#car::Anova(anc)
cat("AIC GLR full:", extractAIC(anc))

#Going for manual simplification (loop takes to much computational time)
nf1 = update(formula, ~ . - Year:Item:Cluster:log(pest))
anc1 = lm(nf1, data = fullscdata)
#anova(anc, anc1)
#car::Anova(anc1)
#cat("AIC GLR1 :", extractAIC(anc1))

nf2 = update(nf1, ~ . - Year:Cluster:log(pest))
anc2 = lm(nf2, data = fullscdata)
#anova(anc1, anc2)
#car::Anova(anc2)
#cat("AIC GLR2 :", extractAIC(anc2))

nf3 = update(nf2, ~ . - Year:Item:log(pest))
anc3 = lm(nf3, data = fullscdata)
#anova(anc2, anc3)
#car::Anova(anc3)
#cat("AIC GLR3 :", extractAIC(anc3))

nf4 = update(nf3, ~ . - Year:Item:Cluster)
anc4 = lm(nf4, data = fullscdata)
#anova(anc3, anc4)
#car::Anova(anc4)
#cat("AIC GLR4 :", extractAIC(anc4))

nf5 = update(nf4, ~ . - Year:Cluster)
anc5 = lm(nf5, data = fullscdata)
#anova(anc4, anc5)
#car::Anova(anc5)
#cat("AIC GLR5 :", extractAIC(anc5))

# Everything is significant for nf6
nf6 = update(nf5, ~ . - Year:Item)
anc6 = lm(nf6, data = fullscdata)
anova(anc5, anc6)
car::Anova(anc6)

summary(anc6)

par(mfrow = c(2,2)) 
plot(anc6)

cat("Test de Anderson-Darling ANC 6:", ad.test(residuals(anc6))$p.value)

cat("AIC GLR6 :", extractAIC(anc6))

# Trial without year (7) + removing the 3-interaction (8)

fwithoutyear = log(yield) ~ Item*Cluster*log(pest)
anc7 = lm(fwithoutyear, data = fullscdata)
#summary(anc7)
#car::Anova(anc7)

#cat("AIC GLR7 :", extractAIC(anc7))

nf8 = update(fwithoutyear, ~ . - Item:Cluster:log(pest))
anc8 = lm(nf8, data = fullscdata)

summary(anc8)
car::Anova(anc8)
par(mfrow = c(2,2))
plot(anc8)

cat("Test de Anderson-Darling ANC 8:", ad.test(residuals(anc8))$p.value)
cat("AIC GLR8 :", extractAIC(anc8))


# # Initialize the full model
# full_model_formula <- log(yield) ~ Year * Item * Cluster * log(pest)
# full_model <- lm(full_model_formula, data = fullscdata)
# full_anova <- car::Anova(full_model)
# 
# # Set the significance threshold
# significance_threshold <- 0.05
# 
# # Loop until no interactions or variables have p-values above the threshold
# while (any(full_anova$Pr(>F) > significance_threshold)) {
#   # Identify the variable or interaction with the highest p-value
#   max_p_value_index <- which.max(full_anova$Pr(>F))
#   max_p_value_variable <- rownames(full_anova)[max_p_value_index]
#   
#   # Remove the variable or interaction with the highest p-value
#   full_model_formula <- update(full_model_formula, . ~ . - eval(parse(text = max_p_value_variable)))
#   full_model <- lm(full_model_formula, data = fullscdata)
#   full_anova <- car::Anova(full_model)
# }
# 
# # Display the final model and ANOVA table
# summary(full_model)
# full_anova



#### iii. ANCOVA 2 factors

# New df with 3 crops compatible with parallel regression slopes

full3cdata = fullyeardata[fullyeardata$Item %in% c("Maize", "Potatoes", "Sorghum"), ]

full3cdata$Item = factor(full3cdata$Item, 
                         levels = c("Maize", "Potatoes", "Sorghum"))

full3cdata$log_yield = log(full3cdata$yield)

# Regression slope

scatter_plot = ggplot(full3cdata, aes(x = log(pest), y = log_yield, color = Item, linetype = Item)) +
  geom_point(aes(color = Item), size = 0.01) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, mapping = aes(color = Item), linetype = 1, size = 1)

scatter_plot

# Test interactions
ancov = lm(log_yield ~ log(pest) + Item + Cluster + Item*Cluster, data = full3cdata)
car::Anova(ancov)
# Some interactions are significant, there is no homogeneity in regression slopes
# but we will ignore this - just remembering this hypothesis isn't respected
# cat("AIC ANCOV :", extractAIC(ancov))

# Residuals are not normal and variances are not homogenous...
ad.test(augment(ancov)$.resid)
car::leveneTest(`.resid` ~ Item*Cluster, data = (augment(ancov)))

# but a graph says it's okay
par(mfrow = c(2,2))
plot(ancov)



# Pairwise comparison test

# Looking for the Items effect
# Even with Bonferroni correction (p must be < 8e-3), everything is statistically significant
full3cdata %>%
  group_by(Cluster) %>%
  rstatix::anova_test(log_yield ~ log(pest) + Item)

pwc_item <- full3cdata %>% 
  group_by(Cluster) %>%
  rstatix::emmeans_test(
    log_yield ~ Item, covariate = pest,
    p.adjust.method = "bonferroni")

pwc_item
emmeans(ancov, pairwise ~ Item, adjust = "bonferroni")


# Looking for the cluster's effect
# Even with Bonferroni correction (p must be < 8e-3), everything is statistically significant
full3cdata %>%
  group_by(Item) %>%
  rstatix::anova_test(log_yield ~ log(pest) + Cluster)

pwc_clust <- full3cdata %>% 
  group_by(Item) %>%
  rstatix::emmeans_test(
    log_yield ~ Cluster, covariate = pest,
    p.adjust.method = "bonferroni")

print(pwc_clust, n = 45)

emmeans(ancov, pairwise ~ Cluster, adjust = "bonferroni")


# Plot
  
pwc_item <- pwc_item %>% rstatix::add_xy_position(x = "Cluster", fun = "mean_se", step.increase = 0.05)

lp_item = ggpubr::ggline(
  rstatix::get_emmeans(pwc_item), x = "Cluster", y = "emmean", 
  color = "Item", palette = "jco") +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, color = Item), 
    width = 0.1) +
  ggpubr::stat_pvalue_manual(
    pwc_item, hide.ns = F, tip.length = 0,
    bracket.size = 0, bracket.nudge.y = - 1.4, 
    label = "{substr(group1, 1, 1)} vs. {substr(group2, 1, 1)} : {p.adj.signif}", label.size = 2.5)

pwc_clust = pwc_clust %>% add_xy_position(x = "Cluster", fun = "min", step.increase = 0.3)
pwc_clust_filtered = pwc_clust %>% dplyr::filter(p.adj.signif == "ns")
lp_clust = ggpubr::ggline(
  rstatix::get_emmeans(pwc_item), x = "Cluster", y = "emmean", 
  color = "Item", palette = "jco") +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, color = Item), 
    width = 0.1) + 
  stat_pvalue_manual(
    pwc_clust_filtered, hide.ns = FALSE, tip.length = 0.02,
    step.group.by = "Item", color = "Item", bracket.nudge.y = 4)

lp_item
lp_clust




#### iv. MANCOVA

mancova(fullscnona, deps = c("Maize", "Potatoes", "Rice..paddy", "Sorghum", "Wheat"), 
        factors = c("Cluster", "Year"), covs = c("rain", "avg_temp", "pest"),
        boxM = T,
        shapiro = T, qqPlot = T)