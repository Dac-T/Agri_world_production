library(ggplot2)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(FactoMineR) #for PCA
library(factoextra)
library(GGally)
library(cluster)

## blablabla

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
  select(-r_rain)


data = merge(data, by_country, by = c("Area", "Year"), all.x = TRUE)
data = data[, -c(9:20)]
data = data %>%
  mutate(rain.x = ifelse(is.na(rain.y), rain.x, ifelse(rain.x != rain.y, rain.y, rain.x))) %>%
  select(-rain.y)

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

# somecultnona represents all the countries where Maize, Potatoes, Rice and Sorghum are cultivated
somecultnona = na.omit(somecult)


## b) Excluding countries where some years are missing

# Extract the countries with less than 23 observations
country_counts = table(by_country$Area)
selected_countries = names(country_counts[country_counts == 23])

# Dataframes with countries that have 23 years of data
fullyeardata = data[data$Area %in% selected_countries, ]
fullyearbycount = by_country[by_country$Area %in% selected_countries, ]


## c) Combining

nona_country_counts = table(somecultnona$Area)
nona_selected_countries = names(nona_country_counts[nona_country_counts == 23])

# by_country for countries that have all the dominant crops for 23 years of data
fullscnona = somecultnona[somecultnona$Area %in% nona_selected_countries, ]



##############################

### 4. Descriptive stats


## a) numerical infos
str(data)
str(by_country)
str(fullscnona)

summary(data)
summary(by_country)
summary(fullscnona)

print(paste("Nombre de pays :", length(unique(data$Area))))
print(paste("Nombre d'années : ", length(unique(data$Year))))

sort(table(data$Area), decreasing = TRUE)
sort(table(data$Year), decreasing = TRUE)
sort(table(data$Item), decreasing = TRUE)

sort(table(by_country$Area), decreasing = TRUE)
sort(table(by_country$Year), decreasing = TRUE)

with(data, table(Item, Year))
with(data, table(Cluster, Item))


#by_country %>% group_by(Area) %>% summarise(Mean = mean(Maize), Sd = sd(Maize))

## b) Plot allowing to observe distribution, correlations and clustering
pairwiseplot = ggpairs(data, columns = c("rain", "temp", "yield", "pest"), 
                       mapping = aes(color = Cluster))


## c) Boxplots

# Boxplot for yield(crop)
mean_yield = by_country %>%
  select(-c(Year, Area)) %>%
  summarise_all(mean, na.rm = TRUE)

mean_yield = mean_yield %>%
  select(-c(Area, rain, pest, avg_temp, Cluster))

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

cropca=PCA(crop_data,scale.unit=T)

# Kaiser criterion 
cat("Nb de valeurs propres supérieures à 1 : ",  length(which(cropcaround[,1] > 1)))

#elbow + empirical mean criteria
cropcaround = round(cropca$eig,2)

fviz_eig(cropca, addlabels = TRUE, ylim = c(0, 50)) +
  geom_hline(yintercept = mean(cropcaround[,2]), linetype = "dashed", color = "red") +
  theme_minimal()

# Visualization
p1=fviz_pca_var(cropca, axes = 1:2)
p2=fviz_pca_var(cropca, axes = 3:4)
grid.arrange(p1,p2,nrow=1)

#for individuals
# ind_1=fviz_pca_ind(cropca, axes = 1:2,col.ind="cos2")
# ind_2=fviz_pca_ind(cropca, axes = 3:4,col.ind="cos2")
# grid.arrange(ind_1,ind_2,nrow=1)



# b) PCA on everything

# Nouvelle colonne contenant toutes les informations catégorielles
data$Indiv = paste(data$Area, data$Year, data$Item, data$Cluster, sep = "_")

data$Area = as.numeric((data$Area))
data$Item = as.numeric((data$Item))
data$Cluster = as.numeric(as.character(data$Cluster))

everypca = PCA(data, scale.unit = T)

indivpca = PCA(data[, c("rain", "temp", "yield", "pest")], scale.unit = T)




# yield_item = PCA(by_country[, -c("Area", "Year", "rain", "pest", "avg_temp")])
columns_to_include <- setdiff(names(by_country), c("Area", "Year", "rain", "pest", "avg_temp"))
yield_item <- PCA(na.omit(by_country[, columns_to_include]))
summary(yield_item)



##############################

# 4. Statistiques inférentielles

# Convert 'Year' and 'Item' columns to factors using as.factor
by_country$Year = as.factor(by_country$Year)
data$Year = as.factor(data$Year)

# Modèle Ancova
ancova_model = lm(yield ~ rain + temp + pest + Year, data = data)
summary(ancova_model)

par(mfrow = c(2,2)) #partitionner la fenêtre graphique en matrice carrée de dimension 2
plot(ancova_model)

# Régression linéaire multiple
regression_model = lm(yield ~ rain + temp + pest, data = data)
summary(regression_model)

par(mfrow = c(2,2)) #partitionner la fenêtre graphique en matrice carrée de dimension 2
plot(regression_model)

# Régression linéaire multiple
logreg_model = lm(log(yield) ~ rain + temp + pest, data = data)
summary(logreg_model)

par(mfrow = c(2,2)) #partitionner la fenêtre graphique en matrice carrée de dimension 2
plot(logreg_model)

# Classification ascendante hiérarchique (CAH) avec 'FactoMineR'
cah_result = HCPC(yield_acp, nb.clust = 3)
print(cah_result)




# factor_columns = sapply(bats, is.factor)
# numeric_columns = !factor_columns
# 

# 
# a = head(bats, 5)
# a %>%
#   kable() %>%
#   kable_styling(position = "left", font_size = 20)
# 
# str(bats) # a été privilégié à summary(bats)
# 
# reg1 = lm(BRW ~ BOW, data = bats_phyto)
# 
# # graphes de diagnostic
# par(mfrow = c(2,2)) #partitionner la fenêtre graphique en matrice carrée de dimension 2
# plot(reg1)
# 
# reg2 = lm(log(BRW) ~ log(BOW), data = bats_phyto)
# 
# # graphes de diagnostic
# par(mfrow = c(2,2)) #partitionner la fenêtre graphique en matrice carrée de dimension 2
# plot(reg2)
# 
# summary(reg2)
# 
# #AIC
# reg0 = lm(BRW ~ 1, data = bats_phyto)
# step(reg0, scope = ~ AUD + MOB + HIP, direction = "both")
