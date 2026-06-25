# model testing

# v1 - as before

ses_mclustBIC1 <- lapply(ses_data_sel, function(df) {
  BIC <- NULL
  LL <- NULL
  for(j in 1:200) {
    rBIC <- mclustBIC(df[,-c(1:2)], G=1:5, verbose=FALSE,
                      initalization = list(hcPairs=hcRandomPairs(df[,-c(1:2)])))
    rModel <- mclustModel(df[,-c(1:2)], BICvalues = rBIC)
    LL <- c(LL, rModel$loglik)
    BIC <- mclustBICupdate(BIC, rBIC)
  }
  list(mod=Mclust(df[,-c(1:2)], x = BIC), LL=LL)
})

# v2 - restrict ellipsoid models with varying orientations (too many params)

ses_mclustBIC2 <- lapply(ses_data_sel, function(df) {
  BIC <- NULL
  LL <- NULL
  for(j in 1:200) {
    rBIC <- mclustBIC(df[,-c(1:2)], G=1:5, verbose=FALSE,
                      initalization = list(hcPairs=hcRandomPairs(df[,-c(1:2)])),
                      modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "XII", "XXI"))
    rModel <- mclustModel(df[,-c(1:2)], BICvalues = rBIC)
    LL <- c(LL, rModel$loglik)
    BIC <- mclustBICupdate(BIC, rBIC)
  }
  list(mod=Mclust(df[,-c(1:2)], x = BIC), LL=LL)
})

# v3 - add prior

ses_mclustBIC3 <- lapply(ses_data_sel, function(df) {
  BIC <- NULL
  LL <- NULL
  for(j in 1:200) {
    rBIC <- mclustBIC(df[,-c(1:2)], G=1:5, verbose=FALSE,
                      initalization = list(hcPairs=hcRandomPairs(df[,-c(1:2)])),
                      prior = priorControl(mean = rep(0, ncol(df[,-c(1:2)]))))
    rModel <- mclustModel(df[,-c(1:2)], BICvalues = rBIC)
    LL <- c(LL, rModel$loglik)
    BIC <- mclustBICupdate(BIC, rBIC)
  }
  list(mod=Mclust(df[,-c(1:2)], x = BIC), LL=LL)
})

# v4 - add prior AND restrict ellipsoid models with varying orientations (too many params)

ses_mclustBIC4 <- lapply(ses_data_sel, function(df) {
  BIC <- NULL
  LL <- NULL
  for(j in 1:200) {
    rBIC <- mclustBIC(df[,-c(1:2)], G=1:5, verbose=FALSE,
                      initalization = list(hcPairs=hcRandomPairs(df[,-c(1:2)])),
                      modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "XII", "XXI"),
                      prior = priorControl(mean = rep(0, ncol(df[,-c(1:2)]))))
    rModel <- mclustModel(df[,-c(1:2)], BICvalues = rBIC)
    LL <- c(LL, rModel$loglik)
    BIC <- mclustBICupdate(BIC, rBIC)
  }
  list(mod=Mclust(df[,-c(1:2)], x = BIC), LL=LL)
})

# v5 - cluster on PCA instead of scaled data

resource_pca <- prcomp(ses_data_list_sc$Resource[,-c(1:2)])
summary(resource_pca)
ses_data_pca <- list(Resource = cbind(ses_data_list_sc$Resource[,c(1:2)], resource_pca$x[,c(1:12)]))

actor_pca <- prcomp(ses_data_list_sc$Actors[,-c(1:2)])
summary(actor_pca)
ses_data_pca$Actors <- cbind(ses_data_list_sc$Actors[,c(1:2)], actor_pca$x[,c(1:10)])

ses_mclustBIC5 <- lapply(ses_data_pca, function(df) {
  BIC <- NULL
  LL <- NULL
  for(j in 1:200) {
    rBIC <- mclustBIC(df[,-c(1:2)], G=1:5, verbose=FALSE,
                      initalization = list(hcPairs=hcRandomPairs(df[,-c(1:2)])))
    rModel <- mclustModel(df[,-c(1:2)], BICvalues = rBIC)
    LL <- c(LL, rModel$loglik)
    BIC <- mclustBICupdate(BIC, rBIC)
  }
  list(mod=Mclust(df[,-c(1:2)], x = BIC), LL=LL)
})

# v6 - PCA data with prior

ses_mclustBIC6 <- lapply(ses_data_pca, function(df) {
  BIC <- NULL
  LL <- NULL
  for(j in 1:200) {
    rBIC <- mclustBIC(df[,-c(1:2)], G=1:5, verbose=FALSE,
                      initalization = list(hcPairs=hcRandomPairs(df[,-c(1:2)])),
                      prior = priorControl(mean = rep(0, ncol(df[,-c(1:2)]))))
    rModel <- mclustModel(df[,-c(1:2)], BICvalues = rBIC)
    LL <- c(LL, rModel$loglik)
    BIC <- mclustBICupdate(BIC, rBIC)
  }
  list(mod=Mclust(df[,-c(1:2)], x = BIC), LL=LL)
})

########################################
# compare versions

resource_mod_list <- list(Original = ses_mclustBIC1$Resource$mod,
                          RestrictMod = ses_mclustBIC2$Resource$mod,
                          Prior = ses_mclustBIC3$Resource$mod,
                          PriorRestrict = ses_mclustBIC4$Resource$mod,
                          PCA = ses_mclustBIC5$Resource$mod,
                          PCAPrior = ses_mclustBIC6$Resource$mod)
actors_mod_list <- list(Original = ses_mclustBIC1$Actors$mod,
                          RestrictMod = ses_mclustBIC2$Actors$mod,
                          Prior = ses_mclustBIC3$Actors$mod,
                          PriorRestrict = ses_mclustBIC4$Actors$mod,
                          PCA = ses_mclustBIC5$Actors$mod,
                          PCAPrior = ses_mclustBIC6$Actors$mod)

# num groups
lapply(resource_mod_list, function(x) {print(x$G)})
lapply(actors_mod_list, function(x) {print(x$G)})

# degrees of freedom
lapply(resource_mod_list, function(x) {print(x$df)})
lapply(actors_mod_list, function(x) {print(x$df)})

# dimension reduction plots
lapply(resource_mod_list, function(x) {plot(MclustDR(x), what="contour")})
lapply(actors_mod_list, function(x) {plot(MclustDR(x), what="contour")})

# maps
clust_dat_R <- mapply(function(df, mod, system_name) {
  left_join(ses_data_list$Resource, aoi[,c("NAME", "geometry")], by = c("county" = "NAME")) %>%
    st_as_sf() %>%
    cbind(., mod$classification) %>%
    mutate(cluster = paste0(substr("Resource",1,1), mod.classification)) %>%
    select(-mod.classification)
}, mod = resource_mod_list, SIMPLIFY = FALSE)

clust_dat_R[[4]] %>%
  ggplot(aes(fill = cluster)) +
  geom_sf() +
  facet_wrap(vars(year), nrow=1) +
  scale_fill_manual(name = NULL, values = rev(viridis::turbo(5, begin=0.075, end=0.8)),
                    labels = function(x) str_wrap(x, width = 20, whitespace_only = FALSE)) +
  labs(x=NULL, y=NULL) +
  theme_void()

clust_dat_A <- mapply(function(df, mod, system_name) {
  left_join(ses_data_list$Actors, aoi[,c("NAME", "geometry")], by = c("county" = "NAME")) %>%
    st_as_sf() %>%
    cbind(., mod$classification) %>%
    mutate(cluster = paste0(substr("Actors",1,1), mod.classification)) %>%
    select(-mod.classification)
}, mod = actors_mod_list, SIMPLIFY = FALSE)

ggplot(clust_dat_A[[4]], aes(fill = cluster)) +
  geom_sf() +
  facet_wrap(vars(year), nrow=1) +
  scale_fill_manual(name = NULL, values = rev(viridis::viridis(5)),
                    labels = function(x) str_wrap(x, width = 35, whitespace_only = FALSE)) +
  labs(x=NULL, y=NULL) +
  theme_void() +
  theme(legend.position="bottom",
        legend.key.spacing.x = unit(20, "pt"),
        legend.key.justification = "center",
        legend.key.width = unit(50, "pt"),
        text = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size=14, face="bold"),
        plot.background = element_rect(fill=NA, color="black"),
        plot.margin = margin(t=4, r=0, b=8, l=8))+
  guides(fill = guide_legend(label.position = "bottom"))

