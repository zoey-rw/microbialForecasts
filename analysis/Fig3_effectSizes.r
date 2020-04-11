library(runjags)
library(coda)
library(forestplot)

#Summarize performance

mod.list1 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS.rds")
mod.list2 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds")

# Running for fungi
rank.list <- mod.list1
sum <- list()
for (j in 1:length(rank.list)){
    
    mod <- rank.list[[j]]
    
    df <- as.data.frame(summary(mod$JAGS, vars=c("beta_temp","beta_precip","beta_pH","beta_IC")))
    df$group <- names(rank.list)[j]
    df$cat <- "Fungi"
    sum[[j]] <- df
}
sum.all <- do.call(rbind, sum)

ci.data <-  cbind.data.frame(mean = as.numeric(sum.all[,4]), 
                             lower = as.numeric(sum.all[,1]), 
                             upper = as.numeric(sum.all[,3]))
ci.data <- rbind(c(NA,NA,NA), ci.data) # create empty top column to match text column titles.
tabletext <- cbind(c("Group", sum.all$kingdom), 
                  c("Parameter", rownames(sum.all)), 
                   c("Taxon", sum.all$group), 
                   c("EffectiveSampleSize", sum.all$SSeff))

# Extract covariate posterior estimates from summary table.
# TEMPERATURE
temp.dat <- ci.data[grep("temp",tabletext[,2]),]
temp.dat <- rbind(c(NA,NA,NA),temp.dat, c(mean(temp.dat[,1]), mean(temp.dat[,2]), mean(temp.dat[,3])))
temp.text <- tabletext[grep("temp",tabletext[,2]),]
temp.text <- rbind(c("Group", NA,"Taxon","SSeff"),
                   temp.text, c(NA,NA,"Summary",NA))
# PRECIPITATION
precip.dat <- ci.data[grep("precip",tabletext[,2]),]
precip.dat <- rbind(c(NA,NA,NA),precip.dat, c(mean(precip.dat[,1]), mean(precip.dat[,2]), mean(precip.dat[,3])))
precip.text <- tabletext[grep("precip",tabletext[,2]),]
precip.text <- rbind(c("Group", NA,"Taxon","SSeff"),
                   precip.text, c(NA, NA,"Summary",NA))
# pH
pH.dat <- ci.data[grep("pH",tabletext[,2]),]
pH.dat <- rbind(c(NA,NA,NA),pH.dat, c(mean(pH.dat[,1]), mean(pH.dat[,2]), mean(pH.dat[,3])))
pH.text <- tabletext[grep("pH",tabletext[,2]),]
pH.text <- rbind(c("Group",NA,"Taxon","SSeff"),
                     pH.text, c(NA,NA,"Summary",NA))

### CREATE FOREST PLOTS SIDE-BY-SIDE
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(temp.text[,c(3)], 
           temp.dat,
           lineheight=unit(1,'cm'),
           #graphwidth = unit(2, "inches"),
           is.summary = c(TRUE, rep(FALSE,10), TRUE),
           hrzl_lines = gpar(col="444444"),
           col=fpColors(box="royalblue",line="darkblue"),
           new_page = FALSE,
           title="Temperature effect")
popViewport()

pushViewport(viewport(layout.pos.col = 2))
forestplot(precip.text[,3], 
           precip.dat,
           lineheight=unit(1,'cm'),
           is.summary = c(TRUE, rep(FALSE, 10), TRUE),
           hrzl_lines = gpar(col="444444"),
           col=fpColors(box="royalblue",line="darkblue"),
           new_page = FALSE,
           title="Precipitation effect")
popViewport(1)

pushViewport(viewport(layout.pos.col = 3))
forestplot(pH.text[,3], 
           pH.dat,
           lineheight=unit(1,'cm'),
           is.summary = c(TRUE, rep(FALSE, 10), TRUE),
           hrzl_lines = gpar(col="444444"),
           col=fpColors(box="royalblue",line="darkblue"),
           new_page = FALSE,
           title="pH effect")
popViewport(2)

grid.text('a)', x=unit(.03, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))
grid.text('b)', x=unit(.36, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))
grid.text('c)', x=unit(.69, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))




# Now run for bacteria
rank.list <- mod.list2

sum <- list()
for (j in 1:length(rank.list)){
  mod <- rank.list[[j]]
  df <- as.data.frame(summary(mod$JAGS, vars=c("beta_temp","beta_precip","beta_pH","beta_IC")))
  df$group <- names(rank.list)[j]
  df$cat <- "Bacteria"
  sum[[j]] <- df
}
sum.all <- do.call(rbind, sum)

ci.data <-  cbind.data.frame(mean = as.numeric(sum.all[,4]), 
                             lower = as.numeric(sum.all[,1]), 
                             upper = as.numeric(sum.all[,3]))
ci.data <- rbind(c(NA,NA,NA), ci.data) # create empty top column to match text column titles.
tabletext <- cbind(c("Group", sum.all$kingdom), 
                   c("Parameter", rownames(sum.all)), 
                   c("Taxon", sum.all$group), 
                   c("EffectiveSampleSize", sum.all$SSeff))


# Extract covariate posterior estimates from summary table.
# TEMPERATURE
temp.dat <- ci.data[grep("temp",tabletext[,2]),]
temp.dat <- rbind(c(NA,NA,NA),temp.dat, c(mean(temp.dat[,1]), mean(temp.dat[,2]), mean(temp.dat[,3])))
temp.text <- tabletext[grep("temp",tabletext[,2]),]
temp.text <- rbind(c("Group", NA,"Taxon","SSeff"),
                   temp.text, c(NA,NA,"Summary",NA))
# PRECIPITATION
precip.dat <- ci.data[grep("precip",tabletext[,2]),]
precip.dat <- rbind(c(NA,NA,NA),precip.dat, c(mean(precip.dat[,1]), mean(precip.dat[,2]), mean(precip.dat[,3])))
precip.text <- tabletext[grep("precip",tabletext[,2]),]
precip.text <- rbind(c("Group", NA,"Taxon","SSeff"),
                     precip.text, c(NA, NA,"Summary",NA))
# pH
pH.dat <- ci.data[grep("pH",tabletext[,2]),]
pH.dat <- rbind(c(NA,NA,NA),pH.dat, c(mean(pH.dat[,1]), mean(pH.dat[,2]), mean(pH.dat[,3])))
pH.text <- tabletext[grep("pH",tabletext[,2]),]
pH.text <- rbind(c("Group",NA,"Taxon","SSeff"),
                 pH.text, c(NA,NA,"Summary",NA))

### CREATE FOREST PLOTS SIDE-BY-SIDE
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(temp.text[,c(3)], 
           temp.dat,
           lineheight=unit(1,'cm'),
           #graphwidth = unit(2, "inches"),
           is.summary = c(TRUE, rep(FALSE,10), TRUE),
           hrzl_lines = gpar(col="444444"),
           col=fpColors(box="royalblue",line="darkblue"),
           new_page = FALSE,
           title="Temperature effect")
popViewport()

pushViewport(viewport(layout.pos.col = 2))
forestplot(precip.text[,3], 
           precip.dat,
           lineheight=unit(1,'cm'),
           is.summary = c(TRUE, rep(FALSE, 10), TRUE),
           hrzl_lines = gpar(col="444444"),
           col=fpColors(box="royalblue",line="darkblue"),
           new_page = FALSE,
           title="Precipitation effect")
popViewport(1)

pushViewport(viewport(layout.pos.col = 3))
forestplot(pH.text[,3], 
           pH.dat,
           lineheight=unit(1,'cm'),
           is.summary = c(TRUE, rep(FALSE, 10), TRUE),
           hrzl_lines = gpar(col="444444"),
           col=fpColors(box="royalblue",line="darkblue"),
           new_page = FALSE,
           title="pH effect")
popViewport(2)


grid.text('d)', x=unit(.03, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))
grid.text('e)', x=unit(.36, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))
grid.text('f)', x=unit(.69, 'npc'), y=unit(.95, 'npc'), gp = gpar(fontsize=16))

