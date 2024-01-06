library(maps)
library(ggplot2)
library(here)
library(phyloseq)
library(usmap)
theme_set(theme_minimal())


zc <- zip_code_db
ps <- readRDS(here('Rdata','following_study','ps.rds'))
meta_data <- sample_data(ps)
us<-map_data('state')
all.states <- us_map()$full %>% unique

zip.data <- zc %>% filter(zipcode %in% meta_data$`Zip Code`)
ggplot() +
    geom_polygon(data=us,aes(x=long, y=lat, group=group),color="black", fill="gray", alpha = 0.2) +
    geom_point(data = zip.data, aes(x=lng, y=lat), color = 'red')
ggsave(here('figures','geo_plot.png'), width = 6, height = 4)    

state.info <- data.frame(meta_data) %>%
    count(State.Infor) %>% 
    data.frame()

for (i in all.states) {
    if (!(i %in% state.info[,1])){
        state.info <- state.info %>% rbind(c(i, 0))
    }
}
colnames(state.info) <- c('state', 'values')
state.info$values <- as.numeric(state.info$values)

plot_usmap(data = state.info) +
    scale_fill_continuous(low = "white", high = "blue") +
    ggtitle('Sample distribution in the US')
ggsave(here('figures','sample_distribution_state.png'), width = 6, height = 4)  
