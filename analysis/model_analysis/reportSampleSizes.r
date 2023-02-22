

hindcast_201511_201801 <- readRDS("./data/summary/beta_hindcast_fg_2015-11_2018-01.rds")


n_newsite_timepoints = hindcast_201511_201801 %>% filter(new_site==T & !is.na(truth)) %>%
	select(dateID) %>% distinct() %>% tally


n_newsite_plot_obs = hindcast_201511_201801 %>% filter(new_site==T & !is.na(truth)) %>%
	select(plotID, dateID) %>% distinct() %>% tally
