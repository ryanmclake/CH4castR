#################################################################
# CH4castR workflow                                             #
# Ryan McClure                                                  #
# COMPILE/ANALYZE/PLOT FORECASTS                                #
#################################################################
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA PATH NEEDS TO BE CHANGED to ./forecast_output/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_path <- "./forecast_output"

trap_all <- list.files(data_path, pattern = "ebu_forecast_wDA_")%>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist(fill = T)%>%
  group_by(forecast_date)%>%
  mutate(days = seq_along(time))%>%
  group_by(forecast_date)%>%
  mutate(weeks = days)%>%
  mutate(weeks = weeks-1)%>%select(-days)

trap_all_static <- list.files(data_path, pattern = "ebu_forecast_nDA_")%>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist(fill = T)%>%
  group_by(forecast_date)%>%
  mutate(days = seq_along(time))%>%
  group_by(forecast_date)%>%
  mutate(days = days-1)

trap_all_per_null <- list.files(data_path, pattern = "null_ebu_forecast_")%>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist() %>%
  group_by(forecast_date)%>%
  mutate(days = seq_along(time))%>%
  group_by(forecast_date)%>%
  mutate(days = days-1)


trap_compare <- left_join(trap_all, trap_all_static, trap_all_per_null, by = "time")

gelman_ebu_model <- list.files(data_path, pattern = "ebu_model_gelman_diagnostics_") %>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist()

gelman_ebu_model <- as.data.frame(gelman_ebu_model)


gelman <- ggplot(gelman_ebu_model, aes(forecast_date, mpsrf))+
  geom_line(size = 1.2)+
  geom_hline(yintercept = 1.005, lty = "dashed")+
  ylim(c(1,1.01))+
  theme_bw()+
  ylab("Multivariate potential scale reduction factors")+
  xlab("Forecast date")
gelman
ggsave(path = ".", filename = "./figures/SI_gelman.tiff", width = 6, height = 6, device='tiff', dpi=500)


### Parameter estimates ###
trap_all_parameters <- list.files(data_path, pattern = "ebullition_parameters") %>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist(fill = T) %>%
  group_by(forecast_date)%>%
  summarize(mean_process = mean(sd.pro, na.rm = T),
            sd_process = sd(sd.pro, na.rm = T),
            mean_intercept = mean(mu2, na.rm = T),
            sd_intercept = sd(mu2, na.rm = T),
            mean_observe = mean(phi, na.rm = T),
            sd_observe = sd(phi, na.rm = T),
            mean_temp = mean(omega, na.rm = T),
            sd_temp = sd(omega, na.rm = T))


##### UNCERTAINTY PARTITIONING DATA ####

trap_all_IC <- list.files(data_path, pattern = "model_initial_") %>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_ic = var)

trap_all_PRO <- list.files(data_path, pattern = "model_process_") %>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_pro = var)

trap_all_DRI <- list.files(data_path, pattern = "model_driver_") %>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_dri = var)

trap_all_PAR <- list.files(data_path, pattern = "model_parameter") %>%
  map(~ readRDS(file.path(data_path, .))) %>% 
  data.table::rbindlist()%>%
  select(time, forecast_date, var)%>%
  rename(var_par = var)

partition <- cbind(trap_all_IC, trap_all_PRO[,3], trap_all_DRI[,3], trap_all_PAR[,3])

all_partitioned <- partition %>%
    mutate(sum_var = var_ic+var_pro+var_dri+var_par)%>%
    mutate(`Initial condition` = var_ic/sum_var)%>%
    mutate(`Driver data` = var_dri/sum_var)%>%
    mutate(`Model parameter` = var_par/sum_var)%>%
    mutate(`Model process` = var_pro/sum_var)%>%
    mutate(sum_check = `Model process`+`Model parameter`+`Driver data`+`Initial condition`)%>%
    select(time, forecast_date, `Initial condition`, `Driver data`, `Model parameter`, `Model process`, sum_check)%>%
  mutate(days = seq_along(sum_check))%>%
  group_by(forecast_date)%>%
  mutate(day_in_future = seq_along(forecast_date))%>%
  mutate(day_in_future = day_in_future-1)%>%
  ungroup(.)

all_partitioned_melt <- all_partitioned%>%
  select(-sum_check)%>%
  melt(., id = c("time","forecast_date","days","day_in_future"))


trap_all_partition <- left_join(trap_all, full_ebullition_model_alltrap, by = c("time"))

trap_static_partition <- left_join(trap_all_static, full_ebullition_model_alltrap, by = c("time"))

trap_null_partition <- left_join(trap_all_per_null, full_ebullition_model_alltrap, by = c("time"))


one_week_forecast_wDA_eval <- trap_all_partition %>%
  select(time, mean, ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(forecast_date)%>%
  filter(row_number(forecast_date) == 2) %>%
  filter(forecast_date < as.Date("2019-11-07"))%>%
  ungroup(.)%>%
  summarize(one_week_nse_flare_model = NSE(mean, ebu_rate),
            one_week_rmse_flare_model = RMSE(mean, ebu_rate))

two_week_forecast_wDA_eval <- trap_all_partition %>%
  select(time, mean, ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(forecast_date)%>%
  filter(row_number(forecast_date) == 3) %>%
  filter(forecast_date < as.Date("2019-11-07"))%>%
  ungroup(.)%>%
  summarize(two_week_nse_flare_model = NSE(mean, ebu_rate),
            two_week_rmse_flare_model = RMSE(mean, ebu_rate))

one_week_forecast_nDA_eval <- trap_static_partition %>%
  select(time, mean, ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(forecast_date)%>%
  filter(row_number(forecast_date) == 2) %>%
  filter(forecast_date < as.Date("2019-11-07"))%>%
  ungroup(.)%>%
  summarize(one_week_nse_flare_model = NSE(mean, ebu_rate),
            one_week_rmse_flare_model = RMSE(mean, ebu_rate))

two_week_forecast_nDA_eval <- trap_static_partition %>%
  select(time, mean, ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(forecast_date)%>%
  filter(row_number(forecast_date) == 3) %>%
  filter(forecast_date < as.Date("2019-11-07"))%>%
  ungroup(.)%>%
  summarize(two_week_nse_flare_model = NSE(mean, ebu_rate),
            two_week_rmse_flare_model = RMSE(mean, ebu_rate))


one_week_forecast_null_eval <- trap_null_partition %>%
  select(time, mean, ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(forecast_date)%>%
  filter(row_number(forecast_date) == 2) %>%
  filter(forecast_date < as.Date("2019-11-07"))%>%
  ungroup(.)%>%
  summarize(one_week_nse_flare_model = NSE(mean, ebu_rate),
            one_week_rmse_flare_model = RMSE(mean, ebu_rate))

two_week_forecast_null_eval <- trap_null_partition %>%
  select(time, mean, ebu_rate, forecast_date)%>%
  na.omit(.)%>%
  group_by(forecast_date)%>%
  filter(row_number(forecast_date) == 3) %>%
  filter(forecast_date < as.Date("2019-11-07"))%>%
  ungroup(.)%>%
  summarize(two_week_nse_flare_model = NSE(mean, ebu_rate),
            two_week_rmse_flare_model = RMSE(mean, ebu_rate))

### FIGURES ###
### visualizations of the full ebullition forecasts (figure 3)

upper_y = 10
lower_y = -3

ebullition_forecasts_wDA <- trap_all %>%
  ggplot(., aes(x = time, y = mean, group = forecast_date)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "purple4", size = 1, alpha = 0.7)+
  geom_pointrange(data = full_ebullition_model_alltrap, aes(x = time, y = ebu_rate, ymin = ebu_rate-ebu_rate_se, ymax = ebu_rate+ebu_rate_se), inherit.aes = FALSE, pch = 21, color = "red", fill = "red", cex = 0.5) +
  theme_bw()+
  labs(title = "A: Forecasts refitted with new data")+
  ylab(expression(paste("log Ebullition Rate (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-25"),as.Date("2019-11-30")),
                  ylim = c(lower_y,upper_y))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

ebullition_forecasts_nDA <- trap_all_static %>%
  ggplot(., aes(x = time, y = mean, group = forecast_date)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "purple4", size = 1, alpha = 0.7)+
  geom_pointrange(data = full_ebullition_model_alltrap, aes(x = time, y = ebu_rate, ymin = ebu_rate-ebu_rate_se, ymax = ebu_rate+ebu_rate_se), inherit.aes = FALSE, pch = 21, color = "red", fill = "red", cex = 0.5) +
  theme_bw()+
  labs(title = "B: Forecasts not refitted with new data")+
  ylab(expression(paste("log Ebullition Rate (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-25"),as.Date("2019-11-30")),
                  ylim = c(lower_y,upper_y))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

null_forecasts <- trap_all_per_null %>%
  ggplot(., aes(x = time, y = mean, group = forecast_date)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "purple4", size = 1, alpha = 0.7)+
  geom_pointrange(data = full_ebullition_model_alltrap, aes(x = time, y = ebu_rate, ymin = ebu_rate-ebu_rate_se, ymax = ebu_rate+ebu_rate_se), inherit.aes = FALSE, pch = 21, color = "red", fill = "red", cex = 0.5) +
  theme_bw()+
  labs(title = "C: Persistence null forecasts")+
  ylab(expression(paste("log Ebullition Rate (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-05-25"),as.Date("2019-11-30")),
                  ylim = c(lower_y,upper_y))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

fig3 <- (ebullition_forecasts_wDA+ebullition_forecasts_nDA)/(null_forecasts+plot_spacer())
fig3
ggsave(path = ".", filename = "./figures/FIGURE3_forecasts.jpg", width = 18, height = 12, device='jpg', dpi=400)

uncertatinty <- trap_all%>% group_by(forecast_date)%>% filter(weeks > 0)%>%
  ggplot(., aes(x = weeks, y = sqrt(var), group = weeks)) +
  geom_boxplot()+
  geom_jitter(aes(color = forecast_date), width = 0.1, size = 2)+
  theme_bw()+
  labs(title = "A: Forecast Uncertainty Between One and Two-Week Forecast Horizon")+
  ylab(expression(paste("Forecast Uncertainty (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("Weeks into future")+
  scale_x_continuous(breaks = c(1,2))+
  coord_cartesian(xlim=c(0.5, 2.5))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = c(0.1, 0.8),
        legend.text = element_text(size = 16, color = "black"))


season_variance <- trap_all %>% group_by(forecast_date) %>%
  rename(`Weeks into future` = weeks)%>%filter(`Weeks into future`>0)%>%
  mutate(`Weeks into future` = as.character(`Weeks into future`))%>%
  ggplot(., aes(x = time, y = sqrt(var), group = `Weeks into future`)) +
  geom_line(aes(color = `Weeks into future`), size = 3)+
  geom_line(aes(x = time, y = sqrt(var), group = forecast_date), size = 0.3, color = "grey50")+
  theme_bw()+
  scale_color_viridis(option = "C", limits = factor(c(1,2)), 
                      breaks = c(1,2), discrete = T)+
  labs(title = "B: Forecast Uncertainty Across 2019 Forecast Period")+
  ylab(expression(paste("Forecast Uncertainty (mg CH "[4]," ",m^-2,"",d^-1,")")))+
  xlab("")+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size=15, color = "black"),
        title = element_text(size = 15),legend.position = c(0.6, 0.9),
        legend.text = element_text(size = 16, color = "black"))


fig4 <- uncertatinty/season_variance
fig4
ggsave(path = ".", filename = "./figures/FIGURE5_variance.jpg", width = 10, height = 12, device='jpg', dpi=400)


# PARAMETER ESTIAMTES FROM FORECASTS
# process <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_process)) +
#   geom_ribbon(aes(ymin = mean_process-sd_process, ymax = mean_process+sd_process), alpha = 0.2, fill = "midnightblue") +
#   geom_line(color = "black")+
#   theme_bw()+
#   labs(title = "D: Process error parameter")+
#   ylab(expression(paste(epsilon[t])))+
#   xlab("")+
#   coord_cartesian(xlim=c(as.Date("2019-06-17"),as.Date("2019-11-08")))+
#   theme(axis.text=element_text(size=15, color = "black"),
#         axis.title=element_text(size=15, color = "black"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         legend.title = element_blank(),
#         title = element_text(size = 15),legend.position = "none",
#         legend.text = element_text(size = 16, color = "black"))
# 
# int <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_intercept)) +
#   geom_ribbon(aes(ymin = mean_intercept-sd_intercept, ymax = mean_intercept+sd_intercept), alpha = 0.2, fill = "midnightblue") +
#   geom_line(color = "black")+
#   theme_bw()+
#   labs(title = "C: Intercept parameter")+
#   ylab(expression(paste(beta[1])))+
#   xlab("")+
#   coord_cartesian(xlim=c(as.Date("2019-06-17"),as.Date("2019-11-08")))+
#   theme(axis.text=element_text(size=15, color = "black"),
#         axis.title=element_text(size=15, color = "black"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         legend.title = element_blank(),
#         title = element_text(size = 15),legend.position = "none",
#         legend.text = element_text(size = 16, color = "black"))

AR <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_observe)) +
  geom_ribbon(aes(ymin = mean_observe-sd_observe, ymax = mean_observe+sd_observe), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "black")+
  theme_bw()+
  labs(title = "A: Autoregressive parameter")+
  ylab(expression(paste(beta[1])))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-06-17"),as.Date("2019-11-08")))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

temp <- ggplot(trap_all_parameters, aes(x = forecast_date, y = mean_temp)) +
  geom_ribbon(aes(ymin = mean_temp-sd_temp, ymax = mean_temp+sd_temp), alpha = 0.2, fill = "midnightblue") +
  geom_line(color = "black")+
  theme_bw()+
  labs(title = "B: Temperature parameter")+
  ylab(expression(paste(beta[2])))+
  xlab("")+
  coord_cartesian(xlim=c(as.Date("2019-06-17"),as.Date("2019-11-08")))+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15),legend.position = "none",
        legend.text = element_text(size = 16, color = "black"))

paramter = (AR+temp)
paramter
ggsave(path = ".", filename = "./figures/FIGURE4_paramters.jpg", width = 10, height = 5, device='jpg', dpi=600)



#  PARTITION UNCERTATINY
  c <- all_partitioned_melt%>%
    group_by(variable)%>%
    mutate(day_in_future = seq_along(forecast_date))%>%
    ggplot(., aes(x = day_in_future, y = value, group=interaction(forecast_date,variable), fill = interaction(forecast_date,variable)))+
           geom_area(aes(fill = variable))+
    scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9"))+
    theme_bw()+
    labs(title = "A: Partitioned uncertainty of each forecast")+
    ylab("Proportion of total variance")+
    xlab("Forecast cycle")+
    theme(axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          title = element_text(size = 15), legend.position = "top",
          legend.text = element_text(size = 10, color = "black"))
  
  e <- all_partitioned_melt %>% filter(day_in_future > 0) %>%
    ggplot(., aes(x = variable, y = value, group=variable, fill = variable))+
    geom_boxplot(aes(fill = variable))+
    scale_fill_manual(values = c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9"))+
    theme_bw()+
    labs(title = "B: Uncertainty sources across 2019 forecast season")+
    ylab("Proportion of total variance")+
    xlab("Source of uncertainty")+
    theme(axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15, color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.title = element_blank(),
          title = element_text(size = 15), legend.position = "none",
          legend.text = element_text(size = 10, color = "black"))

  
  partition = c/e
  partition
  
  ggsave(path = ".", filename = "./figures/FIGURE6_partition.jpg", width = 12, height = 12, device='jpg', dpi=400)

