# wili for all years-------------------------------------------
load(file =  "./Data/usflu.RData")
library(ggplot2)
library(data.table)
data$season <- ifelse(
  data$week <= 30,
  paste0(data$year - 1, "/", data$year),
  paste0(data$year, "/", data$year + 1)
)
## Season week column: week number within season
data$season_week <- sapply(seq_len(nrow(data)), function(row_ind) {
  sum(data$season == data$season[row_ind] & data$time_index <= data$time_index[row_ind])
})

data$season_week[data$season == "1997/1998"] <- data$season_week[data$season == "1997/1998"] + 9
data$InSeason <- data$week %in% c(1:20, 40:53)


data <- data.table(data)
# remove season 1997/1998 and season 2018/2019 because they are incomplete.
data <- data[ !( season %in% c("1997/1998", "2018/2019")),]
# take season 2014/2015, 2015/2016, 2016/2017 and 2017/2018 as test period.
data[, train := !(season %in% c("2014/2015", "2015/2016", "2016/2017", "2017/2018"))]

data[, group := ifelse(time <= "2010-07-25" & time >= "2008-07-27", 2,
                       ifelse(time > "2010-07-25", 3, 1) )]
data <- data.frame(data)

shade_l <- as.Date(c("1999-05-21","2000-05-20","2001-05-21","2002-05-21"))
shade_r <- as.Date(c("1999-09-24","2000-09-23","2001-09-24","2002-09-24"))

pdf(file = "./Results/plots/wILI.pdf", width = 8, height = 2)
data$group <- as.factor(data$group)
ggplot(data = data, aes(x=time, y=weighted_ili, group = group)) +
  geom_line(aes(color = group)) +
  scale_color_manual(values=c("black", "grey", "black")) +
  xlab("Year") +
  ylab("wILI (%)") +
  geom_vline(xintercept = data$time[which(data$train == FALSE)[1]], linetype="dashed", color = "black", size = 1)  +
  annotate("rect", xmin = shade_l[1], xmax = shade_r[1], ymin = -Inf, ymax = Inf,
           alpha = .2) + 
  annotate("rect", xmin = shade_l[2], xmax = shade_r[2], ymin = -Inf, ymax = Inf,
           alpha = .2) + 
  annotate("rect", xmin = shade_l[3], xmax = shade_r[3], ymin = -Inf, ymax = Inf,
           alpha = .2) + 
  annotate("rect", xmin = shade_l[4], xmax = shade_r[4], ymin = -Inf, ymax = Inf,
           alpha = .2)  +
  theme_bw() +
  theme(legend.position = "none")

dev.off()
# wili plots for each season --------------------------------------------------------------------
load(file = "./Data/data_holidays.RData")
data <- data[is.na(data$weighted_ili) == FALSE & data$season_week %in% 10 : 42, ]

seasons <- unique(data$season)
observed_peak_week <- numeric(length(seasons))
observed_peak_height <- numeric(length(seasons))

prediction_target_var <- "weighted_ili"

## get peak timing
for(i in 1 : length(seasons)) {
  analysis_time_season <- seasons[i]
  observed_peak_height[i] <- max(data[data$season == analysis_time_season, prediction_target_var])
  observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                    (data[, prediction_target_var] == observed_peak_height[i]))
  observed_peak_week[i] <- data$season_week[observed_peak_week_ind]
}
peak_points <- data.frame(time = observed_peak_week, value = observed_peak_height, season = seasons)
test_label <- data.frame(week = peak_points$time[15 : 18],
                         value = peak_points$value[15 : 18] + 0.15,
                         label = c("14/15", "15/16","16/17","17/18"))

library(ggplot2)

pdf(file = "./Results/plots/wILIholiday.pdf", width = 8, height = 4)
ggplot() +
  geom_vline(xintercept = 22, linetype="dashed", color = "gray", size = 1) +
  geom_line(data = data, aes(x=season_week, y=weighted_ili, group = season, color = train, linetype=train)) +
  scale_linetype_manual(values=c("solid", "solid")) +
  scale_color_manual(values=c( "black", "grey")) +
  xlab("Season week") +
  ylab("wILI (%)") +
  theme_bw() +
  theme(legend.position = "none") + 
  geom_point(data = peak_points, aes(x = time, y = value), shape=18, size=1) + 
  geom_text(data = test_label, aes(x = week, y = value, label = label), size = 4) 

dev.off()
# short-term scatter plots-------------------------------------------------------------------------
# source(file = "./Code/summary_ph1-4.r")
# 
# Beta_lag_res <- Beta_lag4$result
# Beta_res <- Beta$result
# KCDE_res <- kcde_prediction
# arima_res <- ARIMA$result
# 
# shade_l <- as.Date(c("2014-09-28", "2015-10-04", "2016-10-02", "2017-10-01"))
# shade_r <- as.Date( c("2015-05-17", "2016-05-15", "2017-05-14", "2018-05-13"))
# 
# library(ggplot2)
# library(colorspace)
# col_list <- qualitative_hcl(4, palette = "Dark3")
# 
# l <- length(Beta_lag_res$data_set.week)
# # in result table log score is loglik, here take difference to make that positive value indicate Beta(4) is better.
# LS_diff_T <- data.frame(diff = c( - Beta_res$log_score - (- Beta_lag_res$log_score) ,
#                                   - KCDE_res$log_score - (- Beta_lag_res$log_score),
#                                   - arima_res$log_score - (- Beta_lag_res$log_score)),
#                         week = Beta_lag_res$data_set.week,
#                         time = Beta_lag_res$data_set.time,
#                         ph = Beta_lag_res$prediction_horizon,
#                         season = Beta_lag_res$data_set.season,
#                         model = c(rep("Beta(1)",l), rep("KCDE",l), rep("ARIMA",l)))
# 
# LS_diff_T$phT <- paste("ph", LS_diff_T$ph)
# LS_diff_T$ph <- as.factor(LS_diff_T$ph)
# 
# pdf("./Results/plots/STline.pdf", height = 6, width = 8)
# p <- ggplot(LS_diff_T, aes(x=time, y=diff, color = ph, linetype = ph)) + 
#   geom_line() +
#   scale_color_manual(values = col_list) +
#   scale_linetype_manual(values=c("solid", "dashed","twodash", "longdash")) +
#   xlab("Time") +
#   ylab("Log score difference") +
#   annotate("rect", xmin = shade_l[1], xmax = shade_r[1], ymin = -Inf, ymax = Inf,
#            alpha = .2) + 
#   annotate("rect", xmin = shade_l[2], xmax = shade_r[2], ymin = -Inf, ymax = Inf,
#            alpha = .2) + 
#   annotate("rect", xmin = shade_l[3], xmax = shade_r[3], ymin = -Inf, ymax = Inf,
#            alpha = .2) + 
#   annotate("rect", xmin = shade_l[4], xmax = shade_r[4], ymin = -Inf, ymax = Inf,
#            alpha = .2) +
#   theme_bw()
# p + facet_wrap( ~  model, scales="free_x", shrink = TRUE, ncol = 1)
# dev.off()
# short term box plot --------------------------------------------------------------------------
library(ggplot2)
library(colorspace)
col_list <- sequential_hcl(4, palette = "Sunset")
pdf(file = "./Results/plots/STbox.pdf", height = 6, width = 8)
p <- ggplot(data = LS_diff_T[ LS_diff_T$week %in% c(40:53, 1:20),], aes(x=model, y=diff, fill = season)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(fill = season), outlier.shape = NA) +
  coord_cartesian(ylim = c(-1.5, 1.5))+
  scale_fill_manual(values = col_list) + 
  stat_summary(fun.y = mean, geom="point", position=position_dodge(width=0.75), shape=18) 

p + facet_wrap( ~  phT, scales="free_x", shrink = TRUE) + xlab("Model and prediction horizon") +
  ylab("Log score difference") + guides(fill=guide_legend(title="Season")) +
  theme_bw()

dev.off()

# peak box plot before peak --------------------------------------------------------------------------------
peak_arima <- readRDS(file = "./Results/Peak/peak-week-arima.rds")
peak_sarima <- readRDS(file = "./Results/Peak/peak-week-sarima.rds")
peak_Beta <- readRDS(file = "./Results/Peak/peak-week-beta3-4.rds")
peak_Beta_lags <- readRDS(file = "./Results/Peak/peak-week-beta-lags.rds")
peak_prophet <- readRDS(file = "./Results/Peak/peak-week-prophet.rds")
peak_naive <- readRDS(file = "./Results/Peak/peak-week-naive.rds")
peak_KCDE <- readRDS(file = "./Results/Peak/peak-week-KCDE.rds")

result_table <- list(peak_Beta,
                     peak_Beta_lags,
                     peak_KCDE,
                     peak_arima,
                     peak_sarima,
                     peak_prophet,
                     peak_naive)


model_name <- c("Beta(1)",
                "Beta(4)",
                "KCDE",
                "ARIMA",
                "SARIMA",
                "Prophet",
                "Naive")


result_table0 <- list()
for (i in 1 : length(result_table)){
  result_table0[[i]] <- result_table[[i]][,c("analysis_time_season",
                                             "analysis_time_season_week",
                                             "peak_week_log_score",
                                             "peak_height_log_score")]
}

load(file = "./Data/data_holidays.RData")

analysis_seasons <- c("2014/2015","2015/2016", "2016/2017", "2017/2018")
observed_peak_week <- numeric(length(analysis_seasons))
observed_peak_height <- numeric(length(analysis_seasons))

first_analysis_time_season_week <- 10 # == week 40 of year
last_analysis_time_season_week <- 41 # analysis for 33-week season, consistent with flu competition -- at week 41, we do prediction for a horizon of one week ahead

prediction_target_var <- "weighted_ili"

for(i in 1 : length(analysis_seasons)) {
  analysis_time_season <- analysis_seasons[i]
  ## get observed quantities, for computing log score
  observed_peak_height[i] <- max(data[data$season == analysis_time_season, prediction_target_var])
  observed_peak_week_ind <- which((data$season == analysis_time_season) &
                                    (data[, prediction_target_var] == observed_peak_height[i]))
  observed_peak_week[i] <- data$season_week[observed_peak_week_ind]
}

season_list <- character()
LS_timing <- numeric()
LS_int <- numeric() 
model_list <- numeric()
Diff_row_idx <- 1

for(i in 1: length(model_name)){
  target <- result_table0[[i]]
  
  for(j in 1 : length(analysis_seasons)){
    season <- analysis_seasons[j]
    BP_idx <- target$analysis_time_season == season & (target$analysis_time_season_week < observed_peak_week[j])
    l <- sum(BP_idx)
    season_list[Diff_row_idx : (Diff_row_idx + l - 1)] <- rep(analysis_seasons[j], l)
    LS_timing[Diff_row_idx : (Diff_row_idx + l - 1)] <- target$peak_week_log_score[BP_idx] - log(1/33)
    LS_int[Diff_row_idx : (Diff_row_idx + l - 1)] <- target$peak_height_log_score[BP_idx] - log(1/27)
    model_list[Diff_row_idx : (Diff_row_idx + l - 1)] <- rep(model_name[i], l)
    Diff_row_idx <- Diff_row_idx + l   
  }
}
k <- length(LS_int)
PT_table <- data.frame(season = season_list,
                       Model = model_list,
                       LS_diff = c(LS_int, LS_timing),
                       group = c(rep("Peak Intensity",k), rep("Peak Timing", k)))

library(ggplot2)
library(gridExtra)
pdf(file = "./Results/plots/Peakbox.pdf", height = 6, width = 8)

p <- ggplot(data = PT_table, aes(x = Model, y = LS_diff)) + 
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot( outlier.shape = NA) +
  coord_cartesian(ylim = c(-4, 4)) +
  scale_x_discrete(limits = model_name) +
  stat_summary(fun.y = mean, geom="point", position=position_dodge(width=0.75), shape=18)  +  
  theme_bw() +
  theme(axis.text.x = element_text(angle=40))

x1 <- p + facet_grid(rows = vars( group), cols = vars(season), scales="free_x" , as.table=TRUE) + xlab("Model") +
  ylab("Log score difference") + 
  guides(fill=guide_legend(title="Season"))

q <- ggplot(data = data[data$train == FALSE & data$season_week %in% 10:42,], aes(x = season_week, y = weighted_ili)) +
  geom_line() +
  xlab("Season week") +
  ylab("wILI (%)")
x2 <- q + facet_grid(cols = vars(season), scales="free_x" , as.table=TRUE) + xlab("Season week") +
  ylab("wILI (%)") + 
  guides(fill=guide_legend(title="Season")) + theme_bw()

grid.arrange(x1, x2, nrow = 2, ncol = 1, heights = c(2, 1))
dev.off()
