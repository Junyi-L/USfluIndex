library(cdcfluview)
library(lubridate)
library(dplyr)
library(here)

MMWRweekday <- function(date) {
  factor(strftime(as.Date(date), "%w"), # Sunday is 0
         levels = 0:6,
         labels = c('Sunday','Monday','Tuesday','Wednesday','Thursday','Friday','Saturday'))
}

start_date = function(year) {
  # Finds start state for this calendar year
  jan1 = as.Date(paste(year, '-01-01', sep=''))
  wday = as.numeric(MMWRweekday(jan1))
  jan1 - (wday-1) + 7*(wday>4)
}

MMWRweek2Date <- function(MMWRyear, MMWRweek, MMWRday = NULL) {
  stopifnot(all(is.numeric(MMWRyear)))
  stopifnot(all(is.numeric(MMWRweek)))
  stopifnot(all(0 < MMWRweek & MMWRweek < 54))
  stopifnot(length(MMWRyear) == length(MMWRweek))
  if (is.null(MMWRday)) 
    MMWRday = rep(1, length(MMWRweek))
  stopifnot(all(0 < MMWRday & MMWRday < 8))
  jan1 = start_date(MMWRyear)
  return(jan1 + (MMWRweek - 1) * 7 + MMWRday - 1)
  
}

usflu <- ilinet(region = "national", years = 1997:2018)
# load data date: 11.03.2019
usflu <- as.data.frame(usflu)

data <- transmute(usflu,
                  region_type = region_type,
                  region = region,
                  year = year,
                  week = week,
                  weighted_ili = as.numeric(weighted_ili))
data[data$weighted_ili == 0,]$weighted_ili <- NA

data <- data.table(data)
data[, week_number := max(week), by = year]
data[, InPeriod := week/week_number, by = year]

data[, time := MMWRweek2Date(MMWRyear = year, MMWRweek = week)]
data$time_index <- as.integer(data$time -  ymd(paste("1970", "01", "01", sep = "-")))

data <- data[data$year <= 2018,]

data <- data.frame(data)


save(data,file = here("./Data/usflu.RData"))


