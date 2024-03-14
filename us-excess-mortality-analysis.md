# Using the excessmort package to generate excess mortality plot for USA
2017-2023
Rafael A. Irizarry

This file contains two parts. First the data wrangling part, which is
somewhat complex. Then we use compute expected death counts and plot
along with observed counts. This is relatively simple with the
`excessmort` package.

## Data wrangling

The first step is download data death count data from the CDC.

``` r
library(httr2)
library(jsonlite)
library(data.table)
library(lubridate)
library(excessmort)
mmwr_to_date <- function(mmwr_year, mmwr_week) {
  first_day <- floor_date(make_date(mmwr_year, 1, 4) , unit = "week")
  date <- first_day + weeks(mmwr_week - 1) + days(6) ## plus six to get week_end
  return(date)
}

api <- "https://data.cdc.gov/resource/3yf8-kanr.json"
dt2014_2019 <- request(api) |> 
  req_url_query("$limit" = 10000000) |>
  req_perform() |> resp_body_string() |> 
  fromJSON(flatten = TRUE) |>
  setDT()

dat1 <- dt2014_2019[jurisdiction_of_occurrence == "United States"]
dat1[, date := as_date(weekendingdate)]
dat1 <- dat1[, c("date", "allcause")]
setnames(dat1, "allcause", "outcome")

api <- "https://data.cdc.gov/resource/r8kw-7aab.json"
dt2020_present <- request(api) |> 
  req_url_query("$limit" = 10000000) |>
  req_perform() |> resp_body_string() |> 
  fromJSON(flatten = TRUE) |>
  setDT()

dat2 <- dt2020_present[state == "United States" & group == "By Week"]
dat2 <- dat2[, c("end_date", "total_deaths")]
setnames(dat2, c("end_date", "total_deaths"), c("date", "outcome"))
dat2[, date := as_date(date)]

dat <- rbindlist(list(dat1, dat2))
dat[, outcome := as.numeric(outcome)]
dat <- dat[date <= as_date("2023-12-31")]
dat <- dat[order(date),]
```

Next we get population estimates from the Census. You will need a census
API key which you can get from
<https://api.census.gov/data/key_signup.html>. Create a file called
`census-key.R` in your working directory with the following line in it:

``` r
census_key <- "YOURKEYHERE"
```

Next you can obtain population with this code:

``` r
source("census-key.R") 

api <- "https://api.census.gov/data/2019/pep/population"
pop1 <- request(api) |>  
  req_url_query(get = I("POP,DATE_CODE,DATE_DESC"), 
                `for` = I("us:*"),
                key = census_key) |>
  req_perform() |>
  resp_body_string() |> 
  fromJSON(flatten = TRUE) |>
  as.data.table()

pop1 <- pop1 |> janitor::row_to_names(1) 
pop1 <- pop1[!grepl("4/1", DATE_DESC)]
pop1 <- data.table(year = 2010 + as.numeric(pop1$DATE_CODE) - 3, population = as.numeric(pop1$POP))

pop2 <- fread("https://www2.census.gov/programs-surveys/popest/datasets/2020-2023/state/totals/NST-EST2023-ALLDATA.csv") 
years <- 2020:2023
pop2 <- pop2[NAME == "United States", ]
cols <- paste0("POPESTIMATE", years)
pop2 <- data.table(year = years, population = unlist(pop2[,..cols]))

pop <- rbindlist(list(pop1,pop2))[order(year)]
pop <- approx_demographics(pop, first_day = min(dat$date), last_day = max(dat$date))
```

Finally, merge the population and mortality counts data.

``` r
dat <- merge(dat, pop, by = "date", all.x = TRUE)
```

## Computing expected counts

``` r
library(tidyverse)
```

Because there was an upward trend from 2014 to 2016 that seemed to have
leveled off in 2017 we use only data from 2017 to 2020 to estimated
expected counts.

``` r
dat <- dat[date >= make_date(2017,1,1)]
```

We exclude two flu seasons that were worse than usual along with the
dates since the pandemic started:

``` r
exclude_dates <- c(seq(make_date(2017, 1, 1), make_date(2017, 2, 15), by = "day"),
                   seq(make_date(2018, 1, 1), make_date(2018, 2, 15), by = "day"),
                   seq(make_date(2020, 3, 1), make_date(2023, 12, 31), by = "day"))
```

Now we use the `excessmort` package:

``` r
expected <- compute_expected(dat, exclude = exclude_dates)
```

    No frequency provided, determined to be 52 measurements per year.

    Overall death rate is 9.38.

## Observed versus expected plot

To make the plot we use the `expected_diagnostic` function which we
recommend always using when performing this type of analysis. This is
particularly important when you are extrapolating.

``` r
diagnostics <- expected_diagnostic(expected)
diagnostics$expected + ggtitle("Total deaths in the USA from 2017 to 2023")
```

![](us-excess-mortality-analysis_files/figure-commonmark/unnamed-chunk-9-1.png)

## Excess mortality

We can also use this package to estimate excess mortality **assuming the
expected mortality model extrapolates correctly**. We compute it for
every year in the pandemic

``` r
em <- excess_model(expected, start = make_date(2015, 1, 1), end = make_date(2023, 12, 31), exclude = exclude_dates, 
                   knots.per.year = 18, 
                   intervals = list(seq(make_date(2020, 1, 1), make_date(2020, 12, 31), by = "day"),
                                    seq(make_date(2021, 1, 1), make_date(2021, 12, 31), by = "day"),
                                    seq(make_date(2022, 1, 1), make_date(2022, 12, 31), by = "day"),
                                    seq(make_date(2023, 1, 1), make_date(2023, 12, 31), by = "day")))
```

We can then make a table of excess deaths:

``` r
em$excess |> mutate(year = year(start)) |> mutate(lower = excess - 1.96*sd, upper = excess + 1.96*sd) |>
  select(year, excess, lower, upper) |>
  mutate(across(-year, function(x) paste0(round(x/1000,1),"K"))) |>
  arrange(desc(year)) |>
  knitr::kable()
```

| year | excess | lower  | upper  |
|-----:|:-------|:-------|:-------|
| 2023 | 145.9K | 135.7K | 156.2K |
| 2022 | 378.4K | 368.1K | 388.7K |
| 2021 | 576.7K | 566.6K | 586.9K |
| 2020 | 479.9K | 469.7K | 490K   |

and a plot:

``` r
tibble(date = em$date, obs = em$observed, excess = em$expected*(exp(em$fitted) - 1), 
       upper = em$expected*(exp(em$fitted + 1.96*em$log_expected_se) - 1),
       lower = em$expected*(exp(em$fitted - 1.96*em$log_expected_se) - 1)) |>
  filter(date >= make_date(2022,1,1)) |>
  ggplot(aes(date, excess)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "purple") +
  geom_line(color = "purple") +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw()
```

![](us-excess-mortality-analysis_files/figure-commonmark/unnamed-chunk-12-1.png)