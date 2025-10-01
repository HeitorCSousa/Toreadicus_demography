# Load packages -----------------------------------------------------------
# require(devtools)
# install.packages("INLA",repos=c(getOption("repos"),
#                                 INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

library(sf)
library(lubridate)
library(ncdf4)
library(tidyverse)
library(missForest)
library(gamm4)
library(ggplot2)
library(usdm)
library(corrplot)
library(openCR)
library(reshape)
library(jagsUI)
library(rjags)
library(scales)
library(ggmcmc)
library(INLA)
library(INLAtools)
library(dplyr)
library(mgcv)
library(mgcViz)
library(MODISTools)
library(raster)
library(terra)
library(rts)
library(spacetime)
library(viridis)
library(cowplot)
library(R2ucare)
library(rosm)
library(gimms)
library(tidyr)
library(locfit)
library(MuMIn)
library(AICcmodavg)
library(parallel)
library(RJDemetra)

# Check if other packages are instalee
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

package_vec <- c(
  "tidyr", # for turning rasters into ggplot-dataframes
  "ggplot2", # for plotting
  "viridis", # colour palettes
  "cowplot", # gridding multiple plots
  "rosm", # obtaining satellite maps
  "gimms" # to get some pre-existing data to match in our downscaling
)
sapply(package_vec, install.load.package)


# Read trap geographical coordinates --------------------------------------

traps.pts <- read.table("data/Points_Traps.txt", h = T)
traps.pts

(traps.pts.PEL <- subset(traps.pts, local == "PEL"))
coordinates(traps.pts.PEL) <- c("long", "lat")
proj4string(traps.pts.PEL) <- CRS("+proj=longlat +datum=WGS84")
traps.pts.PEL
traps.pts.PEL <- st_as_sf(traps.pts.PEL)

env.vars <- readRDS("data/macro_micro_clima_PEL_imp.rds")


# Critical temperatures ---------------------------------------------------

# Load data
ct.pel <- read.table("data/CT_PEL.txt", h = T)
ct.Toreadicus <- ct.pel[ct.pel$species == "T_oreadicus", ]
ct.Toreadicus

table(ct.Toreadicus$sex)


# Locomotor performance ---------------------------------------------------

# Load data
des.loc <- read.table("data/desemp_loc_PEL.txt", h = T)

head(des.loc)

table(des.loc$species, des.loc$PEL)
attach(des.loc)

summary(veloc)
summary(acel)
boxplot(veloc[veloc < 20])

# Aggregate data by maximum velocity achieved by each run
data <- aggregate(
  abs(veloc),
  by = list(
    species,
    PEL,
    temp,
    run,
    sex,
    SVL
  ),
  FUN = max,
  na.rm = T
)

head(data)
names(data) <- c(
  "species",
  "PEL",
  "temp",
  "run",
  "sex",
  "SVL",
  "Veloc"
)

head(data)

data.ctmin <- data.frame(
  "species" = ct.pel$species,
  "PEL" = ct.pel$pel,
  "temp" = ct.pel$Ctmin,
  "run" = c(rep(NA, length(ct.pel$Ctmin))),
  "sex" = ct.pel$sex,
  "SVL" = ct.pel$svl,
  "Veloc" = c(rep(0, length(ct.pel$Ctmin)))
)

data.ctmax <- data.frame(
  "species" = ct.pel$species,
  "PEL" = ct.pel$pel,
  "temp" = ct.pel$Ctmax,
  "run" = c(rep(NA, length(ct.pel$Ctmin))),
  "sex" = ct.pel$sex,
  "SVL" = ct.pel$svl,
  "Veloc" = c(rep(0, length(ct.pel$Ctmin)))
)

data.complete <- rbind(
  data,
  data.ctmin,
  data.ctmax
)

head(data.complete)
tail(data.complete)


# Locomotor performance curves --------------------------------------------

detach(des.loc)

# Tropidurus oreadicus
# Generalized additive mixed effect models
m.Toreadicus.sprint <- gamm4(
  Veloc ~ t2(temp) + SVL,
  random = ~ (1 | PEL),
  data = data.complete[data.complete$species == "T_oreadicus", ]
)
summary(m.Toreadicus.sprint$gam)

m.Toreadicus.sprint <- gamm4(
  Veloc ~ t2(temp),
  random = ~ (1 | PEL),
  data = data.complete[data.complete$species == "T_oreadicus", ]
)
summary(m.Toreadicus.sprint$gam)

plot(
  m.Toreadicus.sprint$gam,
  main = expression(italic("Tropidurus oreadicus")),
  ylab = "Locomotor performance",
  xlab = "Temperature (°C)"
)

# Predict locomotor performance
preddata_tpc_Toreadicus <- data.frame(temp = seq(10, 50, 0.1))

# Make prediction
pred_tpc_Toreadicus <- predict(
  m.Toreadicus.sprint$gam,
  preddata_tpc_Toreadicus,
  se.fit = T
)

## Attach mean predictions and errors to dataframe
preddata_tpc_Toreadicus$predicted <- pred_tpc_Toreadicus$fit
preddata_tpc_Toreadicus$se <- pred_tpc_Toreadicus$se.fit
preddata_tpc_Toreadicus$lower <- pred_tpc_Toreadicus$fit -
  pred_tpc_Toreadicus$se.fit
preddata_tpc_Toreadicus$upper <- pred_tpc_Toreadicus$fit +
  pred_tpc_Toreadicus$se.fit

ggplot(preddata_tpc_Toreadicus, aes(x = temp, y = predicted)) +
  geom_line() +
  geom_point(
    data = data.complete[data.complete$species == "T_oreadicus", ],
    aes(x = temp, y = Veloc),
    alpha = 0.5
  ) +
  geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 2, alpha = 0.1) +
  labs(
    title = "Tropidurus oreadicus",
    y = "Predicted Speed (cm/s)",
    x = "Temperature (°C)"
  ) +
  lims(x = c(10, 45), y = c(0, 2.25)) +
  theme(plot.title = element_text(hjust = 0.5))

# Optimal temperature
preddata_tpc_Toreadicus[
  preddata_tpc_Toreadicus$predicted == max(preddata_tpc_Toreadicus$predicted),
]
# 28.1°C

# Predict locomotor performance for microclimate data
env.vars$Toreadicus_perf <- predict.gam(
  m.Toreadicus.sprint$gam,
  newdata = data.frame(temp = env.vars$temp_med)
)


# Preferred temperatures --------------------------------------------------
# Load data
tpref_PEL <- read.table("data_Tpref.txt", h = T)
tpref_Toreadicus <- dplyr::filter(tpref_PEL, sp == "T_oreadicus")
table(tpref_Toreadicus$pel, tpref_Toreadicus$sp)

# Function to calculate hours of activity
hvtFUN <- function(temp.envr, temp.lab, quantiles, radiation) {
  vtmin <- quantile(temp.lab, quantiles[1], na.rm = TRUE)
  vtmax <- quantile(temp.lab, quantiles[2], na.rm = TRUE)
  hv <- ifelse(temp.envr > vtmin & temp.envr < vtmax, 1, 0)
  hv[radiation == 0] <- 0
  hv
}


# Hours of activity considering 90% percentile
env.vars$Toreadicus_ha90 <- hvtFUN(
  env.vars$temp_med,
  tpref_Toreadicus$temp,
  c(0.05, 0.95),
  env.vars$solrad.y
)

# Hours of activity considering 50% percentile
env.vars$Toreadicus_ha50 <- hvtFUN(
  env.vars$temp_med,
  tpref_Toreadicus$temp,
  c(0.25, 0.75),
  env.vars$solrad.y
)

# Aggregate data by trap and month
env.vars.month <-
  env.vars %>%
  dplyr::group_by(trap, year.GMT_3, month.GMT_3) %>%
  dplyr::summarise(
    tmed = mean(temp_med),
    tmin = min(temp_med),
    tmax = max(temp_med),
    tmed2m = mean(temp_2m.y),
    tmin2m = min(temp_2m.y),
    tmax2m = max(temp_2m.y),
    tmedskin = mean(temp_skin.y),
    tminskin = min(temp_skin.y),
    tmaxskin = max(temp_skin.y),
    rh = mean(rh_med),
    rhmin = min(rh_med),
    rhmax = max(rh_med),
    sol = mean(solrad.y),
    solmax = max(solrad.y),
    precip = sum(precip.y),
    Toreadicus_perf = mean(Toreadicus_perf),
    Toreadicus_ha50 = sum(Toreadicus_ha50),
    Toreadicus_ha90 = sum(Toreadicus_ha90)
  )

summary(env.vars.month)

# Find collinear variables
vifstep(as.data.frame(env.vars.month[, 4:21]))
vifcor(as.data.frame(env.vars.month[, 4:21]), th = .8)

# Correlation matrix plot
m <- cor(env.vars.month[, c(4:21)])
corrplot.mixed(m, order = "AOE", upper = "ellipse")

env.vars.month <- subset(
  env.vars.month,
  select = -c(tmaxskin, tmed2m, rh, tmin, tmax2m, rhmin, sol, tminskin, tmax)
)

# Exclude variables with collinearity problems: tmaxskin tmed2m rh tmin tmax2m rhmin sol tminskin tmax

# Demographic analyses ----------------------------------------------------

# Clean data
Toreadicus.data <- readRDS("data/Toreadicus_PEL_imp.rds")
head(Toreadicus.data)
tail(Toreadicus.data)

# Remove dead animals
Toreadicus.data$dead[is.na(Toreadicus.data$dead)] <- "N"
data.demography <- Toreadicus.data[Toreadicus.data$dead == "N", ]
head(data.demography)
table(data.demography$campaign)
table(data.demography$id)

# Remove NAs
complete <- complete.cases(data.demography[, c("fieldtag", "trap")])
Toreadicus.df <- droplevels(data.demography[complete, ])
head(Toreadicus.df)
str(Toreadicus.df)
table(Toreadicus.df$campaign)
table(Toreadicus.df$id)

## openCR ------------------------------------------------------------------

# Load trap locations
pts.traps <- read.table("Points_Traps.txt", h = T)
pts.traps.PEL <- pts.traps[pts.traps$local == "PEL", ]

(traps.pts.PEL <- subset(traps.pts, local == "PEL"))
coordinates(traps.pts.PEL) <- c("long", "lat")
proj4string(traps.pts.PEL) <- CRS("+proj=longlat +datum=WGS84")
traps.pts.PEL.utm <- spTransform(traps.pts.PEL, CRS = CRS("+init=epsg:32722"))
pts.traps.PEL$X <- coordinates(traps.pts.PEL.utm)[, 1]
pts.traps.PEL$Y <- coordinates(traps.pts.PEL.utm)[, 2]


Toreadicus.df.XY <- merge(Toreadicus.df, pts.traps.PEL, all.x = T)
head(Toreadicus.df.XY)

capts <- data.frame(
  Session = Toreadicus.df.XY$local,
  ID = Toreadicus.df.XY$fieldtag,
  occasion = Toreadicus.df.XY$campaign,
  X = Toreadicus.df.XY$X,
  Y = Toreadicus.df.XY$Y,
  sex = as.factor(Toreadicus.df.XY$sex),
  year = as.factor(Toreadicus.df.XY$year)
)

summary(capts)

# Immature individuals (sex unknown)
capts$sex[capts$sex == "I"] <- NA

# remove NAs
complete <- complete.cases(capts[, c("ID", "X", "Y", "sex")])
capts <- droplevels(capts[complete, ])
summary(capts)

traps.xy <- data.frame(
  trapID = pts.traps.PEL$local,
  x = pts.traps.PEL$X,
  y = pts.traps.PEL$Y
)
traps.xy <- read.traps(data = traps.xy[, -1], detector = "multi")

Toreadicus.CHxy <- make.capthist(capts, traps.xy, fmt = "XY", noccasions = 46)

write.capthist(Toreadicus.CHxy)

mesh <- make.mask(traps.xy)

summary(Toreadicus.CHxy)
summary(traps(Toreadicus.CHxy))

par(mfrow = c(1, 1))
plot(mesh)
plot(Toreadicus.CHxy, tracks = T, add = T)

# CJS
cjs.0 <- openCR.fit(Toreadicus.CHxy)
summary(cjs.0)

cjs.t <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + year, p ~ -1 + year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  trace = 2,
  start = cjs.0
)

summary(cjs.t)


cjs.full <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex * year, p ~ -1 + sex * year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.full)

cjs.fulladd <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex + year, p ~ -1 + sex + year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.fulladd)

cjs.psextphisex <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex, p ~ -1 + sex + year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.psextphisex)

cjs.psexphisext <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex + year, p ~ -1 + sex),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.psexphisext)

cjs.ptphisext <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex + year, p ~ -1 + year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.ptphisext)

cjs.psextphit <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + year, p ~ -1 + sex + year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.psextphit)

cjs.psexphit <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + year, p ~ -1 + sex),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.psexphit)

cjs.ptphisex <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex, p ~ -1 + year),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.t
)
summary(cjs.ptphisex)


cjs.sex <- openCR.fit(
  Toreadicus.CHxy,
  model = list(phi ~ -1 + sex, p ~ -1 + sex),
  method = "BFGS",
  details = list(
    control = list(
      maxit = 100000,
      reltol = 1e-5
    )
  ),
  ncores = 6,
  start = cjs.0
)
summary(cjs.sex)


AIC(
  cjs.0,
  cjs.psexphisext,
  cjs.psexphit,
  cjs.psextphisex,
  cjs.psextphit,
  cjs.ptphisex,
  cjs.ptphisext,
  cjs.sex,
  cjs.t,
  cjs.full,
  cjs.fulladd
)

summary(cjs.psextphit)
summary(cjs.fulladd)

pred.fullcjs.df <- modelAverage(
  openCRlist(
    cjs.0,
    cjs.psexphisext,
    cjs.psexphit,
    cjs.psextphisex,
    cjs.psextphit,
    cjs.ptphisex,
    cjs.ptphisext,
    cjs.sex,
    cjs.t,
    cjs.full,
    cjs.fulladd
  ),
  newdata = expand.grid(
    sex = as.factor(c("F", "M")),
    year = as.factor(c(2018:2021))
  )
)

pred.fullcjs.df <- list(
  p = pred.fullcjs.df[,, "p"],
  phi = pred.fullcjs.df[,, "phi"]
)

pred.fullcjs.df$p <- cbind(
  expand.grid(
    sex = as.factor(c("F", "M")),
    year = as.factor(c(2018:2021))
  ),
  pred.fullcjs.df$p
)

pred.fullcjs.df$phi <- cbind(
  expand.grid(
    sex = as.factor(c("F", "M")),
    year = as.factor(c(2018:2021))
  ),
  pred.fullcjs.df$phi
)


ggplot(pred.fullcjs.df$p, aes(x = year:sex, y = estimate, color = sex)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl)) +
  labs(x = "", y = "Capture") +
  theme_classic()

ggplot(pred.fullcjs.df$phi, aes(x = year:sex, y = estimate, color = sex)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl)) +
  labs(x = "", y = "Survival") +
  theme_classic()

# JAGS Pradel model ----------------------------------------------

# Use monthly data

# Filter variables of interest
table1 <- Toreadicus.df.XY[
  Toreadicus.df.XY$recapture != "(s)",
  c(
    "fieldtag",
    "campaign",
    "sex",
    "svl",
    "mass",
    "recapture",
    "trap"
  )
]
str(table1)


# Identify captures without ID
table2 <- complete.cases(table1[, c("fieldtag")])

# Remove captures without ID and SVL
table3 <- table1[table2, ]

# Order the data by ID and campaign (month)
table4 <- table3[order(table3$fieldtag, table3$campaign), ]
str(table4)
summary(table4)

# Convert ID from factor to character
table5 <- droplevels(table4)
table5$fieldtag <- as.character(table5$fieldtag)
table5$trap <- as.character(table5$trap)
str(table5)


# Calculate recapture frequencies
recap.table <- data.frame(table(table5$fieldtag))
names(recap.table) <- c("id", "captures")
recap.table
table(recap.table$captures)
218 +
  (2 * 69) +
  (3 * 37) +
  (4 * 17) +
  (5 * 13) +
  (6 * 7) +
  (7 * 4) +
  (8 * 3) +
  11 +
  12 +
  13
# two captures=1, 3 captures=2, 4 captures=3 ...

# Subset the dataframe to use in the analysis
head(table5)

Age <- c(rep(NA, nrow(table5)))
Age

datA <- data.frame(
  table5$campaign,
  table5$sex,
  table5$fieldtag,
  table5$svl,
  table5$mass,
  Age,
  table5$trap
)
names(datA) <- c("Year", "Sex", "TrueID", "SVL", "Massa", "Age", "Trap")
datA$Year <- datA$Year
datA$Age[datA$SVL <= 35] <- 0
head(datA)
tail(datA)
str(datA)

###  del is the time period since the first capture (0 in the first capture)

del <- c() ### months since the first capture

for (i in 1:nrow(datA)) {
  del[i] <- datA$Year[i] - min(datA$Year[datA$TrueID == datA$TrueID[i]])
}

table(del)


trap <- as.integer(datA$Trap)
ind <- as.numeric(factor(datA$TrueID))
(time <- datA$Year)
y <- datA$SVL
m <- nrow(datA) ### number of observations

age <- c() ## age at first capture
for (a in 1:m) {
  age[a] <- datA$Age[ind == a][1]
}

year <- c()
for (a in 1:m) {
  year[a] <- datA$Year[ind == a][1]
}

head(datA)
tail(datA)
str(datA)

known.states.cjs <- function(ch) {
  state <- ch
  for (i in 1:dim(ch)[1]) {
    n1 <- min(which(ch[i, ] == 1))
    n2 <- max(which(ch[i, ] == 1))
    state[i, n1:n2] <- 1
    state[i, n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

cjs.init.z <- function(ch, f) {
  for (i in 1:dim(ch)[1]) {
    if (sum(ch[i, ]) == 1) {
      next
    }
    n2 <- max(which(ch[i, ] == 1))
    ch[i, f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]) {
    ch[i, 1:f[i]] <- NA
  }
  return(ch)
}

eh <- cast(
  datA,
  TrueID ~ Year,
  fun.aggregate = function(x) as.numeric(length(x) > 0),
  value = "SVL"
)
eh <- eh[, 2:ncol(eh)]
eh.all <- seq(min(datA$Year), max(datA$Year)) # fill all the years ignored
missing <- eh.all[!(eh.all %in% names(eh))]
col <- matrix(0, nrow = nrow(eh), ncol = length(missing))
colnames(col) <- missing + 100
eh <- cbind(eh, col)
eh <- eh[, sort.int(as.integer(colnames(eh)))]
head(eh)

table(rowSums(eh)) # Number of captures per individual
sum(eh) # Number of captures
nrow(eh) # Number of individuals

# Mark-recapture statistics
table(rowSums(eh)) # Number of captures for each individual
sum(eh) # Number of captures

(sum(eh) - nrow(eh)) / nrow(eh) # Mean recapture rate
mean(rowSums(eh)) # Mean capture rate
round(mean(rowSums(eh) - 1), 2) # Mean recapture rate
round(sd(rowSums(eh) - 1), 2) # SD capture rate
nc <- rowSums(eh) - 1
nc[nc == 0] <- NA
round(sd(nc, na.rm = T), 2) # SD recapture rate

table(rowSums(known.states.cjs(eh), na.rm = T)) # time period between first and last captures
table(del)
mean(table(rowSums(known.states.cjs(eh), na.rm = T))[-1]) # Mean time period between first and last captures
sd(table(rowSums(known.states.cjs(eh), na.rm = T))[-1]) # SD time period between first and last captures

(f <- apply(eh, 1, function(x) which(x == 1)[1]))
(nind <- nrow(eh))
(n.occasions <- ncol(eh))
m
nind

datA$Sex[datA$Sex == "I"] <- NA
table(datA$TrueID, datA$Sex)
colSums(table(datA$TrueID, datA$Sex))


sex <- cast(datA, TrueID ~ ., value = "Sex", fun.aggregate = function(x) {
  tail(x, 1)
})
table(sex$`(all)`)
round(prop.table(table(sex$`(all)`)) * 100, 1)
table(is.na(sex$`(all)`))
round(prop.table(table(is.na(sex$`(all)`))) * 100, 1)
sex.m <- eh * as.numeric(as.factor(sex$`(all)`))

sex.months <- t(as.matrix(as.data.frame(apply(sex.m, 2, tabulate))))
sex.months <- cbind(
  sex.months,
  c(rep(2018, 11), rep(2019, 12), rep(2020, 12), rep(2021, 11))
)
females <- tapply(sex.months[, 1], sex.months[, 3], sum)
males <- tapply(sex.months[, 2], sex.months[, 3], sum)
(sex.year <- rbind(females, males))
chisq.test(sex.year, simulate.p.value = T)

pred.fullcjs.df$p$estimate[pred.fullcjs.df$p$sex == "F"]

(sex.year.cjs <- rbind(
  round(females / pred.fullcjs.df$p$estimate[pred.fullcjs.df$p$sex == "F"]),
  round(males / pred.fullcjs.df$p$estimate[pred.fullcjs.df$p$sex == "M"])
))
chisq.test(sex.year.cjs, simulate.p.value = T)
chisq.test(sex.year.cjs[, 1], simulate.p.value = T)
chisq.test(sex.year.cjs[, 2], simulate.p.value = T)
chisq.test(sex.year.cjs[, 3], simulate.p.value = T)
chisq.test(sex.year.cjs[, 4], simulate.p.value = T)

# Derive data for the model
# e = index of the oldest observation
get.first <- function(x) min(which(x != 0))
e <- apply(eh, 1, get.first)

# l = index of the last observation
get.last <- function(x) max(which(x != 0))
l <- apply(eh, 1, get.last)

table(l - e) # Difference between first and last capture
boxplot(l - e)

# u = number of animals observed for the first time at i
u <- data.frame(table(e))


u.all <- seq(1, n.occasions) # Fill all the months ignored
missing <- u.all[!(u.all %in% u$e)]
df <- data.frame(e = as.factor(missing), Freq = rep(0, length(missing)))
names(u) <- names(df)
u <- rbind(u, df)
u$e <- as.numeric(levels(u$e))[u$e]
u <- u[order(u$e), ]
u <- u$Freq
u

# n = number of animals observed at i
n <- colSums(eh)

# v = number of animals observed for the last time at i
v <- data.frame(table(l))

v.all <- seq(1, n.occasions) # Fill all the months ignored
missing <- v.all[!(v.all %in% v$l)]
df <- data.frame(l = as.factor(missing), Freq = rep(0, length(missing)))
names(v) <- names(df)
v <- rbind(v, df)
v$l <- as.numeric(levels(v$l))[v$l]
v <- v[order(v$l), ]
v <- v$Freq
v

# d = number of animals removed from the population at i (can add, but not in our case)
d <- rep(0, dim(eh)[2])

# Environmental variables

env.vars.month$t <- rep(1:46, 25)

tmed <- cast(env.vars.month, trap ~ t, value = "tmed")
tmed <- as.matrix(tmed[, -1])
tmin2m <- cast(env.vars.month, trap ~ t, value = "tmin2m")
tmin2m <- as.matrix(tmin2m[, -1])
tmedskin <- cast(env.vars.month, trap ~ t, value = "tmedskin")
tmedskin <- as.matrix(tmedskin[, -1])
rhmax <- cast(env.vars.month, trap ~ t, value = "rhmax")
rhmax <- as.matrix(rhmax[, -1])
solmax <- cast(env.vars.month, trap ~ t, value = "solmax")
solmax <- as.matrix(solmax[, -1])
precip <- cast(env.vars.month, trap ~ t, value = "precip")
precip <- as.matrix(precip[, -1])
Toreadicus_perf <- cast(env.vars.month, trap ~ t, value = "Toreadicus_perf")
Toreadicus_perf <- as.matrix(Toreadicus_perf[, -1])
Toreadicus_ha90 <- cast(env.vars.month, trap ~ t, value = "Toreadicus_ha90")
Toreadicus_ha90 <- as.matrix(Toreadicus_ha90[, -1])

time.env.cov <- data.frame(
  tmedJS = colMeans(tmed),
  tmin2mJS = colMeans(tmin2m),
  tmedskinJS = colMeans(tmedskin),
  rhmaxJS = colMeans(rhmax),
  solmaxJS = colMeans(solmax),
  precipJS = colMeans(precip),
  Toreadicus_perfJS = colMeans(Toreadicus_perf),
  Toreadicus_ha90JS = colMeans(Toreadicus_ha90)
)

as.matrix(time.env.cov)
# Package data
bugs.data <- list(
  u = u,
  n = n,
  v = v,
  d = d,
  first = f,
  ind = ind,
  nind = dim(eh)[1],
  n.occasions = dim(eh)[2],
  y = eh,
  z = known.states.cjs(eh),
  time = time,
  env.JS = as.matrix(scale(time.env.cov)),
  trap = trap
)

saveRDS(bugs.data, "Toreadicus_data.rds")

# Initial values
inits <- function() {
  list(
    tauphib = 1,
    betaTphi = rep(0, 4),
    varphi = rep(0, 4),
    tauphibJS = 1,
    betaTphiJS = rep(0, 8),
    varphiJS = rep(0, 8),
    sigma.phiJS = runif(1, 1, 2),
    sigma.f = runif(1, 0.75, 1.5),
    taufb = 1,
    betaTf = rep(0, 8),
    varf = rep(0, 8),
    taupb = 1,
    betaTp = rep(0, 8),
    varp = rep(0, 8),
    taupbJS = 1,
    betaTpJS = rep(0, 8),
    varpJS = rep(0, 8),
    sigma.pJS = runif(1, 0.1, 0.5)
  )
}

# Define the parameters to be monitored
parameters <- c(
  "phiJS",
  "alpha.phiJS",
  "sigma.phiJS",
  "betaphiJS",
  "varphiJS",
  "f",
  "alpha.f",
  "sigma.f",
  "betaf",
  "varf",
  "pJS",
  "alpha.pJS",
  "sigma.pJS",
  "betapJS",
  "varpJS",
  "rho"
)


# Specify model in BUGS language
sink("cjs-Pradel-Toreadicus.jags")
cat(
  "

data{
C<-10000
zeros<-0
}

model {

   #################
   #Pradel JS model#
   #################



###########PRIORn.occasions#######################
gamma[1]<-0
phiJS[n.occasions]<-0

#logit constraint for survival probability(phiJS)
for(t in 1:(n.occasions-1)){
phiJS[t]<-1/(1+exp(-logit.phiJS[t]))
logit.phiJS[t]<- alpha.phiJS + eps.phiJS[t] + inprod(env.JS[t,],betaphiJS)
eps.phiJS[t]~dnorm(0,tau.phiJS)
}

for(j in 1:8){
varphiJS[j]~dbern(0.5)
betaTphiJS[j]~dnorm(0,tauphibJS)
betaphiJS[j]<-varphiJS[j]*betaTphiJS[j]
}

alpha.phiJS~dnorm(0,0.01)

#alpha.phiJS on prob scale
mean.phiJS<-1/(1+exp(-alpha.phiJS))

#environmental parameters
tauphibJS~dgamma(1,0.001)


#log constraint for recruitment rate(f)

for(t in 1:(n.occasions-1)){
f[t]<-exp(log.f[t])
log.f[t]<- alpha.f + eps.f[t]+ inprod(env.JS[t,],betaf)
eps.f[t]~dnorm(0,tau.f)
}

for(j in 1:8){
varf[j]~dbern(0.5)
betaTf[j]~dnorm(0,taufb)
betaf[j]<-varf[j]*betaTf[j]
}

alpha.f~dnorm(0,0.01)


#alpha.f on prob scale
mean.f<-exp(alpha.f)

#environmental parameters
taufb~dgamma(1,0.001)


for(t in 1:n.occasions){
#logit constraint for detectability (pJS)
pJS[t]<-1/(1+exp(-logit.pJS[t]))
logit.pJS[t]<- alpha.pJS + eps.p[t]+ inprod(env.JS[t,],betapJS)
eps.p[t]~dnorm(0,tau.pJS)

#no constraints formula
mu[t]~dunif(0,1)
}

for(j in 1:8){
varpJS[j]~dbern(0.5)
betaTpJS[j]~dnorm(0,taupbJS)
betapJS[j]<-varpJS[j]*betaTpJS[j]
}

alpha.pJS~dnorm(0,0.01)

#alpha.pJS on prob scale
mean.pJS<-1/(1+exp(-alpha.pJS))

#environmental parameters
taupbJS~dgamma(1,0.001)

#temporal random variation
tau.f<-1/(sigma.f*sigma.f)
sigma.f~dunif(0,2)

tau.phiJS<-1/(sigma.phiJS*sigma.phiJS)
sigma.phiJS~dunif(0,2)

tau.pJS<-1/(sigma.pJS*sigma.pJS)
sigma.pJS~dunif(0,2)

###########LIKELIHOOD(ZERO-TRICK)######
zeros~dpois(zero.mean)
zero.mean<--L+C
L<-sum(l.num[1:n.occasions])-l.denom

#####log-likelihood for the first occasion
l.num[1]<-(u[1]*log(xi[1]))+(n[1]*log(pJS[1]))+(secondexpo[1]*log(1-pJS[1]))+
(thirdexpo[1]*log(phiJS[1]))+(fourthexpo[1]*log(mu[1]))+
(d[1]*log(1-mu[1]))+(fifthexpo[1]*log(1-(pJS[1]*(1-mu[1]))))+
(sixthexpo[1]*log(chi[1]))
xi[1]<-1
secondexpo_a[1]<-sum(u[1:1])
secondexpo_b[1]<-0
secondexpo[1]<-secondexpo_a[1]-secondexpo_b[1]-n[1]
thirdexpo[1]<-sum(v[2:n.occasions])
fourthexpo[1]<-n[1]-d[1]
fifthexpo[1]<-sum(u[2:n.occasions])
sixthexpo[1]<-v[1]-d[1]

#####log-likelihood for the last occasion
l.num[n.occasions]<-(u[n.occasions]*log(xi[n.occasions]))+(firstexpo[n.occasions]*(log(phiJS[n.occasions-1])-log(phiJS[n.occasions-1]+f[n.occasions-1])))+
(n[n.occasions]*log(pJS[n.occasions]))+(secondexpo[n.occasions]*log(1-pJS[n.occasions]))+
(fourthexpo[n.occasions]*log(mu[n.occasions]))+(d[n.occasions]*log(1-mu[n.occasions]))+
(fifthexpo[n.occasions]*log(1-(pJS[n.occasions]*(1-mu[n.occasions]))))+
(sixthexpo[n.occasions]*log(chi[n.occasions]))
chi[n.occasions]<-1

firstexpo[n.occasions]<-sum(u[1:(n.occasions-1)])
secondexpo_a[n.occasions]<-sum(u[1:n.occasions])
secondexpo_b[n.occasions]<-sum(v[1:(n.occasions-1)])
secondexpo[n.occasions]<-secondexpo_a[n.occasions]-secondexpo_b[n.occasions]-n[n.occasions]
fourthexpo[n.occasions]<-n[n.occasions]-d[n.occasions]
fifthexpo[n.occasions]<-0
sixthexpo[n.occasions]<-v[n.occasions]-d[n.occasions]

#####likelihood from occasion 2 to n.occasions-1
for(i in 2:(n.occasions-1)){
l.num[i]<-(u[i]*log(xi[i]))+(firstexpo[i]*(log(phiJS[i-1])-log(phiJS[i-1]+f[i-1])))+
(n[i]*log(pJS[i]))+(secondexpo[i]*log(1-pJS[i]))+
(thirdexpo[i]*log(phiJS[i]))+(fourthexpo[i]*log(mu[i]))+
(d[i]*log(1-mu[i]))+(fifthexpo[i]*log(1-(pJS[i]*(1-mu[i]))))+
(sixthexpo[i]*log(chi[i]))

#first exponent
firstexpo[i]<-sum(u[1:(i-1)])

#second exponent
secondexpo_a[i]<-sum(u[1:i])
secondexpo_b[i]<-sum(v[1:(i-1)])
secondexpo[i]<-secondexpo_a[i]-secondexpo_b[i]-n[i]

#third exponent
thirdexpo[i]<-sum(v[(i+1):n.occasions])

#fourth exponent
fourthexpo[i]<-n[i]-d[i]

#fifth exponent
fifthexpo[i]<-sum(u[(i+1):n.occasions])

#sixth exponent
sixthexpo[i]<-v[i]-d[i]
}

#####likelihood denominator
#1st product
PROD1[1]<-1
for(j in 1:(n.occasions-1)){
PROD1_tmp[1,j]<-0
}

#fill part of PROD1_tmp
for(i in 2:(n.occasions-1)){
for(j in i:(n.occasions-1)){
PROD1_tmp[i,j]<-0
}
}

for(i in 2:n.occasions){
for(j in 1:(i-1)){
PROD1_tmp[i,j]<-phiJS[j]*(1-(pJS[j]*(1-mu[j])))
}
}

PROD1[2]<-PROD1_tmp[2,1]
for(i in 3:n.occasions){
PROD1[i]<-prod(PROD1_tmp[i,1:(i-1)])
}

#2nd product
PROD2[n.occasions]<-1
for(i in 1:(n.occasions-1)){
for(j in (i+1):n.occasions){
PROD2_tmp[i,j]<-gamma[j]
}
}
#fill part of PROD2_tmp
for(i in 1:(n.occasions-1)){
for(j in 1:i){
PROD2_tmp[i,j]<-0
}
}

PROD2[n.occasions-1]<-PROD2_tmp[(n.occasions-1),n.occasions]
for(i in 1:(n.occasions-2)){
PROD2[i]<-prod(PROD2_tmp[i,(i+1):n.occasions])
}
for(i in 1:n.occasions){
denom_base_tmp[i]<-xi[i]*PROD1[i]*PROD2[i]*pJS[i]
}

denom_base <- sum(denom_base_tmp[])
denom_expo <- sum(u[1:n.occasions])
l.denom <- denom_expo * log(denom_base)

#################Define xi and chi
for(i in 2:n.occasions){
xi.tmp[i]<-(1-gamma[i])+
(gamma[i]*((1-pJS[i-1])/(1-(pJS[i-1]*(1-mu[i-1]))))*xi[i-1])
xi[i]<-max(xi.tmp[i],0.00001)
}
for(i in 1:(n.occasions-1)){
chi[i]<-(1-phiJS[i])+(phiJS[i]*(1-pJS[i+1])*chi[i+1])
}
#################Gamma and rho as derived parameter
for(i in 2:n.occasions){
rho[i]<-phiJS[i-1]+f[i-1]
}

for(i in 2:n.occasions){
gamma[i]<-phiJS[i-1]/(phiJS[i-1]+f[i-1])
}


}
",
  fill = TRUE
)
sink()

# MCMC settings
ni <- 150000
nt <- 1
nb <- 50000
nc <- 3

# Call JAGS from R (This can take a long time)
# pradel.Toreadicus<- jags(bugs.data, inits, parameters, "cjs-Pradel-Toreadicus.jags",
#                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
#                        parallel = T)#, codaOnly = parameters
#
# print(pradel.Toreadicus, digits = 3)
# saveRDS(pradel.Toreadicus, "results_pradel_Toreadicus.rds")
# summary(pradel.Toreadicus)
# pradel.Toreadicus.df <- pradel.Toreadicus$summary
# View(pradel.Toreadicus.df)
# write.table(pradel.Toreadicus.df,"results_pradel_Toreadicus.txt",sep="\t")

# Please, load the results below
pradel.Toreadicus <- readRDS("output/results_pradel_Toreadicus.rds")
pradel.Toreadicus.df <- pradel.Toreadicus$summary
View(pradel.Toreadicus.df)

PJS.df <- as.data.frame(pradel.Toreadicus.df[
  -c(47:64, 110:127, 174:191, 237),
])
head(PJS.df)
tail(PJS.df)
PJS.df$time <- c(1:46, 1:45, 1:46, 2:46)
PJS.df$date <- c(
  seq(as.Date("2018-02-01"), as.Date("2021-11-30"), "month"),
  seq(as.Date("2018-02-01"), as.Date("2021-10-30"), "month"),
  seq(as.Date("2018-02-01"), as.Date("2021-11-30"), "month"),
  seq(as.Date("2018-03-01"), as.Date("2021-11-30"), "month")
)
PJS.df$date <- as.character(PJS.df$date)
PJS.df$date <- as.Date(PJS.df$date)

PJS.df$month <- month(PJS.df$date, label = T, abbr = T)
PJS.df$year <- year(PJS.df$date)
PJS.df$param <- c(
  rep("phi", 46),
  rep("f", 45),
  rep("p", 46),
  rep("rho", 45)
)

names(PJS.df) <- c(
  "mean",
  "sd",
  "low2.5",
  "low.25",
  "median",
  "hi75",
  "hi97.5",
  "Rhat",
  "n.eff",
  "overlap0",
  "f",
  "time",
  "date",
  "month",
  "year",
  "param"
)

PJS.df <- na.omit(PJS.df)

rect1 <- data.frame(xmin = 1, xmax = 3, ymin = -Inf, ymax = Inf) # set rainy season
rect2 <- data.frame(xmin = 9, xmax = 15, ymin = -Inf, ymax = Inf)
rect3 <- data.frame(xmin = 21, xmax = 27, ymin = -Inf, ymax = Inf)
rect4 <- data.frame(xmin = 32, xmax = 39, ymin = -Inf, ymax = Inf)
rect5 <- data.frame(xmin = 45, xmax = 46, ymin = -Inf, ymax = Inf)
meses <- month(
  seq(as.POSIXlt("2018-02-01"), as.POSIXlt("2021-11-30"), "month"),
  label = T,
  abbr = T
)

# Plot results

# Survival (phi)
ggplot(PJS.df[PJS.df$param == "phi", ], aes(time)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = mean), size = 1.5) +
  geom_ribbon(aes(ymin = low2.5, ymax = hi97.5), alpha = .4) +
  geom_ribbon(aes(ymin = low.25, ymax = hi75), alpha = .8) +
  scale_x_continuous(
    breaks = c(1:45),
    labels = meses[-46],
    expand = c(0, 0.6)
  ) +
  labs(x = NULL, y = "Survival") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10)
  )

# Recruitment (f)
ggplot(PJS.df[PJS.df$param == "f", ], aes(time, mean)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = mean), size = 1.5) +
  scale_x_continuous(
    breaks = c(1:45),
    labels = meses[-46],
    expand = c(0, 0.6)
  ) +
  geom_ribbon(aes(ymin = low2.5, ymax = hi97.5), alpha = .4) +
  geom_ribbon(aes(ymin = low.25, ymax = hi75), alpha = .8) +
  labs(x = NULL, y = "Recruitment") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Capture probability (p)
ggplot(PJS.df[PJS.df$param == "p", ], aes(time, mean)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = mean), size = 1.5) +
  scale_x_continuous(breaks = c(1:46), labels = meses, expand = c(0, 0.6)) +
  geom_ribbon(aes(ymin = low2.5, ymax = hi97.5), alpha = .4) +
  geom_ribbon(aes(ymin = low.25, ymax = hi75), alpha = .8) +
  labs(x = NULL, y = "Capture") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Realized population growth (rho)
ggplot(PJS.df[PJS.df$param == "rho", ], aes(time, log(mean))) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = log(mean)), size = 1.5) +
  scale_x_continuous(
    breaks = c(1:45),
    labels = meses[-46],
    expand = c(0, 0.6)
  ) +
  geom_ribbon(aes(ymin = log(low2.5), ymax = log(hi97.5)), alpha = .4) +
  geom_ribbon(aes(ymin = log(low.25), ymax = log(hi75)), alpha = .8) +
  labs(x = NULL, y = "Population growth") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Locomotor performance mediated by temperature
ggplot(time.env.cov, aes(x = 1:46, y = Toreadicus_perfJS)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = Toreadicus_perfJS), size = 1.5) +
  scale_x_continuous(breaks = c(1:46), labels = meses, expand = c(0, 0.6)) +
  labs(x = NULL, y = "Locomotor performance") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

#
# Mean of population growth
rho.df <- PJS.df[PJS.df$param == "rho", ]
tapply(log(rho.df$mean), rho.df$year, mean)
tapply(log(rho.df$low2.5), rho.df$year, mean)
tapply(log(rho.df$hi97.5), rho.df$year, mean)
tapply(log(rho.df$low.25), rho.df$year, mean)
tapply(log(rho.df$hi75), rho.df$year, mean)

round(mean(log(rho.df$mean)), 3)
round(mean(log(rho.df$low2.5)), 3)
round(mean(log(rho.df$hi97.5)), 3)
round(mean(log(rho.df$low.25)), 3)
round(mean(log(rho.df$hi75)), 3)


tapply(rho.df$mean, rho.df$year, function(x) {
  exp(mean(log(x)))
})
exp(mean(log(rho.df$mean)))
exp(mean(log(rho.df$low2.5)))
exp(mean(log(rho.df$hi97.5)))
exp(mean(log(rho.df$low.25)))
exp(mean(log(rho.df$hi75)))

# Mean of survival
pradel.Toreadicus.df["alpha.phiJS", "mean"]
round(
  exp(pradel.Toreadicus.df["alpha.phiJS", "mean"]) /
    (1 + exp(pradel.Toreadicus.df["alpha.phiJS", "mean"])),
  3
) # Mean
round(
  exp(pradel.Toreadicus.df["alpha.phiJS", "2.5%"]) /
    (1 + exp(pradel.Toreadicus.df["alpha.phiJS", "2.5%"])),
  3
) # Low CI
round(
  exp(pradel.Toreadicus.df["alpha.phiJS", "97.5%"]) /
    (1 + exp(pradel.Toreadicus.df["alpha.phiJS", "97.5%"])),
  3
) # High CI

phi.df <- PJS.df[PJS.df$param == "phi", ]

# Mean of capture
pradel.Toreadicus.df["alpha.pJS", "mean"]
round(
  exp(pradel.Toreadicus.df["alpha.pJS", "mean"]) /
    (1 + exp(pradel.Toreadicus.df["alpha.pJS", "mean"])),
  3
) # Mean
round(
  exp(pradel.Toreadicus.df["alpha.pJS", "2.5%"]) /
    (1 + exp(pradel.Toreadicus.df["alpha.pJS", "2.5%"])),
  3
) # Low CI
round(
  exp(pradel.Toreadicus.df["alpha.pJS", "97.5%"]) /
    (1 + exp(pradel.Toreadicus.df["alpha.pJS", "97.5%"])),
  3
) # High CI

# Mean of Recruitment
pradel.Toreadicus.df["alpha.f", "mean"]
round(exp(pradel.Toreadicus.df["alpha.f", "mean"]), 3) # Mean
round(exp(pradel.Toreadicus.df["alpha.f", "2.5%"]), 3) # Low CI
round(exp(pradel.Toreadicus.df["alpha.f", "97.5%"]), 3) # High CI

f.df <- PJS.df[PJS.df$param == "f", ]
round(tapply(f.df$mean, f.df$year, sum), 3)
round(tapply(f.df$low2.5, f.df$year, sum), 3)
round(tapply(f.df$hi97.5, f.df$year, sum), 3)
round(tapply(f.df$low.25, f.df$year, sum), 3)
round(tapply(f.df$hi75, f.df$year, sum), 3)

tapply(f.df$mean, f.df$year, mean)
tapply(f.df$mean, f.df$year, function(x) {
  sd(x) / mean(x)
})

tapply(f.df$mean, f.df$month, mean)
tapply(f.df$mean, f.df$month, function(x) {
  sd(x) / mean(x)
})


# GOF tests ---------------------------------------------------------------

overall_CJS(as.matrix(bugs.data$y), rep(1, nrow(bugs.data$y)))

# Test 3.SR - transients
# H0 = Newly encountered individuals have the same chance to be later reobserved  as  recaptured  (previously  encountered)  individuals
test3sr(as.matrix(bugs.data$y), rep(1, nrow(bugs.data$y)))$test3sr

# Test 2.CT - trap dependence
# H0 = Missed individuals have the same chance to be recaptured at the next occasion as currently captured individuals
test2ct(as.matrix(bugs.data$y), rep(1, nrow(bugs.data$y)))$test2ct

# Test 3.SM
# H0 = Among those individuals seen again, when they were seen does not differ among previously and newly marked individuals
test3sm(as.matrix(bugs.data$y), rep(1, nrow(bugs.data$y)))$test3sm

# Test 2.CL
# H0 = There is no difference in the timing of reencounters between the individuals encountered and not encountered at occasion i, conditional on presence at both occasions i and i +2
test2cl(as.matrix(bugs.data$y), rep(1, nrow(bugs.data$y)))$test2cl


# Sensitivity analysis ----------------------------------------------------

(sens.phi <- colMeans(
  pradel.Toreadicus$sims.list$phiJS[, -46] /
    pradel.Toreadicus$sims.list$rho[, -1]
))

(sens.phi2.5 <- apply(
  (pradel.Toreadicus$sims.list$phiJS[, -46] /
    pradel.Toreadicus$sims.list$rho[, -1]),
  2,
  quantile,
  probs = 0.025
))
(sens.phi97.5 <- apply(
  (pradel.Toreadicus$sims.list$phiJS[, -46] /
    pradel.Toreadicus$sims.list$rho[, -1]),
  2,
  quantile,
  probs = 0.975
))


(sens.f <- colMeans(
  pradel.Toreadicus$sims.list$f[, -46] / pradel.Toreadicus$sims.list$rho[, -1]
))

(sens.f2.5 <- apply(
  (pradel.Toreadicus$sims.list$f[, -46] /
    pradel.Toreadicus$sims.list$rho[, -1]),
  2,
  quantile,
  probs = 0.025
))
(sens.f97.5 <- apply(
  (pradel.Toreadicus$sims.list$f[, -46] /
    pradel.Toreadicus$sims.list$rho[, -1]),
  2,
  quantile,
  probs = 0.975
))

round(mean(sens.phi), 3) * 100
round(mean(sens.f), 3) * 100


tapply(sens.phi, phi.df$year, mean)
tapply(sens.f, f.df$year, mean)

sens.df <- data.frame(
  sens.phi,
  sens.phi2.5,
  sens.phi97.5,
  sens.f,
  sens.f2.5,
  sens.f97.5,
  month = phi.df$month,
  time = phi.df$time
)

boxplot(scale(mean) ~ month, data = phi.df)
boxplot(scale(mean) ~ month, data = f.df)

ggplot(sens.df, aes(time, sens.phi)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = sens.phi), size = 1.5, col = "gold") +
  geom_line(aes(y = sens.f), size = 1.5, col = "blue") +
  scale_x_continuous(
    breaks = c(1:45),
    labels = meses[-46],
    expand = c(0, 0.6)
  ) +
  geom_ribbon(
    aes(ymin = sens.phi2.5, ymax = sens.phi97.5),
    fill = "gold",
    alpha = .2
  ) +
  geom_ribbon(
    aes(ymin = sens.f2.5, ymax = sens.f97.5),
    fill = "blue",
    alpha = .2
  ) +
  labs(x = NULL, y = "Elasticity") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )


# Diagnostics
S <- ggs(pradel.Toreadicus$samples[, c(47:64, 110:127, 174:191, 237)])
str(S)
levels(S$Parameter)

ggmcmc(S)

# SVL analyses ------------------------------------------------------------

# Load data
Toreadicus.data <- readRDS("Toreadicus_PEL_imp.rds")

# Total captures and recaptures
length(Toreadicus.data$svl)


# SVL temporal variation plots

(meses <- c(
  rep(
    c(
      "Feb",
      "Mar",
      "Apr",
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Dec",
      "Jan"
    ),
    3
  ),
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct",
  "Nov"
))

head(data)

plot(
  svl ~ campaign,
  data = Toreadicus.data,
  ylab = "Snout-vent length (mm)",
  col = rgb(0, 0, 0, .25),
  pch = 19,
  axes = F,
  cex = 1.5,
  bty = "n",
  ylim = c(20, 100),
  xlab = "Months"
)
axis(1, 1:46, labels = meses)
axis(2)
abline(h = 50, lty = 3)

boxplot(
  svl ~ month,
  data = Toreadicus.data,
  axes = F,
  ylab = "Snout-vent length (mm)",
  ylim = c(20, 100),
  xlab = "Months"
)
axis(
  1,
  1:12,
  labels = c(
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
  )
)
axis(2)
abline(h = 50, lty = 3)

# SVL between sexes
head(data)
boxplot(svl ~ sex, data = Toreadicus.data, ylab = "Snout-vent length (mm)")


cap.month <- function(x4, recap) {
  cap <- subset(x4, recapture == recap)
  cap.month <- as.data.frame(table(cap$month.n))
  names(cap.month) <- c("month", "captures")

  month.n <- data.frame(1:max(x4$month.n, na.rm = T))
  names(month.n) <- c("month")

  cap.month.1 <- merge(month.n, cap.month, all.x = T)
  cap.month.1$captures[is.na(cap.month.1$captures)] <- 0
  cap.month.1

  recap <- subset(x4, recapture != recap)
  recap.month <- as.data.frame(table(recap$month.n))
  names(recap.month) <- c("month", "recaptures")

  recap.month.1 <- merge(month.n, recap.month, all.x = T)
  recap.month.1$recaptures[is.na(recap.month.1$recaptures)] <- 0
  recap.month.1

  cap.recap.months <- cbind(cap.month.1, recap.month.1[, -1])
  names(cap.recap.months) <- c("month", "captures", "recaptures")
  cap.recap.months
}


(cap.recap.months.To.PEL <- cap.month(Toreadicus.data, "N"))
barplot(
  as.matrix(t(cap.recap.months.To.PEL[, -1])),
  names.arg = meses,
  legend.text = T,
  ylim = c(0, 50),
  ylab = "Frequency",
  xlab = "Months",
  las = 1,
  args.legend = list(
    bty = "n",
    x = 50,
    y = 50,
    legend = c("Recaptures", "Captures"),
    cex = 1
  ),
  main = expression(italic("Tropidurus oreadicus") - PEL)
)

arms.caps <- function(x4, narms) {
  tab.caps <- as.data.frame(table(x4$trap))
  names(tab.caps) <- c("trap", "captures") # Rename columns
  armadilha.n <- data.frame(1:narms)
  names(armadilha.n) <- c("trap")
  arms.caps.sp <- merge(armadilha.n, tab.caps, all.x = T)
  arms.caps.sp$captures[is.na(arms.caps.sp$captures)] <- 0
  print(arms.caps.sp)
}

arms.caps.To.PEL <- arms.caps(Toreadicus.data, narms = 25)
summary(arms.caps.To.PEL)

barplot(
  arms.caps.To.PEL$captures,
  names.arg = 1:25,
  col = c(
    rep("#2c7bb6", 6),
    rep("#abd9e9", 6),
    rep("#fdae61", 3),
    rep("#d7191c", 5),
    rep("#fdae61", 5)
  ),
  las = 1,
  ylab = "Captures",
  xlab = "Trap",
  ylim = c(0, 70),
  main = expression(italic("Tropidurus oreadicus") - PEL)
)
text(4, 68, "Forest")
text(11, 68, "Cerradão")
text(16, 68, "Cerrado sensu stricto")
text(21, 68, "Open Cerrado")
text(27, 68, "Cerrado sensu stricto")

# Using INLA to relate captures with environmental variables-------------------------------------------------------------------------

## INLA------------------------------------------------------------------------
# Rectify the camapaign number of one individual
Toreadicus.df.XY$campaign[Toreadicus.df.XY$sampling_day > 6] <- 43

Toreadicus.df.XY <- Toreadicus.df.XY %>%
  group_by(campaign) %>%
  mutate(sampling_day = dense_rank(date)) %>% # Rank the dates within each month
  ungroup()

summary(Toreadicus.df.XY$sampling_day)

Toreadicus.captures.day <- as.data.frame(table(
  Toreadicus.df.XY$trap,
  Toreadicus.df.XY$campaign,
  Toreadicus.df.XY$sampling_day
))
names(Toreadicus.captures.day) <- c("trap", "t", "sampling_day", "freq")

Toreadicus.captures.day <- Toreadicus.captures.day %>%
  pivot_wider(
    names_from = sampling_day,
    values_from = freq,
    names_prefix = "day_"
  )

summary(Toreadicus.captures.day)

Toreadicus.captures$x <- rep(traps.xy$x, 46)
Toreadicus.captures$y <- rep(traps.xy$y, 46)

Toreadicus.captures.day.env <- left_join(
  Toreadicus.captures.day,
  env.vars.month,
  by = c("trap", "t")
)

summary(Toreadicus.captures.day.env)

pts.traps.PEL$x_km <- pts.traps.PEL$x / 1000
pts.traps.PEL$y_km <- pts.traps.PEL$y / 1000

Toreadicus.captures.day.env$t <- as.integer(Toreadicus.captures.day.env$t)
Toreadicus.captures.day.env$x_km <- Toreadicus.captures.day.env$x / 1000
Toreadicus.captures.day.env$y_km <- Toreadicus.captures.day.env$y / 1000

(y.mat <- as.matrix(Toreadicus.captures.day.env[, c(3:8)]))
counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  scale(Toreadicus.captures.day.env$precip),
  scale(Toreadicus.captures.day.env$Toreadicus_perf)
)

print(counts.and.count.covs)

## Run model with Additive effects of precip and perf on lambda and p, no spatio-temporal random effects --------

out.inla.0 <- inla(
  counts.and.count.covs[, c(1:7)] ~ 1,
  data = list(counts.and.count.covs = counts.and.count.covs[, c(1:7)]),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(
        prior = "flat",
        param = numeric()
      )
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE
  )
)
summary(out.inla.0, digits = 3)


plot(
  out.inla.0,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.0,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.0,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)


## Run model with Additive effects of precip and perf on lambda and p, no spatio-temporal random effects-------------------------------------------------------------------------

out.inla.1 <- inla(
  counts.and.count.covs ~ 1 + x1.p + x2.p,
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    x1.p = scale(Toreadicus.captures.day.env$precip),
    x2.p = scale(Toreadicus.captures.day.env$Toreadicus_perf)
  ),
  family = "nmixnb",
  # Priors for detection parameters
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  # Priors for abundance and overdispersion parameters
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(param = c(0, 0.01)),
      theta3 = list(param = c(0, 0.01)),
      theta5 = list(
        prior = "flat",
        param = numeric()
      )
    )
  ),
  verbose = TRUE,
  control.inla = list(int.strategy = "eb"),
  control.compute = list(
    config = TRUE,
    waic = TRUE
  )
)
summary(out.inla.1, digits = 3)


plot(
  out.inla.1,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.1,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = TRUE,
  plot.q = F,
  plot.cpo = F
)

# Create SPDE mesh for spatial effects (traps along a transect)
coords <- cbind(pts.traps.PEL$x, pts.traps.PEL$y)
bnd <- inla.nonconvex.hull(traps.pts.PEL.utm, convex = 100)
mesh.inla <- inla.mesh.2d(
  boundary = bnd,
  traps.pts.PEL.utm,
  max.edge = c(30, 100),
  cutoff = 30
)

par(mfrow = c(1, 1))
plot(mesh.inla)
points(traps.pts.PEL.utm, col = "red", pch = 19)

spde <- inla.spde2.pcmatern(
  mesh = mesh.inla,
  prior.range = c(300, 0.5),
  prior.sigma = c(0.5, 0.01)
)

A <- inla.spde.make.A(
  mesh.inla,
  loc = cbind(
    Toreadicus.captures.day.env$x,
    Toreadicus.captures.day.env$y
  ),
  group = Toreadicus.captures.day.env$t
)

## Run model with Additive effects of precip and perf on lambda and p, smooth temporal and spatio-temporal effects --------

out.inla.2 <- inla(
  counts.and.count.covs ~
    1 +
      x1.p +
      x2.p +
      f(
        spatial,
        model = spde,
        A.local = A,
        group = Toreadicus.captures.day.env$t,
        control.group = list(
          model = "ar1",
          hyper = list(
            rho = list(
              prior = "pc.cor0",
              param = c(0.5, 0.3)
            )
          )
        )
      ) +
      f(t, model = "rw2", cyclic = T),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    x1.p = scale(Toreadicus.captures.day.env$precip),
    x2.p = scale(Toreadicus.captures.day.env$Toreadicus_perf),
    spatial = rep(NA, nrow(Toreadicus.captures.day.env)),
    t = Toreadicus.captures.day.env$t
  ),
  family = "nmixnb",
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(param = c(0, 0.01)),
      theta3 = list(param = c(0, 0.01)),
      theta4 = list(prior = "gamma", param = c(2, 0.5)),
      theta5 = list(param = c(0, 0.01)),
      theta6 = list(param = c(0, 0.01)),
      theta7 = list(param = c(0, 0.01))
    )
  ),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.compute = list(
    config = TRUE,
    waic = TRUE
  )
)
summary(out.inla.2, digits = 3)

plot(
  out.inla.2,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.2,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)


## Run model with Additive effects of precip and perf on lambda and p, no temporal effect, only spatio-temporal random effects -------------------------------------------------------------------------

out.inla.3 <- inla(
  counts.and.count.covs ~
    1 +
      x1.p +
      x2.p +
      f(
        spatial,
        model = spde,
        A.local = A,
        group = Toreadicus.captures.day.env$t,
        control.group = list(
          model = "ar1",
          hyper = list(
            rho = list(
              prior = "pc.cor0",
              param = c(0.5, 0.3)
            )
          )
        )
      ),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    x1.p = scale(Toreadicus.captures.day.env$precip),
    x2.p = scale(Toreadicus.captures.day.env$Toreadicus_perf),
    spatial = rep(NA, nrow(Toreadicus.captures.day.env)),
    t = Toreadicus.captures.day.env$t
  ),
  family = "nmixnb",
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(param = c(0, 0.01)),
      theta3 = list(param = c(0, 0.01)),
      theta4 = list(prior = "gamma", param = c(2, 0.5)),
      theta5 = list(param = c(0, 0.01)),
      theta6 = list(param = c(0, 0.01)),
      theta7 = list(param = c(0, 0.01))
    )
  ),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.compute = list(
    config = TRUE,
    waic = TRUE
  )
)

summary(out.inla.3, digits = 3)

plot(
  out.inla.3,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.3,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

## Run model with Interaction between precip and perf only for p -----------

out.inla.4 <- inla(
  counts.and.count.covs ~
    1 +
      x1.p +
      x2.p +
      x1.p:x2.p +
      f(
        spatial,
        model = spde,
        A.local = A,
        group = Toreadicus.captures.day.env$t,
        control.group = list(
          model = "ar1",
          hyper = list(
            rho = list(
              prior = "pc.cor0",
              param = c(0.5, 0.3)
            )
          )
        )
      ),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    x1.p = scale(Toreadicus.captures.day.env$precip),
    x2.p = scale(Toreadicus.captures.day.env$Toreadicus_perf),
    spatial = rep(NA, nrow(Toreadicus.captures.day.env)),
    t = Toreadicus.captures.day.env$t
  ),
  family = "nmixnb",
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(param = c(0, 0.01)),
      theta3 = list(param = c(0, 0.01)),
      theta4 = list(prior = "loggamma", param = c(2, 0.5))
    )
  ),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)

summary(out.inla.4, digits = 3)

out.inla.4$residuals$deviance.residuals
table(out.inla.4$cpo$failure)

plot(
  out.inla.4,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.4,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = T
)


## Run model with No environmental variables, only spatio-temporal random effects -------

out.inla.5 <- inla(
  counts.and.count.covs[, c(1:7)] ~
    1 +
      f(
        spatial,
        model = spde,
        A.local = A,
        group = Toreadicus.captures.day.env$t,
        control.group = list(
          model = "ar1",
          hyper = list(
            rho = list(
              prior = "pc.cor0",
              param = c(0.5, 0.3)
            )
          )
        )
      ),
  data = list(
    counts.and.count.covs = counts.and.count.covs[, c(1:7)],
    spatial = rep(NA, nrow(Toreadicus.captures.day.env)),
    t = Toreadicus.captures.day.env$t
  ),
  family = "nmixnb",
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(prior = "gamma", param = c(2, 0.5))
    )
  ),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)

summary(out.inla.5, digits = 3)

plot(
  out.inla.5,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.5,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.5,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = T
)


counts.and.count.covs <- inla.mdata(
  y.mat,
  1,
  scale(Toreadicus.captures.day.env$precip),
  scale(Toreadicus.captures.day.env$Toreadicus_perf),
  scale(Toreadicus.captures.day.env$precip) *
    scale(Toreadicus.captures.day.env$Toreadicus_perf)
)

print(counts.and.count.covs)

## Run model with Interaction between precip and perf for N and p --------

out.inla.6 <- inla(
  counts.and.count.covs ~
    1 +
      x1.p +
      x2.p +
      x1.p:x2.p +
      f(
        spatial,
        model = spde,
        A.local = A,
        group = Toreadicus.captures.day.env$t,
        control.group = list(
          model = "ar1",
          hyper = list(
            rho = list(
              prior = "pc.cor0",
              param = c(0.5, 0.3)
            )
          )
        )
      ),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    x1.p = scale(Toreadicus.captures.day.env$precip),
    x2.p = scale(Toreadicus.captures.day.env$Toreadicus_perf),
    spatial = rep(NA, nrow(Toreadicus.captures.day.env)),
    t = Toreadicus.captures.day.env$t
  ),
  family = "nmixnb",
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(param = c(0, 0.01)),
      theta3 = list(param = c(0, 0.01)),
      theta4 = list(param = c(0, 0.01)),
      theta5 = list(prior = "gamma", param = c(2, 0.5)),
      theta6 = list(param = c(0, 0.01)),
      theta7 = list(param = c(0, 0.01)),
      theta8 = list(param = c(0, 0.01)),
      theta9 = list(param = c(0, 0.01))
    )
  ),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.compute = list(
    config = TRUE,
    waic = TRUE
  )
)

summary(out.inla.6, digits = 3)

plot(
  out.inla.6,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.6,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)


## Run model with Interaction between precip and perf for N and additive for p --------

out.inla.7 <- inla(
  counts.and.count.covs ~
    1 +
      x1.p +
      x2.p +
      f(
        spatial,
        model = spde,
        A.local = A,
        group = Toreadicus.captures.day.env$t,
        control.group = list(
          model = "ar1",
          hyper = list(
            rho = list(
              prior = "pc.cor0",
              param = c(0.5, 0.3)
            )
          )
        )
      ),
  data = list(
    counts.and.count.covs = counts.and.count.covs,
    x1.p = scale(Toreadicus.captures.day.env$precip),
    x2.p = scale(Toreadicus.captures.day.env$Toreadicus_perf),
    spatial = rep(NA, nrow(Toreadicus.captures.day.env)),
    t = Toreadicus.captures.day.env$t
  ),
  family = "nmixnb",
  control.fixed = list(
    mean = 0,
    mean.intercept = 0,
    prec = 0.01,
    prec.intercept = 0.01
  ),
  control.family = list(
    hyper = list(
      theta1 = list(param = c(0, 0.01)),
      theta2 = list(param = c(0, 0.01)),
      theta3 = list(param = c(0, 0.01)),
      theta4 = list(param = c(0, 0.01)),
      theta5 = list(prior = "gamma", param = c(2, 0.5)),
      theta6 = list(param = c(0, 0.01)),
      theta7 = list(param = c(0, 0.01)),
      theta8 = list(param = c(0, 0.01)),
      theta9 = list(param = c(0, 0.01))
    )
  ),
  control.inla = list(int.strategy = "eb"),
  verbose = TRUE,
  control.compute = list(
    config = TRUE,
    waic = TRUE,
    residuals = TRUE,
    cpo = TRUE
  )
)

summary(out.inla.7, digits = 3)

plot(
  out.inla.7,
  plot.fixed.effects = TRUE,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.7,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = TRUE,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = F
)

plot(
  out.inla.7,
  plot.fixed.effects = F,
  plot.lincomb = F,
  plot.random.effects = F,
  plot.hyperparameters = F,
  plot.predictor = F,
  plot.q = F,
  plot.cpo = TRUE
)


## Evaluate models ---------------------------------------------------------

out.inla.0$waic$waic
out.inla.1$waic$waic
out.inla.2$waic$waic
out.inla.3$waic$waic
out.inla.4$waic$waic
out.inla.5$waic$waic
out.inla.6$waic$waic
out.inla.7$waic$waic


inla_models <- list(
  out.inla.0,
  out.inla.1,
  out.inla.2,
  out.inla.3,
  out.inla.4,
  out.inla.5,
  out.inla.6,
  out.inla.7
)

waic_values <- sapply(inla_models, function(model) model$waic$waic)
peff_values <- sapply(inla_models, function(model) model$waic$p.eff)
mlik_values <- sapply(inla_models, function(model) model$mlik[1])


delta_waic <- waic_values - min(waic_values)

w_akaike <- exp(-0.5 * delta_waic) / sum(exp(-0.5 * delta_waic))

waic_table <- data.frame(
  Model = paste0("Model ", 0:7),
  mlik = mlik_values,
  WAIC = waic_values,
  deltaWAIC = delta_waic,
  pWAIC = peff_values,
  Weight = round(w_akaike, 2)
)

waic_table$formula[1] <- "(lambda ~ 1), (p ~ 1)"
waic_table$formula[2] <- "(lambda ~ precip + perf), (p ~ precip + perf)"
waic_table$formula[
  3
] <- "(lambda ~ precip + perf), (p ~ precip + perf + spatiotemporal + temporal)"
waic_table$formula[
  4
] <- "(lambda ~ precip + perf), (p ~ precip + perf + spatiotemporal)"
waic_table$formula[
  5
] <- "(lambda ~ precip + perf), (p ~ precip * perf + spatiotemporal)"
waic_table$formula[6] <- "(lambda ~ 1), (p ~ spatiotemporal)"
waic_table$formula[
  7
] <- "(lambda ~ precip * perf), (p ~ precip * perf + spatiotemporal)"
waic_table$formula[
  8
] <- "(lambda ~ precip * perf), (p ~ precip + perf + spatiotemporal)"


# Sort by WAIC (best model first)
waic_table <- waic_table[order(waic_table$WAIC), ]

# Print table
print(waic_table)

summary(out.inla.7)

summary(out.inla.4)
summary(out.inla.3)


summary(plogis(as.matrix(out.inla.7$summary.linear.predictor[, c(
  "mean",
  "sd",
  "0.025quant",
  "0.5quant",
  "0.975quant",
  "mode"
)])))
summary(out.inla.7$summary.fitted.values)

save(
  out.inla.0,
  out.inla.1,
  out.inla.2,
  out.inla.3,
  out.inla.4,
  out.inla.5,
  out.inla.6,
  out.inla.7,
  file = "inla_models.RData"
)


## Diagnostics ------------------------------------------------------------

# Compute fitted values
fitted_values <- out.inla.7$summary.fitted.values$`0.5quant`

y_long <- as.numeric(y.mat) # Flatten matrix column-wise
fitted_long <- rep(fitted_values, ncol(y.mat)) # Repeat fitted values for each survey
dispersion(
  observed = y_long,
  fitted = fitted_long,
  variance = fitted_long
)

dispersion(
  observed = rowMeans(y.mat),
  fitted = fitted_values,
  variance = fitted_values
)

df_diagnostics <- data.frame(
  t = Toreadicus.captures.day.env$t,
  month = Toreadicus.captures.day.env$month.GMT_3,
  year = Toreadicus.captures.day.env$year.GMT_3,
  trap = Toreadicus.captures.day.env$trap,
  res = out.inla.7$residuals$deviance.residuals,
  cpo = out.inla.7$cpo$cpo
)

summary(df_diagnostics$res)

summary(df_diagnostics$cpo)

# RMSE
sqrt(mean(df_diagnostics$res)^2)

# MAE
mean(abs(df_diagnostics$res))

# cpo
-sum(log(df_diagnostics$cpo))

ggplot(df_diagnostics, aes(x = t, group = t, y = res)) +
  geom_boxplot() +
  labs(x = "Occasions", y = "Residuals") +
  theme_classic()

ggplot(df_diagnostics, aes(x = month, group = month, y = res)) +
  geom_boxplot() +
  labs(x = "Months", y = "Residuals") +
  scale_x_continuous(breaks = 1:12) +
  theme_classic()

ggplot(df_diagnostics, aes(x = year, group = year, y = res)) +
  geom_boxplot() +
  labs(x = "Years", y = "Residuals") +
  scale_x_continuous(breaks = 2018:2021) +
  theme_classic()

ggplot(df_diagnostics, aes(x = trap, group = trap, y = res)) +
  geom_boxplot() +
  labs(x = "Traps", y = "Residuals") +
  # scale_x_continuous(breaks = 2018:2021)+
  theme_classic()

ggplot(df_diagnostics, aes(x = t, group = t, y = cpo)) +
  geom_boxplot() +
  labs(x = "Occasions", y = "CPO-values") +
  theme_classic()

ggplot(df_diagnostics, aes(x = month, group = month, y = cpo)) +
  geom_boxplot() +
  labs(x = "Months", y = "CPO-values") +
  scale_x_continuous(breaks = 1:12) +
  theme_classic()

ggplot(df_diagnostics, aes(x = year, group = year, y = cpo)) +
  geom_boxplot() +
  labs(x = "Years", y = "CPO-values") +
  scale_x_continuous(breaks = 2018:2021) +
  theme_classic()

ggplot(df_diagnostics, aes(x = trap, group = trap, y = cpo)) +
  geom_boxplot() +
  labs(x = "Traps", y = "CPO-values") +
  theme_classic()

par(mfrow = c(1, 1))
hist(df_diagnostics$cpo)

(theta <- 1 /
  out.inla.7$summary.hyperpar["overdispersion for NMixNB observations", "mean"])

plot(
  y_long,
  fitted_long
)

plot(
  log(y_long),
  log(fitted_long)
)

par(mfrow = c(1, 1))
plot(
  jitter(rowMeans(y.mat)),
  out.inla.7$residuals$deviance.residuals,
  bty = "n",
  main = "Residuals vs. Observed",
  xlab = "Observed Values",
  ylab = "Residuals"
)
abline(h = 0, col = "red", lty = 2)

plot(
  fitted_values,
  out.inla.7$residuals$deviance.residuals,
  bty = "n",
  main = "Residuals vs. Fitted",
  xlab = "Fitted Values",
  ylab = "Residuals"
)
abline(h = 0, col = "red", lty = 2)

plot(
  fitted_values,
  out.inla.7$cpo$cpo,
  bty = "n",
  main = "CPO vs. Fitted",
  xlab = "Fitted Values",
  ylab = "CPO-values"
)

plot(
  jitter(rowMeans(y.mat)),
  out.inla.7$cpo$cpo,
  bty = "n",
  main = "CPO vs. Observed",
  xlab = "Observed Values",
  ylab = "CPO-values"
)

plot(
  jitter(rowMeans(y.mat)),
  fitted_values,
  bty = "n",
  main = "Fitted vs. Observed",
  xlab = "Observed Values",
  ylab = "Fitted values"
)

qqnorm(out.inla.7$residuals$deviance.residuals)

hist(out.inla.7$residuals$deviance.residuals, xlab = "Residuals", main = "")

# Generate theoretical quantiles from a negative binomial
theoretical_nb <- qnbinom(
  ppoints(length(out.inla.4$residuals$deviance.residuals)),
  mu = mean(y.mat),
  size = theta
)

# QQ-plot for negative binomial residuals
qqplot(
  theoretical_nb,
  out.inla.7$residuals$deviance.residuals,
  main = "QQ-Plot vs Negative Binomial",
  xlab = "Theoretical Quantiles",
  ylab = "Pearson Residuals"
) +
  theme_classic()


# Predictive R² (squared Pearson correlation)
cor(rowMeans(y.mat), fitted_values)^2
cor(rowMeans(y.mat), fitted_values, method = "spearman")^2


out.spatial <- inla.spde2.result(
  inla = out.inla.7,
  name = "spatial",
  spde = spde,
  do.transf = TRUE
)

# Extract posterior means for time t = 1

inla.nmix.lambda.fitted <- function(
  result,
  sample.size = 1000,
  return.posteriors = FALSE,
  scale = "exp"
) {
  fam <- result$.args$family
  if (length(grep(pattern = "nmix", x = fam)) == 0) {
    stop("This function is only for models with 'nmix' or 'nmixnb' likelihoods")
  }
  if (missing(result)) {
    stop("Please specify a model result")
  }
  s.check <- as.numeric(scale == "exp") + as.numeric(scale == "log")
  if (s.check == 0) {
    stop("Scale must be set to 'exp' or 'log'")
  }
  if (sample.size < 500) {
    warning("Please increase the sample size")
  }
  mdata.obj <- result$.args$data[[1]]
  counts <- as.data.frame(mdata.obj[grep(pattern = "Y", names(mdata.obj))])
  lambda.covs <- as.data.frame(mdata.obj[grep(
    pattern = "X",
    names(mdata.obj)
  )])
  lambda.covs <- as.matrix(lambda.covs)
  n.counts <- ncol(counts)
  n.lambda.covs <- ncol(lambda.covs)
  n.data <- nrow(counts)
  hyperpar.samples <- inla.hyperpar.sample(sample.size, result)
  s.names <- rownames(hyperpar.samples)
  if (fam == "nmixnb") {
    hyperpar.samples <- hyperpar.samples[, -(5:ncol(hyperpar.samples))]
  }
  n.samp.covs <- ncol(hyperpar.samples)
  fitted.posteriors <- matrix(-1, nrow = n.data, ncol = sample.size)
  for (i in 1:n.data) {
    obs <- lambda.covs[i, ]
    for (j in 1:sample.size) {
      post <- hyperpar.samples[j, ]
      fitted <- sum(obs * post)
      fitted.posteriors[i, j] <- fitted
    }
  }
  index <- 1:n.data
  fitted.posteriors <- as.data.frame(fitted.posteriors)
  row.names(fitted.posteriors) <- NULL
  names(fitted.posteriors) <- s.names
  if (scale == "exp") {
    fitted.posteriors <- exp(fitted.posteriors)
  }
  fitted.meds <- round(
    apply(fitted.posteriors, 1, median),
    4
  )
  fitted.means <- round(
    apply(fitted.posteriors, 1, mean),
    4
  )
  fitted.sds <- round(apply(fitted.posteriors, 1, sd), 4)
  fitted.q025 <- round(apply(fitted.posteriors, 1, quantile, probs = 0.025), 4)
  fitted.q500 <- round(apply(fitted.posteriors, 1, quantile, probs = 0.5), 4)
  fitted.q975 <- round(apply(fitted.posteriors, 1, quantile, probs = 0.975), 4)
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  fitted.modes <- round(
    apply(fitted.posteriors, 1, Mode),
    4
  )
  fitted.summary <- data.frame(
    mean.lambda = fitted.means,
    sd.lambda = fitted.sds,
    quant025.lambda = fitted.q025,
    median.lambda = fitted.q500,
    quant975.lambda = fitted.q975,
    mode.lambda = fitted.modes
  )
  fitted.summary <- cbind(index, fitted.summary)
  fitted.posteriors <- cbind(index, fitted.posteriors)
  if (return.posteriors == TRUE) {
    out <- list(
      fitted.summary = fitted.summary,
      fitted.posteriors = fitted.posteriors
    )
  } else {
    out <- list(fitted.summary = fitted.summary)
  }
  return(out)
}

out.inla.7.lambda.fits <- inla.nmix.lambda.fitted(
  result = out.inla.7,
  sample.size = 10000,
  return.posteriors = FALSE
)$fitted.summary

summary(out.inla.7.lambda.fits)
table(round(out.inla.7.lambda.fits$median.lambda))
hist(round(out.inla.7.lambda.fits$median.lambda))

table(round(out.inla.7.lambda.fits$mean.lambda))
hist(round(out.inla.7.lambda.fits$mean.lambda))

table(round(out.inla.7.lambda.fits$mode.lambda))
hist(round(out.inla.7.lambda.fits$mode.lambda))


ggplot(
  stack(out.inla.7.lambda.fits[, c(
    "mean.lambda",
    "median.lambda",
    "mode.lambda"
  )]),
  aes(x = values, fill = ind)
) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "log(Abundance)") +
  theme_classic()

ggplot(
  stack(log(out.inla.7.lambda.fits[, c(
    "mean.lambda",
    "median.lambda",
    "mode.lambda"
  )])),
  aes(x = values, fill = ind)
) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "log(Abundance)") +
  theme_classic()


N <- mesh.inla$n
proj <- inla.mesh.projector(mesh.inla)

pred <- NULL
field_t <- NULL
for (i in 1:46) {
  ii <- (i - 1) * N + 1
  jj <- i * N
  t <- i
  pred[[i]] <- out.inla.7$summary.random$spatial$mean[ii:jj]

  # Create projector and project posterior means
  field_t[[i]] <- cbind(
    x = proj$lattice$loc[, 1],
    y = proj$lattice$loc[, 2],
    z = as.vector(inla.mesh.project(proj, pred[[i]])),
    t = i
  )
}

field_t <- do.call("rbind.data.frame", field_t) %>% filter(!is.na(z))
pred_df <- data.frame(
  t = Toreadicus.captures.day.env$t,
  month = Toreadicus.captures.day.env$month.GMT_3,
  precip = scale(Toreadicus.captures.day.env$precip),
  perf = scale(Toreadicus.captures.day.env$Toreadicus_perf),
  x = Toreadicus.captures.day.env$x,
  y = Toreadicus.captures.day.env$y,
  trap = Toreadicus.captures.day.env$trap,
  N = out.inla.7.lambda.fits$median.lambda,
  Nlow = out.inla.7.lambda.fits$quant025.lambda,
  Nhigh = out.inla.7.lambda.fits$quant975.lambda,
  p = plogis(out.inla.7$summary.linear.predictor$mean),
  ncaps = rowSums(y.mat)
)

summary(pred_df)

# Plot
ggplot(field_t, aes(x = x / 1000, y = y / 1000, fill = exp(z))) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10") +
  geom_point(
    data = pred_df,
    aes(x = x / 1000, y = y / 1000, colour = N),
    fill = NA
  ) +
  scale_color_viridis_c(trans = "log10") +
  facet_wrap(~t) +
  labs(x = "Latitude (km)", y = "Longitude (km)", fill = "Spatial field") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(field_t, aes(x = x / 1000, y = y / 1000, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c() +
  geom_point(
    data = pred_df,
    aes(x = x / 1000, y = y / 1000, colour = p),
    fill = NA
  ) +
  scale_color_viridis_c(trans = "log10") +
  facet_wrap(~t) +
  labs(x = "Latitude (km)", y = "Longitude (km)", fill = "Spatial field") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(pred_df, aes(x = t, y = trap, z = N)) +
  geom_raster(aes(fill = N)) +
  scale_fill_viridis_c(trans = "log10", name = "Abundance") +
  labs(x = "Occasions", y = "Traps", z = "Abundance") +
  theme_minimal()

ggplot(pred_df, aes(x = t, y = trap, z = p)) +
  geom_raster(aes(fill = p)) +
  scale_fill_viridis_c(trans = "log10", name = "Capture") +
  labs(x = "Occasions", y = "Traps", z = "Capture") +
  theme_minimal()

ggplot(pred_df, aes(x = t, y = trap, z = ncaps)) +
  geom_raster(aes(fill = ncaps)) +
  scale_fill_viridis_c(name = "Number of captures", transform = "log1p") +
  labs(x = "Occasions", y = "Traps", z = "Number of captures") +
  theme_minimal()


precip_xx <- seq(min(pred_df$precip), max(pred_df$precip), by = 0.1)
perf_xx <- seq(min(pred_df$perf), max(pred_df$perf), by = 0.1)

rel_p_precip <- data.frame(precip = precip_xx)

rel_p_precip$p <- plogis(
  out.inla.7$summary.fixed$`0.5quant`[1] +
    out.inla.7$summary.fixed$`0.5quant`[2] * rel_p_precip$precip
)

rel_p_precip$plower <- plogis(
  (out.inla.7$summary.fixed$`0.975quant`[1]) +
    (out.inla.7$summary.fixed$`0.025quant`[2]) * rel_p_precip$precip
)

rel_p_precip$pupper <- plogis(
  (out.inla.7$summary.fixed$`0.025quant`[1]) +
    (out.inla.7$summary.fixed$`0.975quant`[2]) * rel_p_precip$precip
)

ggplot(pred_df, aes(x = precip, y = p)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_line(
    data = rel_p_precip,
    aes(x = precip, y = p),
    col = "blue",
    linewidth = 2
  ) +
  geom_ribbon(
    data = rel_p_precip,
    aes(ymin = plower, ymax = pupper),
    alpha = 0.3
  ) +
  labs(x = "Precipitation", y = "Capture") +
  scale_y_log10() +
  theme_classic()

rel_p_perf <- data.frame(perf = perf_xx)

rel_p_perf$p <- exp(
  out.inla.7$summary.fixed$`0.5quant`[1] +
    out.inla.7$summary.fixed$`0.5quant`[3] * rel_p_perf$perf
)

rel_p_perf$plower <- exp(
  (out.inla.7$summary.fixed$`0.025quant`[1]) +
    (out.inla.7$summary.fixed$`0.025quant`[3] * rel_p_perf$perf)
)

rel_p_perf$pupper <- exp(
  (out.inla.7$summary.fixed$`0.975quant`[1]) +
    (out.inla.7$summary.fixed$`0.975quant`[3] * rel_p_perf$perf)
)

ggplot(pred_df, aes(x = perf, y = p)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_line(
    data = rel_p_perf,
    aes(x = perf, y = p),
    col = "blue",
    linewidth = 2
  ) +
  geom_ribbon(
    data = rel_p_perf,
    aes(ymin = plower, ymax = pupper),
    alpha = 0.3
  ) +
  labs(x = "Locomotor performance", y = "Capture") +
  theme_classic()


rel_N <- expand.grid(
  precip = precip_xx,
  perf = perf_xx
)

rel_N$interaction <- rel_N$precip * rel_N$perf

rel_N$N <- exp(
  out.inla.7$summary.hyperpar$`0.5quant`[1] + # Intercept
    out.inla.7$summary.hyperpar$`0.5quant`[2] * rel_N$precip + # Precipitation effects
    out.inla.7$summary.hyperpar$`0.5quant`[3] * rel_N$perf + # Performance effects
    out.inla.7$summary.hyperpar$`0.5quant`[4] * rel_N$interaction
) # Interaction


ggplot(pred_df, aes(x = precip, y = N)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Precipitation", y = "Abundance") +
  theme_classic()

ggplot(pred_df, aes(x = perf, y = N)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Locomotor performance", y = "Abundance") +
  theme_classic()

ggplot(rel_N, aes(x = precip, y = perf, z = N)) +
  geom_raster(aes(fill = N)) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(trans = "log10", name = "Abundance") +
  labs(x = "Precipitation", y = "Locomotor performance", z = "Abundance") +
  theme_minimal()

# Heatmap plot with Abundance x traps x t ------------------------------------------
abund.locfit <- locfit(N ~ lp(month, trap), data = pred_df)
plot(abund.locfit)

plot(abund.locfit, type = "persp", theta = 35, phi = 30, shade = 0.5)

ys <- seq(1, 25, by = 1)
xs <- seq(1, 12, by = 1)

nrz <- length(xs)
ncz <- length(ys)
z <- predict(
  abund.locfit,
  newdata = expand.grid(
    month = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  )
)

zs <- matrix(
  z,
  nrow = length(seq(1, 12, by = 1)),
  ncol = length(seq(1, 25, by = 1))
)

pal <- viridis(20)

# Compute the z-value at the facet centres
zfacet <- zs[-1, -1] + zs[-1, -ncz] + zs[-nrz, -1] + zs[-nrz, -ncz]

facetcol <- cut(zfacet, 20)

persp(
  x = xs,
  y = ys,
  z = zs,
  col = pal[facetcol],
  theta = 45,
  phi = 30,
  xlab = "Traps",
  ylab = "Months",
  zlab = "Abundance"
)

df <- data.frame(
  expand.grid(
    month = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  ),
  abund = z
)

ggplot(df, aes(x = month, y = trap, z = abund)) +
  geom_raster(aes(fill = abund)) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(trans = "log10", name = "Abundance") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = 1:25) +
  labs(x = "Months", y = "Traps", z = "Abundance") +
  theme_minimal()


# Heatmap plot with Capture x traps x t ------------------------------------------
p.locfit <- locfit(p ~ lp(month, trap), data = pred_df)
plot(p.locfit)

plot(p.locfit, type = "persp", theta = 35, phi = 30, shade = 0.5)

ys <- seq(1, 25, by = 1)
xs <- seq(1, 12, by = 1)

nrz <- length(xs)
ncz <- length(ys)
z <- predict(
  p.locfit,
  newdata = expand.grid(
    month = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  )
)

zs <- matrix(
  z,
  nrow = length(seq(1, 12, by = 1)),
  ncol = length(seq(1, 25, by = 1))
)

pal <- viridis(20)

# Compute the z-value at the facet centres
zfacet <- zs[-1, -1] + zs[-1, -ncz] + zs[-nrz, -1] + zs[-nrz, -ncz]

facetcol <- cut(zfacet, 20)

persp(
  x = xs,
  y = ys,
  z = zs,
  col = pal[facetcol],
  theta = 45,
  phi = 30,
  xlab = "Traps",
  ylab = "Months",
  zlab = "Capture"
)

df <- data.frame(
  expand.grid(
    month = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  ),
  p = z
)

ggplot(df, aes(x = month, y = trap, z = p)) +
  geom_raster(aes(fill = p)) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(name = "Capture") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = 1:25) +
  labs(x = "Months", y = "Traps", z = "Capture") +
  theme_minimal()

# Heatmap plot with Locomotor performance x traps x t ------------------------------------------
perf.locfit <- locfit(
  Toreadicus_perf ~ lp(month.GMT_3, trap),
  data = Toreadicus.captures.day.env
)
plot(perf.locfit)

plot(perf.locfit, type = "persp", theta = 35, phi = 30, shade = 0.5)

ys <- seq(1, 25, by = 1)
xs <- seq(1, 12, by = 1)

nrz <- length(xs)
ncz <- length(ys)
z <- predict(
  perf.locfit,
  newdata = expand.grid(
    month.GMT_3 = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  )
)

zs <- matrix(
  z,
  nrow = length(seq(1, 12, by = 1)),
  ncol = length(seq(1, 25, by = 1))
)

pal <- viridis(20)

# Compute the z-value at the facet centres
zfacet <- zs[-1, -1] + zs[-1, -ncz] + zs[-nrz, -1] + zs[-nrz, -ncz]

facetcol <- cut(zfacet, 20)

persp(
  x = xs,
  y = ys,
  z = zs,
  col = pal[facetcol],
  theta = 45,
  phi = 30,
  xlab = "Traps",
  ylab = "Months",
  zlab = "Capture"
)

df <- data.frame(
  expand.grid(
    month = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  ),
  perf = z
)

ggplot(df, aes(x = month, y = trap, z = perf)) +
  geom_raster(aes(fill = perf)) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(name = "Locomotor performance") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = 1:25) +
  labs(x = "Months", y = "Traps", z = "Locomotor performance") +
  theme_minimal()


# Heatmap plot with N captures x traps x t ------------------------------------------
ncaps.locfit <- locfit(
  freq ~ lp(month.GMT_3, trap),
  data = Toreadicus.captures.env
)
plot(ncaps.locfit)

plot(ncaps.locfit, type = "persp", theta = 35, phi = 30, shade = 0.5)

ys <- seq(1, 25, by = 1)
xs <- seq(1, 12, by = 1)

nrz <- length(xs)
ncz <- length(ys)
z <- predict(
  ncaps.locfit,
  newdata = expand.grid(
    month.GMT_3 = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  )
)

zs <- matrix(
  z,
  nrow = length(seq(1, 12, by = 1)),
  ncol = length(seq(1, 25, by = 1))
)

pal <- viridis(20)

# Compute the z-value at the facet centres
zfacet <- zs[-1, -1] + zs[-1, -ncz] + zs[-nrz, -1] + zs[-nrz, -ncz]

facetcol <- cut(zfacet, 20)

persp(
  x = xs,
  y = ys,
  z = zs,
  col = pal[facetcol],
  theta = 45,
  phi = 30,
  xlab = "Traps",
  ylab = "Months",
  zlab = "Capture"
)

df <- data.frame(
  expand.grid(
    month = seq(1, 12, by = 1),
    trap = seq(1, 25, by = 1)
  ),
  ncaps = z
)

ggplot(df, aes(x = month, y = trap, z = ncaps)) +
  geom_raster(aes(fill = ncaps)) +
  geom_contour(colour = "white") +
  scale_fill_viridis_c(name = "Number of captures", transform = "log1p") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = 1:25) +
  labs(x = "Months", y = "Traps", z = "Number of captures") +
  theme_minimal()


out.inla.7.lambda.fits$t <- Toreadicus.captures.day.env$t
out.inla.7.lambda.fits$trap <- Toreadicus.captures.day.env$trap
summary(out.inla.7.lambda.fits)

summary(out.inla.7$summary.linear.predictor)
summary(out.inla.7$summary.fitted.values)


ggplot(out.inla.7.lambda.fits, aes(t, mode.lambda)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = mode.lambda), size = 1.5) +
  scale_x_continuous(breaks = c(1:46), labels = meses, expand = c(0, 0.6)) +
  labs(x = NULL, y = "Abundance") +
  facet_wrap(~trap) +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Density dependence ------------------------------------------------------

(abund <- tapply(Toreadicus.captures.env$freq, Toreadicus.captures.env$t, sum))

plot(
  abund,
  type = "b",
  pch = 19,
  bty = "n",
  axes = F,
  xlab = "",
  ylab = "Captures and recaptures",
  ylim = c(0, 50)
)
axis(1, 1:46, labels = meses)
axis(2)

pred_df_t <- pred_df %>%
  group_by(t) %>%
  summarise(
    N = mean(N),
    Nlow = mean(Nlow),
    Nhigh = mean(Nhigh),
    p = mean(p),
    ncaps = sum(ncaps)
  )


abund.df <- data.frame(
  time = 1:46,
  abund = abund,
  Nmean = abund / PJS.df$mean[PJS.df$param == "p"],
  Nlow2.5 = abund / PJS.df$low2.5[PJS.df$param == "p"],
  Nhi97.5 = abund / PJS.df$hi97.5[PJS.df$param == "p"],
  Nlow.25 = abund / PJS.df$low.25[PJS.df$param == "p"],
  Nhi75 = abund / PJS.df$hi75[PJS.df$param == "p"],
  Nmix = pred_df_t$N,
  Nmixlow = pred_df_t$Nlow,
  Nmixhigh = pred_df_t$Nhigh
)

plot(
  abund.df$Nmean[-46],
  phi.df$mean,
  ylab = "Survival",
  xlab = "Abundance",
  bty = "n",
  xlim = c(0, 200),
  ylim = c(0.6, 0.8),
  pch = 19,
  col = rgb(0, 0, 0, .5)
)
abline(v = 115, col = "red")

ddsurv.df <- data.frame(
  phi = phi.df$mean,
  N = abund.df$Nmean[-46],
  Nmix = abund.df$Nmix[-46]
)
m1 <- lm(phi ~ N + I(N^2), data = ddsurv.df)

summary(m1)
lines(
  1:200,
  predict(m1, newdata = data.frame(N = 1:200), type = "response"),
  col = "red"
)

plot(m1)

plot(log(abund.df$Nmean[-1]), log(rho.df$mean))
cor.test(log(abund.df$Nmean[-1]), log(rho.df$mean))
cor.test(log(abund.df$Nmean[-1]), log(rho.df$mean), method = "spearman")


plot(log(abund.df$Nmix[-1]), log(rho.df$mean))
cor.test(log(abund.df$Nmix[-1]), log(rho.df$mean))
cor.test(log(abund.df$Nmix[-1]), log(rho.df$mean), method = "spearman")

plot(log(abund.df$Nmix), PJS.df$mean[PJS.df$param == "p"])
cor.test(log(abund.df$Nmix), log(PJS.df$mean[PJS.df$param == "p"]))
cor.test(
  log(abund.df$Nmix),
  log(PJS.df$mean[PJS.df$param == "p"]),
  method = "spearman"
)


ddrecruit.df <- data.frame(
  f = f.df$mean,
  N = abund.df$Nmean[-46],
  Nmix = abund.df$Nmix[-46]
)

plot(
  log(f) ~ Nmix,
  data = ddrecruit.df,
  ylab = "Recruitment",
  xlab = "Abundance",
  bty = "n",
  ylim = c(-3, 2),
  pch = 19,
  col = rgb(0, 0, 0, .5)
)

cor.test(log(ddrecruit.df$f), ddrecruit.df$Nmix)
cor.test(log(ddrecruit.df$f), log(ddrecruit.df$Nmix))

cor.test(log(ddrecruit.df$f), log(ddrecruit.df$Nmix), method = "spearman")

plot(
  log(f) ~ N,
  data = ddrecruit.df,
  ylab = "Recruitment",
  xlab = "Abundance",
  bty = "n",
  xlim = c(0, 200),
  ylim = c(-3, 2),
  pch = 19,
  col = rgb(0, 0, 0, .5)
)

cor.test(log(ddrecruit.df$f), ddrecruit.df$N)
cor.test(log(ddrecruit.df$f), log(ddrecruit.df$N))

cor.test(log(ddrecruit.df$f), log(ddrecruit.df$N), method = "spearman")


m2 <- lm(log(f) ~ log(N), data = ddrecruit.df)
summary(m2)
lines(
  1:200,
  predict(m2, newdata = data.frame(N = 1:200), type = "response"),
  col = "red"
)

ggplot(ddrecruit.df, aes(N, log(f))) +
  stat_smooth(method = "lm", formula = y ~ log(x)) +
  geom_point(aes(N, log(f))) +
  labs(x = "Abundance", y = "log(recruitment)")

ggplot(ddrecruit.df, aes(log(Nmix), log(f))) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_point(aes(log(Nmix), log(f))) +
  labs(x = "Abundance", y = "log(recruitment)")

ggplot(ddrecruit.df, aes(log(Nmix), log(f))) +
  stat_smooth(method = "lm", formula = y ~ x) +
  geom_point(aes(log(Nmix), log(f))) +
  labs(x = "Abundance", y = "log(recruitment)")


ggplot(ddsurv.df, aes(N, phi)) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_point(aes(N, phi), alpha = 0.5) +
  labs(x = "Abundance", y = "Survival")

ggplot(ddsurv.df, aes(log(Nmix), phi)) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_point(aes(log(Nmix), phi), alpha = 0.5) +
  labs(x = "Abundance", y = "Survival")

# Comparison of abundance estimates ---------------------------------------
ggplot(
  stack(abund.df[, c("abund", "Nmean", "Nmix")]),
  aes(x = values, fill = ind)
) +
  geom_density(alpha = 0.5, col = NA) +
  labs(y = "Density", x = "Abundance values") +
  theme_classic()

ggplot(
  stack(abund.df[, c("abund", "Nmean", "Nmix")]),
  aes(x = values, fill = ind)
) +
  geom_density(alpha = 0.5, col = NA) +
  scale_x_continuous(transform = "log1p") +
  labs(y = "Density", x = "Abundance values") +
  theme_classic()

ggplot(abund.df, aes(time, Nmean)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = Nmean), size = 1.5) +
  scale_x_continuous(breaks = c(1:46), labels = meses, expand = c(0, 0.6)) +
  geom_ribbon(aes(ymin = Nlow2.5, ymax = Nhi97.5), alpha = .4) +
  geom_ribbon(aes(ymin = Nlow.25, ymax = Nhi75), alpha = .8) +
  # scale_y_log10()+
  # geom_point(aes(y=abund),size=3) +
  # geom_line(aes(y=abund),size=1) +
  # geom_line(aes(y = Nmix), size = 1, col = "blue", linetype = 2) +
  # geom_ribbon(aes(ymin = Nmixlow, ymax = Nmixhigh),alpha=.4, fill = "blue") +
  labs(x = NULL, y = "Abundance") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

ggplot(abund.df, aes(time, Nmean)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = Nmean), size = 1.5) +
  scale_x_continuous(breaks = c(1:46), labels = meses, expand = c(0, 0.6)) +
  scale_y_log10() +
  geom_ribbon(aes(ymin = Nlow2.5, ymax = Nhi97.5), alpha = .4) +
  geom_ribbon(aes(ymin = Nlow.25, ymax = Nhi75), alpha = .8) +
  geom_line(aes(y = Nmix), size = 1, col = "blue", linetype = 2) +
  geom_ribbon(aes(ymin = Nmixlow, ymax = Nmixhigh), alpha = .4, fill = "blue") +
  labs(x = NULL, y = "Abundance") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )


# Comparison between PJS and number of captures
ggplot(abund.df, aes(x = log(abund), y = log(Nmean))) +
  geom_point(alpha = 0.5) +
  labs(x = "log(Number of captures)", y = "log(Abundance)") +
  theme_classic()

cor.test(log(abund.df$abund), log(abund.df$Nmean))
cor.test(log(abund.df$abund), log(abund.df$Nmean), method = "spearman")


# Comparison between N-mixture and number of captures
ggplot(abund.df, aes(x = log(abund), y = log(Nmix))) +
  geom_point(alpha = 0.5) +
  labs(x = "log(Number of captures)", y = "log(Abundance)") +
  theme_classic()

cor.test(log(abund.df$abund), log(abund.df$Nmix))
cor.test(log(abund.df$abund), log(abund.df$Nmix), method = "spearman")

# Comparison between N-mixture and PJS
ggplot(abund.df, aes(x = log(Nmix), y = log(Nmean))) +
  geom_point(alpha = 0.5) +
  labs(x = "log(Abundance) - N-mixture", y = "log(Abundance) - PJS") +
  theme_classic()

cor.test(log(abund.df$Nmean), log(abund.df$Nmix))
cor.test(abund.df$Nmean, abund.df$Nmix)

cor.test(log(abund.df$Nmean), log(abund.df$Nmix), method = "spearman")

# Capture probability comparison between PJS and N-mixture models
p.df <- data.frame(
  time = 1:46,
  abund = abund,
  p = PJS.df$mean[PJS.df$param == "p"],
  pmix = pred_df_t$p
)

ggplot(p.df, aes(time, p)) +
  geom_rect(
    data = rect1,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect2,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect3,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect4,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_rect(
    data = rect5,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "#bdbdbd",
    inherit.aes = F
  ) +
  geom_line(aes(y = p), size = 1.5) +
  scale_x_continuous(breaks = c(1:46), labels = meses, expand = c(0, 0.6)) +
  geom_line(aes(y = pmix), size = 1, col = "blue", linetype = 2) +
  labs(x = NULL, y = "Capture") +
  theme(
    plot.margin = unit(c(0, 0, 2, 0), "lines"),
    text = element_text(size = 10),
    legend.position = c(0.1, 0.5)
  )

# Capture (PJS) x Number of captures
ggplot(p.df, aes(x = log(abund), y = p)) +
  geom_point(alpha = 0.5) +
  labs(x = "log(Number of captures)", y = "Capture - PJS") +
  theme_classic()

ggplot(p.df, aes(x = log(abund), y = pmix)) +
  geom_point(alpha = 0.5) +
  labs(x = "log(Number of captures)", y = "Capture - N-mixture") +
  theme_classic()

ggplot(p.df, aes(x = p, y = pmix)) +
  geom_point(alpha = 0.5) +
  labs(x = "Capture - PJS", y = "Capture - N-mixture") +
  theme_classic()

cor.test(p.df$p, p.df$pmix)
cor.test(log(p.df$p), log(p.df$pmix))
cor.test(log(p.df$p), log(p.df$pmix), method = "spearman")


# Temporal autocorrelation ------------------------------------------------

acfncaps <- acf(
  tapply(Toreadicus.captures.env$freq, Toreadicus.captures.env$t, sum),
  lag.max = 30,
  plot = F
)

plot(acfncaps, bty = "n", xlim = c(0, 30), main = "Captures")

acfnmean <- acf(log(abund.df$Nmean), lag.max = 30, plot = F)
which(acfncaps$acf[-c(1:2)] == max(acfncaps$acf[-c(1:2)])) + 1

plot(acfnmean, bty = "n", xlim = c(0, 30), main = "Abundance - PJS")
which(acfnmean$acf[-c(1:2)] == max(acfnmean$acf[-c(1:2)])) + 1


acfnmix <- acf(log(abund.df$Nmix), lag.max = 30, plot = F)

plot(acfnmix, bty = "n", xlim = c(0, 30), main = "Abundance - N-mixture")
which(acfnmix$acf[-c(1:2)] == max(acfnmix$acf[-c(1:2)])) + 1


ccf.nmixpjs <- ccf(log(abund.df$Nmean), log(abund.df$Nmix))
plot(
  ccf.nmixpjs,
  bty = "n",
  main = "Cross-correlation function between abundance estimates"
)

which.max(ccf.nmixpjs$acf)
ccf.nmixpjs$lag[which.max(ccf.nmixpjs$acf)]

cor(abund.df[, c(2, 3, 8)])
cor(log(abund.df[, c(2, 3, 8)]))

cor(abund.df[, c(2, 3, 8)], method = "spearman")


# Time-series decomposition -----------------------------------------------

## Abundance (PJS) ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_NPJS <- ts(log(abund.df$Nmean), start = c(2018, 2), freq = 12)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_NPJS <- RJDemetra::x13(ts_NPJS)
x13_NPJS

# Plots
par(mfrow = c(2, 1))
plot(x13_NPJS)
par(mfrow = c(1, 1))

# Residuals
par(mfrow = c(3, 2))
plot(x13_NPJS$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_NPJS$decomposition)

## Abundance (N-mixture) ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_Nmix <- ts(log(abund.df$Nmix), start = c(2018, 2), freq = 12)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_Nmix <- RJDemetra::x13(ts_Nmix)
x13_Nmix

# Plots
par(mfrow = c(2, 1))
plot(x13_Nmix)
par(mfrow = c(1, 1))

# Residuals
par(mfrow = c(3, 2))
plot(x13_Nmix$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_Nmix$decomposition)

## Population growth ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_rho <- ts(
  log(PJS.df[PJS.df$param == "rho", "mean"]),
  start = c(2018, 3),
  freq = 12
)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_rho <- RJDemetra::x13(ts_rho)
x13_rho

# Plots
par(mfrow = c(2, 1))
plot(x13_rho)
par(mfrow = c(1, 1))

# Residuals
par(mfrow = c(3, 2))
plot(x13_rho$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_rho$decomposition)

## Recruitment ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_f <- ts(
  log(PJS.df[PJS.df$param == "f", "mean"]),
  start = c(2018, 3),
  freq = 12
)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_f <- RJDemetra::x13(ts_f)
x13_f

# Plots
par(mfrow = c(2, 1))
plot(x13_f)
par(mfrow = c(1, 1))

# Residuals
par(mfrow = c(3, 2))
plot(x13_f$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_f$decomposition)

## Survival ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_phi <- ts(
  qlogis(PJS.df[PJS.df$param == "phi", "mean"]),
  start = c(2018, 3),
  freq = 12
)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_phi <- RJDemetra::x13(ts_phi)
x13_phi

# Plots
par(mfrow = c(2, 1))
plot(x13_phi)
par(mfrow = c(1, 1))

# Residuals

par(mfrow = c(3, 2))
plot(x13_phi$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_phi$decomposition)

## Capture (p) - PJS ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_pPJS <- ts(
  PJS.df[PJS.df$param == "p", "mean"],
  start = c(2018, 2),
  freq = 12
)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_pPJS <- RJDemetra::x13(ts_pPJS)
x13_pPJS

# Plots
par(mfrow = c(2, 1))
plot(x13_pPJS)
par(mfrow = c(1, 1))

# Residuals

par(mfrow = c(3, 2))
plot(x13_pPJS$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_pPJS$decomposition)

## Capture (p) - N-mixture ------------------------------------------------------------

# Cria objeto de série temporal (ts)
ts_pmix <- ts(p.df$pmix, start = c(2018, 2), freq = 12)

# Seasonal adjustment with X-13ARIMA-SEATS
x13_pmix <- RJDemetra::x13(ts_pmix)
x13_pmix

# Plots
par(mfrow = c(2, 1))
plot(x13_pmix)
par(mfrow = c(1, 1))

# Residuals

par(mfrow = c(3, 2))
plot(x13_pmix$regarima)
par(mfrow = c(1, 1))

# Plot S-I ratio
plot(x13_pmix$decomposition)


# Other plots -------------------------------------------------------------

ggplot(
  Toreadicus.captures.env,
  aes(x = trap, y = Toreadicus_perf)
) +
  ggdist::stat_halfeye(
    adjust = 5,
    width = .5,
    .width = .9,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .05,
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Trap", y = "Locomotor performance")


ggplot(
  Toreadicus.captures.env,
  aes(x = trap, y = freq)
) +
  ggdist::stat_halfeye(
    adjust = 5,
    width = .5,
    .width = .9,
    justification = -.2,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .05,
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l",
    ## control range of jitter
    range_scale = .5,
    ## add some transparency
    alpha = .3
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Trap", y = "Abundance")

plot(
  jitter(freq) ~ precip,
  pch = 19,
  col = rgb(0, 0, 0, alpha = .2),
  data = Toreadicus.captures.env,
  bty = "n",
  xlab = "Precipitation (m)",
  ylab = "Abundance",
  xlim = c(0, 200),
  ylim = c(0, 14)
)

plot(
  freq ~ Toreadicus_perf,
  pch = 19,
  col = rgb(0, 0, 0, alpha = .2),
  data = Toreadicus.captures.env,
  bty = "n",
  xlab = "Locomotor performance",
  ylab = "Abundance",
  xlim = c(0.8, 1.8),
  ylim = c(0, 14)
)

plot(
  Toreadicus_perf ~ Toreadicus_ha90,
  pch = 19,
  col = rgb(0, 0, 0, alpha = .2),
  data = Toreadicus.captures.env,
  bty = "n",
  xlab = "Hours of activity",
  ylab = "Locomotor performance",
  ylim = c(0.8, 1.8)
)

plot(
  Toreadicus_perf ~ tmed,
  data = Toreadicus.captures.env,
  pch = 19,
  col = rgb(0, 0, 0, alpha = .2),
  bty = "n",
  xlab = "Mean temperature (ºC)",
  ylab = "Locomotor performance",
  ylim = c(0.8, 1.8),
  xlim = c(20, 28)
)

plot(
  Toreadicus_perf ~ rhmax,
  data = Toreadicus.captures.env,
  pch = 19,
  col = rgb(0, 0, 0, alpha = .2),
  bty = "n",
  xlab = "Relative humidity (%)",
  ylab = "Locomotor performance",
  ylim = c(0.8, 1.8),
  xlim = c(65, 100)
)


# Map location ------------------------------------------------------------
# Load packages
library(tidyverse)
library(readxl)
library(geobr)
library(grid)
library(gt)
library(sp)
library(sf)
library(maptiles)
library(tidyterra)
library(spData)

# Read and plot points coordinates
traps.pts <- read.table("Points_Traps.txt", h = T)
traps.pts

(traps.pts.PEL <- subset(traps.pts, local == "PEL"))
coordinates(traps.pts.PEL) <- c("long", "lat")
proj4string(traps.pts.PEL) <- CRS("+proj=longlat +datum=WGS84")
traps.pts.PEL
traps.pts.PEL <- st_as_sf(traps.pts.PEL)

# Read biomes
biomes <- read_biomes()
biomes

Cerrado <- biomes[biomes$name_biome == "Cerrado", ]

# Read Conservation Units
CUs <- read_conservation_units()

PEL <- CUs[CUs$name_conservation_unit == "PARQUE ESTADUAL DO LAJEADO", ]
APASL <- CUs[
  CUs$name_conservation_unit == "ÁREA DE PROTEÇÃO AMBIENTAL SERRA DO LAJEADO",
]

# Read states
states <- read_state()

# Read South America
world <- data(world)
americas <- filter(world, region_un == "Americas")

ggplot() +
  geom_sf(data = biomes, color = NA, aes(fill = name_biome)) +
  geom_sf(data = states, color = "black", fill = NA) +
  scale_fill_manual(values = viridis::turbo(7)) +
  theme_void()

ggplot() +
  geom_sf(data = Cerrado, color = NA) +
  geom_sf(data = states, color = "black", fill = NA) +
  theme_void()

inset.cerrado <- ggplot() +
  geom_sf(data = americas) +
  geom_sf(data = Cerrado, color = NA, fill = "purple", alpha = 0.5) +
  geom_sf(data = states, fill = NA) +
  geom_sf(data = PEL, color = "forestgreen", fill = NA) +
  geom_rect(
    aes(xmin = -49, xmax = -47.5, ymin = -11, ymax = -9),
    color = "red",
    fill = NA
  ) +
  theme_test() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.background = element_rect(fill = "lightblue")
  ) +
  coord_sf(xlim = c(-74, -35), ylim = c(-34, 6), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

inset.pel <- ggplot() +
  geom_sf(data = PEL, color = "forestgreen", fill = NA, linewidth = 1.25) +
  geom_sf(data = traps.pts.PEL, size = 0.8, color = "red") +
  geom_rect(
    aes(xmin = -48.215, xmax = -48.200, ymin = -10.152, ymax = -10.140),
    color = "red",
    fill = NA
  ) +
  theme_test() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    panel.background = element_rect(fill = "lightgray")
  ) +
  coord_sf(xlim = c(-48.29, -48.15), ylim = c(-10.21, -10), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude")

# dowload tiles and compose raster (SpatRaster)
bbox.arms <- st_bbox(traps.pts.PEL)
bbox.arms[1:4] <- c(-48.215, -10.150, -48.200, -10.140)

imagery <- get_tiles(
  bbox.arms,
  crop = TRUE,
  provider = "Esri.WorldImagery",
  zoom = 16
)
plot(imagery)
plot(st_geometry(traps.pts.PEL), col = "red", pch = 21, add = T)

map.arms.pel <- ggplot() +
  geom_spatraster_rgb(data = imagery, interpolate = F) +
  geom_sf(data = traps.pts.PEL, size = 2.5, color = "red") +
  geom_sf_text(
    data = traps.pts.PEL,
    aes(label = trap),
    nudge_x = 0.0003,
    nudge_y = 0.0003,
    size = 2.5,
    color = "red"
  ) +
  xlim(-48.212, -48.202) +
  ylim(-10.150, -10.142) +
  theme_test() +
  guides(size = "none") +
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("grey60", "white"),
    pad_y = unit(0.2, "in")
  ) +
  ggspatial::annotation_north_arrow(
    location = "bl",
    which_north = "true",
    pad_x = unit(0.53, "in"),
    pad_y = unit(0.5, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) +
  coord_sf(
    xlim = c(-48.212, -48.203),
    ylim = c(-10.150, -10.143),
    expand = FALSE
  ) +
  labs(x = "Longitude", y = "Latitude")


map.arms.pel

# Combining both maps
print(inset.cerrado, vp = viewport(0.896, 0.773, width = 0.2, height = 0.2))
print(inset.pel, vp = viewport(0.737, 0.773, width = 0.2, height = 0.2))
