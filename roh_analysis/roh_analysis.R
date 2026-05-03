# =========================================================
# ROH Analysis & Inbreeding (FROH)
# POP1 vs POP2 comparison
# =========================================================

library(detectRUNS)
library(dplyr)
library(tidyr)
library(ggplot2)

# =========================================================
# 1. INPUT FILES
# =========================================================
pop1_ped <- "POP1.ped"
pop1_map <- "POP1.map"

pop2_ped <- "POP2.ped"
pop2_map <- "POP2.map"

# =========================================================
# 2. ROH DETECTION FUNCTION
# =========================================================
run_roh_analysis <- function(ped, map) {

  roh_1_2 <- consecutiveRUNS.run(
    genotypeFile = ped,
    mapFile = map,
    minSNP = 15,
    ROHet = FALSE,
    minLengthBps = 1e6,
    maxGap = 1e6,
    maxOppRun = 0,
    maxMissRun = 0
  )

  roh_2_4 <- consecutiveRUNS.run(ped, map, minSNP = 15, minLengthBps = 2e6)
  roh_4_8 <- consecutiveRUNS.run(ped, map, minSNP = 15, minLengthBps = 4e6)
  roh_8p  <- consecutiveRUNS.run(ped, map, minSNP = 15, minLengthBps = 8e6)

  roh_all <- rbind(roh_1_2, roh_2_4, roh_4_8, roh_8p)

  return(list(
    roh_all = roh_all,
    roh_1_2 = roh_1_2,
    roh_2_4 = roh_2_4,
    roh_4_8 = roh_4_8,
    roh_8p  = roh_8p
  ))
}

# =========================================================
# 3. RUN ANALYSIS
# =========================================================
pop1_roh <- run_roh_analysis(pop1_ped, pop1_map)
pop2_roh <- run_roh_analysis(pop2_ped, pop2_map)

# =========================================================
# 4. FROH SUMMARY
# =========================================================
extract_froh <- function(rohs) {
  rohs$roh_all %>%
    group_by(id) %>%
    summarise(FROH = sum(lengthBps, na.rm = TRUE))
}

froh_pop1 <- extract_froh(pop1_roh)
froh_pop2 <- extract_froh(pop2_roh)

froh_pop1$Population <- "POP1"
froh_pop2$Population <- "POP2"

froh_all <- bind_rows(froh_pop1, froh_pop2)

# =========================================================
# 5. ROH LENGTH DISTRIBUTION
# =========================================================
roh_length_summary <- function(rohs, pop_name) {
  rohs$roh_all %>%
    mutate(Population = pop_name) %>%
    mutate(Class = case_when(
      lengthBps < 2e6 ~ "1-2Mb",
      lengthBps < 4e6 ~ "2-4Mb",
      lengthBps < 8e6 ~ "4-8Mb",
      TRUE ~ ">8Mb"
    )) %>%
    count(Population, Class)
}

roh_dist <- bind_rows(
  roh_length_summary(pop1_roh, "POP1"),
  roh_length_summary(pop2_roh, "POP2")
)

# =========================================================
# 6. VISUALIZATION
# =========================================================
ggplot(froh_all, aes(x = Population, y = FROH, fill = Population)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Inbreeding (FROH) comparison",
    y = "FROH"
  )

ggplot(roh_dist, aes(x = Class, y = n, fill = Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "ROH length distribution",
    x = "ROH class",
    y = "Count"
  )

# =========================================================
# 7. EXPORT RESULTS
# =========================================================
write.csv(froh_all, "FROH_results.csv", row.names = FALSE)
write.csv(roh_dist, "ROH_distribution.csv", row.names = FALSE)
