## load packages
library(xtable)
library(tidytext)
library(DT)
library(psych)
library(ggrepel)
library(factoextra)
library(tm)
library(reshape2)
library(scales)
library(kableExtra)
library(tidyverse)
library(ggpattern)
library(text2vec)
library(data.table)
library(rsample)
library(data.table)
library(magrittr)
library(glmnet)
library(dplyr)
library(ggplot2)
library(labelled)

# seed for simulations
set.seed(1389)

# load data
cand_songs <- read_csv("data/cand_songs_bycampaign.csv") |> 
  mutate(cand_year = fct_reorder(cand_year, year)) |>
  mutate(across(starts_with("NRC_"), ~ 100 * .x / WC), #deal with NRC wordcount issues for scaling
         i_less_we = i-we, fem_male = female-male)  

#### label prep #####
front = c( "Clinton 2016", "Biden 2020",
           "Harris 2024", "Trump 2016", "Trump 2020", "Trump 2024")

emot_colors <- c(
  "gray9", "gray25", "gray40", "gray60",  "#3f1f4d", "#6b3a7d", "#9c6fb0", "#c1a0d4") 

emot_levels <- c("NRC_disgust", "NRC_anger", "NRC_fear", "NRC_sadness", "NRC_trust",
                 "NRC_anticipation", "NRC_surprise", "NRC_joy")
emot_labels <- c("disgust", "anger", "fear", "sadness", "trust", "anticipation", "surprise", "joy")

## graph prep ##
theme_set(theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))

# Set a global ggplot2 theme for all plots
theme_set(
  theme_minimal(base_size = 15)) #+  # base font size


# rename vars for pretty outputs:
var_labels <- c(
  "relig"        = "Religion",
  "politic"      = "Politics",
  "focuspast"    = "Past Focus",
  "focuspresent" = "Present Focus",
  "focusfuture"  = "Future Focus",
  "Drives"       = "Drives",
  "affiliation"  = "Affiliation",
  "Authentic"    = "Authenticity",
  "Clout"        = "Clout",
  "male"         = "Male Terms",
  "female"       = "Female Terms",
  "conflict"     = "Conflict",
  "NRC_trust"    = "Trust (NRC)",
  "NRC_joy"      = "Joy (NRC)",
  "avg_dominance"= "Avg Dominance",
  "i_less_we"    = "I use less We use",
  "fem_male"     = "Female refs less Male refs"
)

for (var in names(var_labels)) {
  if (var %in% names(cand_songs)) {
    var_label(cand_songs[[var]]) <- var_labels[[var]]
  }
}

source("liwc_functions.R")
    
########################### ANALYSIS ########

###### table 1: data information #####
# combine
summary_block <- function(d) {
  d |> summarise(
    "Distinct Candidates" = n_distinct(cand_full),  
    "Total songs" = n(), 
    "Distinct songs" =  n_distinct(title, performer)
  )
}

summary_table <- cand_songs |> group_by(campaign) |> summary_block()
summary_row   <- cand_songs |> summary_block() |> mutate(campaign = "Total (distinct)")

bind_rows(summary_table, relocate(summary_row, campaign)) |>
  xtable(digits = 0) |>
  print(file = "img/summary_camp.tex", type = "latex", include.rownames = FALSE, floating = FALSE)


### summary
summary_table <- cand_songs |> group_by(cand_year) |> summary_block()
summary_row   <- cand_songs |> summary_block() |> mutate(cand_year = "Total (distinct)")

bind_rows(summary_table, relocate(summary_row, cand_year)) |>
  xtable(digits = 0) |>
  print(file = "img/summary_camp_cy.tex", type = "latex", include.rownames = FALSE, floating = FALSE)


########################### SUMMARY STATS: LIWC & NRC ###########################

# 1. Gather all LIWC and NRC variables referenced in analysis
liwc_nrc_vars <- c( "relig", "politic", "focuspast", "focuspresent",
    "focusfuture", "Drives", "affiliation", "Authentic", "Clout", "male", "female",
    "conflict", "NRC_trust", "NRC_joy", "avg_dominance", "i_less_we", "fem_male")


# 2. Extract stats using psych::describe() and select relevant columns
summary_liwc_nrc <- cand_songs |> 
  select(all_of(liwc_nrc_vars)) |> 
  psych::describe() |> 
  as.data.frame() |> 
  rownames_to_column(var = "Variable") |> 
  mutate(Variable = coalesce(var_labels[Variable], Variable)) |>
  select(Variable, n, mean, sd, median, min, max)

# Print preview to console
print(summary_liwc_nrc)

# 3. Export clean table without floating environment (\begin{table})
summary_liwc_nrc |> 
  xtable(
    digits = c(0, 0, 0, 2, 2, 2, 2, 2)  # Formatting decimal places per column
  ) |> 
  print(file = "img/summary_liwc_nrc.tex", type = "latex",  include.rownames = FALSE, floating = FALSE )

###### genre ########
## ---- 1. Contingency table + adjusted standardized residuals ---------------
gp <- cand_songs |>
  filter(party %in% c("D", "R"), !is.na(genre)) |>
  mutate(genre = fct_lump_min(factor(genre), min = 10))   # chi-sq needs expected >= 5

tab <- table(gp$genre, gp$party)
cs  <- chisq.test(tab)

lean <- as.data.frame.table(cs$stdres, responseName = "stdres") |>
  rename(genre = Var1, party = Var2) |>
  filter(party == "D") |>
  mutate(
    # Ensure factors/character types behave correctly
    genre = as.character(genre),
    
    # Compute 2-tailed p-values from stdres (Z-scores)
    p_val = 2 * pnorm(-abs(stdres)),
    
    # Map p-values to significance stars
    stars = case_when(
      p_val < 0.001 ~ "***", p_val < 0.01  ~ "**",  p_val < 0.05  ~ "*",TRUE          ~ ""),
    genre = fct_reorder(genre, stdres), # Sort genres by residual value
    side  = if_else(stdres > 0, "D", "R"), sig   = p_val < 0.05,
    star_x     = stdres + if_else(stdres >= 0, 0.2, -0.2), # Offset star labels 
    star_hjust = if_else(stdres >= 0, 0, 1),
    
    # Add inner bar label text and placement (fixes missing bar_label / label_y error)
    bar_label  = paste0(round(stdres, 2)),
    label_x    = stdres / 2  # Places the number in the middle of the bar
  )

ggplot(lean, aes(x = stdres, y = genre, fill = side)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey80") +
  geom_col(aes(alpha = sig), width = 0.7) + 
  geom_text( # Significance stars at the tip of each bar
    aes(x = star_x, label = stars, hjust = star_hjust),
    vjust = 0.75, size = 4.5, fontface = "bold", colour = "grey20" ) +
  geom_text(  # Inner bar labels showing exact residual values
    aes(x = label_x, label = bar_label), colour = "white",
    size = 3.5, fontface = "bold", vjust = 0.5) +
  scale_fill_manual(values = c(D = "#2a5d9c", R = "#a83232"), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) + # Give room for star labels
  labs(  x = "Adjusted standardized residual", y = NULL,
    title = "Which genres each party over-uses",
    subtitle = "Blue = over-represented among Democrats, red = among Republicans",
    caption = sprintf(
      "Chi-square test of independence: X\u00b2 = %.1f, df = %d, p = %s. Dashed lines mark |residual| = 1.96 (p < .05).\n* p < .05   ** p < .01   *** p < .001 (based on adjusted standardized residuals)",
      cs$statistic, cs$parameter, format.pval(cs$p.value, digits = 3) ) ) +
  theme_minimal(base_size = 12) +
  theme( panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, colour = "grey30", size = 8.5))

ggsave("genre_party_lean.png", path = "img/", width = 7, height = 5.5, dpi = 300)

################# LIWC #########################

###### LIWC: collapsed ####
## 
## CALLS
## 
## H1a: individual vs collective
plot_party(cand_songs, "i_less_we", title = "Individual vs collective language, by party", filename = "liwc_party_iwe.png")
plot_cand_single(cand_songs, "i_less_we", title_prefix = "Individual vs collective language, by campaign", filename_prefix = "liwc_cand_iwe.png")
plot_cand_single(cand_songs, "i_less_we", only_front = TRUE, title_prefix = "Individual vs collective language (Front Runners)", filename_prefix = "liwc_cand_iwe_fr.png")

## H1a: religion
plot_party(cand_songs, c("relig"), title = "Religious & political references, by party", filename = "liwc_party_relig.png", w = 5)
plot_party(cand_songs, c("politic"), title = "Religious & political references, by party", filename = "liwc_party_politic.png", w = 5)
plot_cand_single(cand_songs, c("relig", "politic"), title_prefix = "Religious & political references, by campaign", filename_prefix = "liwc_cand_relig.png", w = 6)

## Gender
plot_party(cand_songs, c("male", "female"), title = "Gender references, by party", filename = "liwc_party_gender.png")
plot_cand_single(cand_songs, c("male", "female", "fem_male"), title_prefix = "Gender references, by campaign", filename_prefix = "liwc_cand_gender.png")

## H1c: temporal orientation
plot_party(cand_songs, c("focuspast", "focuspresent", "focusfuture"), title = "Temporal focus, by party", filename = "liwc_party_focus.png")
plot_cand_single(cand_songs, c("focuspast", "focuspresent", "focusfuture"), title_prefix = "Temporal focus, by campaign", filename_prefix = "liwc_cand_focus.png")


## Cand-specific hypotheses

# 2a: Clinton: Affiliation, Drives
plot_party(cand_songs, c("affiliation", "Drives"), title = "Clinton: Care and Wonk outcomes, by party", filename = "liwc_party_clinton.png")
plot_cand_single(cand_songs, c("affiliation", "Drives"), title_prefix = "Clinton: Care and Wonk outcomes, by campaign", filename_prefix = "liwc_cand_clinton.png")

# 2b: Biden: trust, Authentic
plot_party(cand_songs, c("NRC_trust", "Authentic"), title = "Biden: Trust and Authenticity, by party", filename = "liwc_party_biden.png")
plot_cand_single(cand_songs, c("NRC_trust", "Authentic"), title_prefix = "Biden: Trust and Authenticity, by campaign", filename_prefix = "liwc_cand_biden.png")

# 3: Trump: fem\_male and Clout
plot_party(cand_songs, c("fem_male", "avg_dominance"), title = "Trump: Gender pronoun and Clout, by party", filename = "liwc_party_trump.png")
plot_cand_single(cand_songs, c("fem_male", "avg_dominance"), title_prefix = "Trump: Gender pronoun and Clout, by campaign", filename_prefix = "liwc_cand_trump.png")

## export all these to a table (party then cand):
# Generate the summary table (default unit = "song")
party_stats_df <- summary_party_table(cand_songs, liwc_nrc_vars)

# Quick preview of key comparative columns
party_stats_df |> 
  select(variable, dem_mean, rep_mean, diff_rep_dem, cohen_d, p_adj, sig_adj) |>
  mutate(variable = coalesce(var_labels[variable], variable)) |>
  xtable() |> print("img/tests_table_party.tex", type = "latex", include.rownames = FALSE, floating = FALSE )

# Generate the summary data frame CANDS
# Compute summary statistics
cand_stats_df <- summary_cand_table(cand_songs, liwc_nrc_vars)

# Export individual tables by variable (creates img/cand_table_*.tex and img/all_cand_tables.tex)
export_cand_tables_by_var(cand_stats_df, var_labels)

# 3. Export single combined master file (img/tests_table.tex)
cand_stats_df |> 
  mutate(
    variable_lab = if (!is.null(var_labels)) coalesce(var_labels[variable], variable) else variable
  ) |> 
  select(variable_lab, cand_year, cand_mean, corpus_mean, diff, p_adj, sig_adj) |> 
  xtable(
    caption = "Candidate-level statistical test results",
    label = "tab:cand_tests"
  ) |> 
  print(
    file = "img/tests_table.tex", 
    type = "latex", 
    include.rownames = FALSE, floating = FALSE
  )

###################### METHODS: PCA ##########
keep_vars <- c( liwc_nrc_vars, "id") 

num_data <- cand_songs |> 
  select(where(is.numeric)) |> select(-song_ct) |>
  select(all_of(keep_vars)) |> filter(if_all(all_of(keep_vars), is.finite)) 

nrow(cand_songs) - nrow(num_data) 

# variable in the set.
pca_res <- prcomp(num_data|> select(-id), center = TRUE, scale. = TRUE)

summary(pca_res)
round(pca_res$rotation[, 1:4], 3)

# Scree, to pick how many components to carry forward
round(pca_res$sdev^2 / sum(pca_res$sdev^2), 3)


# # keep id for later merge
 pca_res_id <- data.frame(num_data$id, pca_res$x)
 names(pca_res_id)[1]<-"id"


# prop var explained and variable loadings
fviz_eig(
  pca_res, 
  addlabels = TRUE, 
  labelsize = 5,          # larger label text
  #linecolor = "white"     # variance line white
 )  
ggsave(path = "img/", "pca_eig.png", width = 8)

ylabdim <- c("Dim 2: Time, future to joyful, male, trust ")
xlabdim <- c("Dim 1: Capability and Support to Authentic" )

fviz_pca_var(
  pca_res, 
  col.var = "contrib", 
  axes = c(1, 2), 
  labelsize = 6,                 # larger variable labels
  repel = TRUE,
  gradient.cols = c("gray75", "black"),
  ylab = ylabdim,    
  xlab = xlabdim) #+
ggsave(path = "img/", "pca_var.png", height = 7, width = 9, dpi = 300)


#descriptive table
num_data |> select(-id) |> describe() |> as.data.frame() |>
  rownames_to_column(var = "Variable") |> 
  mutate(Variable = coalesce(var_labels[Variable], Variable)) |>
  select(Variable, min, median, max, n) |>
  xtable() |> print("img/pca_var.tex", type = "latex", include.rownames = FALSE, floating = FALSE )


## Now: PCA and focus on grouping by candidate. #########
# now filter for regular analysis: add relevant elements back in
pca_out <- cand_songs |> 
  select(id, campaign, cand_full, party, genre, title) |>
  inner_join(as.data.frame(pca_res_id), by = "id") |> group_by(party, title) |> 
  mutate(num_party_plays = n()) |> ungroup() 

### CAND LANDSCAPE ##
cand_avg <- pca_out |>
  group_by(campaign, cand_full, party) |>
  summarise(cand_avg_1 = mean(PC1), cand_avg_2 = mean(PC2), .groups = "drop") |>
  group_by(party) |>
  mutate(party_avg_1 = mean(cand_avg_1), party_avg_2 = mean(cand_avg_2)) |> 
  mutate(cand_year = paste(word(cand_full, -1), str_extract(campaign, "\\d{4}"))) |> ungroup()

# Plot only the candidate group labels #
cand_avg |> 
  ggplot() +
  geom_label_repel(aes(x = cand_avg_1, y = cand_avg_2, 
                       label = cand_year, color = party, fill = party), force = 0.05,      
                   force_pull = 3, size = 4.5) +
  geom_label(aes(x = party_avg_1, y = party_avg_2, label = party, color = party), size = 5) + 
  labs(  y = ylabdim,    
         x = xlabdim,
        title = "Candidate Averages (PCA)") + #facet_grid(vars(party)) +
  scale_color_manual(values = party_colors) +
  scale_fill_manual(values = c(  "R" = "lightpink", "D" = "lightblue" ))  
ggsave(filename = "img/pca_cand.png",, width = 9, height = 7)


## PCA BOOTSTRAPPING
## 
## Are campaign centroids stable, or are they noise at high resolution?
## Resample songs WITHIN each campaign, recompute the centroid, repeat.
## 
B <- 500

pca_out <- pca_out |>  
  mutate(cand_year = paste(word(cand_full, -1), str_extract(campaign, "\\d{4}")))

## observed centroids
cent <- pca_out |>
  group_by(cand_year, party) |>
  summarise(n_songs = n(), c1 = mean(PC1), c2 = mean(PC2), .groups = "drop")

## bootstrap
boot <- map_dfr(seq_len(B), function(b) {
  pca_out |>
    group_by(cand_year) |>
    slice_sample(prop = 1, replace = TRUE) |>
    summarise(m1 = mean(PC1), m2 = mean(PC2), .groups = "drop")
})

stab <- boot |>
  group_by(cand_year) |>
  summarise(se1 = sd(m1), se2 = sd(m2),
            lo1 = quantile(m1, .025), hi1 = quantile(m1, .975),
            lo2 = quantile(m2, .025), hi2 = quantile(m2, .975),
            .groups = "drop") |>
  left_join(cent, by = "cand_year")

## 
## The decision numbers
## 
between1 <- sd(cent$c1)
between2 <- sd(cent$c2)

stab  |>
  mutate(ratio1 = se1 / between1, ratio2 = se2 / between2) |>
  select(cand_year, party, n_songs, se1, se2, ratio1, ratio2) |>
  mutate(across(where(is.numeric), ~ round(.x, 3))) |>
  arrange(desc(ratio1)) |> print(n = Inf)

cat("\nbetween-campaign SD:  PC1 =", round(between1, 3),
    " PC2 =", round(between2, 3), "\n")
cat("median noise/signal:  PC1 =", round(median(stab$se1) / between1, 2),
    " PC2 =", round(median(stab$se2) / between2, 2), "\n")

## Independent check: how much of song-level variance does campaign explain?
for (k in c("PC1", "PC2")) {
  a  <- anova(lm(reformulate("cand_year", k), data = pca_out))
  cat(k, " eta^2 =", round(a$`Sum Sq`[1] / sum(a$`Sum Sq`), 3),
      "  p =", format.pval(a$`Pr(>F)`[1], digits = 3), "\n")
}


keep_c <- stab |> pull(cand_year)
dat    <- pca_out |> filter(cand_year %in% keep_c)
camps  <- dat |> distinct(cand_year, party)
by_camp <- split(dat, dat$cand_year)

cat("campaigns kept:", nrow(camps), "\n")
#print(count(camps, party))

## observed party centroid = mean of campaign means
party_cent <- dat |>
  group_by(party, cand_year) |>
  summarise(m1 = mean(PC1), m2 = mean(PC2), .groups = "drop_last") |>
  summarise(n_camp = n(), c1 = mean(m1), c2 = mean(m2), .groups = "drop")

## nested bootstrap: campaigns within party, then songs within campaign
party_boot <- map_dfr(seq_len(B), function(b) {
  map_dfr(unique(camps$party), function(pp) {
    cs <- sample(camps$cand_year[camps$party == pp], replace = TRUE)
    mm <- vapply(cs, function(cc) {
      s <- by_camp[[cc]]
      i <- sample.int(nrow(s), nrow(s), replace = TRUE)
      c(mean(s$PC1[i]), mean(s$PC2[i]))
    }, numeric(2))
    tibble(party = pp, p1 = mean(mm[1, ]), p2 = mean(mm[2, ]))
  })
})

party_stab <- party_boot |>
  group_by(party) |>
  summarise(lo1 = quantile(p1, .025), hi1 = quantile(p1, .975),
            lo2 = quantile(p2, .025), hi2 = quantile(p2, .975),
            .groups = "drop") |>
  left_join(party_cent, by = "party") |>
  mutate(lab = paste0(party, " (", n_camp, " campaigns)"))

print(party_stab)

var_explained <- function(pc, data, group = "cand_year") {
  a  <- anova(lm(reformulate(group, pc), data = data))
  ss <- a$`Sum Sq`; ms_e <- a$`Mean Sq`[2]
  tibble(
    pc    = pc,
    eta2  = ss[1] / sum(ss),
    omega2 = (ss[1] - a$Df[1] * ms_e) / (sum(ss) + ms_e),
    p     = a$`Pr(>F)`[1] )
}

ve <- map_dfr(c("PC1", "PC2"), var_explained, data = dat, group = "cand_year") # dat == plotted data

pct <- function(x) paste0(formatC(100 * x, format = "f", digits = 1), "%")

cap <- paste0(
  "Campaigns with more than five songs. Small points are campaign means ",
  "(crossbars: 95% bootstrap CI over songs).\n",
  "Diamonds are party means of campaign means (crossbars: nested bootstrap over ",
  "campaigns and songs).\n",
  "Candidate accounts for ", pct(ve$eta2[1]), " and ", pct(ve$eta2[2]),
  " of song-level variance in PC1 and PC2 (p = ",
  format.pval(ve$p[1], digits = 2), ", p = ", format.pval(ve$p[2], digits = 2),
  ").")

## making the boostrap plot ####
ggplot(stab, aes(c1, c2, colour = party)) +
  geom_vline(xintercept = 0, colour = "grey88") +
  geom_hline(yintercept = 0, colour = "grey88") +
  
  # campaigns: thin crossbars, small points
  geom_segment(aes(x = lo1, xend = hi1, y = c2, yend = c2), alpha = 0.30) +
  geom_segment(aes(x = c1, xend = c1, y = lo2, yend = hi2), alpha = 0.30) +
  geom_point(size = 1.9) +
  ggrepel::geom_text_repel(aes(label = cand_year), size = 4.5, seed = 13,
                           min.segment.length = 0.2, box.padding = 0.35) +
  # parties: heavy crossbars, diamond
  geom_segment(data = party_stab, linewidth = 1.1,
               aes(x = lo1, xend = hi1, y = c2, yend = c2)) +
  geom_segment(data = party_stab, linewidth = 1.1,
               aes(x = c1, xend = c1, y = lo2, yend = hi2)) +
  geom_point(data = party_stab, aes(c1, c2), shape = 18, size = 6) +
  geom_text(data = party_stab, aes(c1, c2, label = lab),
            vjust = -1.4, fontface = "bold", size = 5.5, show.legend = FALSE) +
  
  scale_colour_manual(values = party_colors, name = "Party") +
  labs(y = ylabdim, x = xlabdim,
       title = "Campaign and party centroids with bootstrap uncertainty",
       caption = cap) + coord_cartesian(ylim = c(-1.1, 1.7), xlim = c(-1.5, 1.55)) +
  theme(plot.caption = element_text(hjust = 0, colour = "grey30", size = 10), legend.position = "bottom")

ggsave(path = "img/", "pca_cand_boot.png", width = 12, height = 10, dpi = 300)

## does the campaign-level t-test agree?
cm <- pca_out |> group_by(cand_year, party) |>
  summarise(m = mean(PC2), .groups = "drop")
t.test(m ~ party, data = cm)

stab |>
print(file = "img/summary_pca_stab.tex", type = "latex",  include.rownames = FALSE, floating = FALSE )


# ####### testing PCA vs genre
library(lmerTest)
stab <- stab |>
  mutate(cand = sub("\\d{4}$", "", cand_year))

# baseline: campaign PC2 by party, Trump's campaigns sharing an intercept
summary(lmer(c2 ~ party + (1 | cand), data = stab))

# song level, genre controlled (lump rare genres so the model is estimable)
song_g <- dat |> mutate(genre = fct_lump_min(factor(genre), min = 10))

m_base  <- lmer(PC2 ~ party         + (1 | cand_year), data = song_g)
m_genre <- lmer(PC2 ~ party + genre + (1 | cand_year), data = song_g)

summary(m_base)
summary(m_genre)

# does party's effect survive once genre is in the model, and does genre
# improve fit at all?
anova(m_base, m_genre)

count(song_g, genre)          # confirm reference level and what fell into "Other"
levels(song_g$genre)[1]

## fixed-effects comparison: does the party effect survive adding genre? ##
fixed_tbl <- function(m, suffix) {
  co <- summary(m)$coefficients
  tibble(
    term = rownames(co),
    est  = co[, "Estimate"],
    se   = co[, "Std. Error"],
    p    = co[, "Pr(>|t|)"]
  ) |>
    rename_with(~ paste0(.x, "_", suffix), c(est, se, p))
}

genre_table <- full_join(
  fixed_tbl(m_base,  "base"),
  fixed_tbl(m_genre, "genre"),
  by = "term"
) |>
  mutate(term = recode(term,
                       "(Intercept)"   = "Intercept",
                       "partyR"        = "Party (R)",
                       "genrepop/rock" = "Genre: Pop/Rock",
                       "genrer&b"      = "Genre: R&B",
                       "genreOther"    = "Genre: Other"
  )) |>
  rename(
    Term              = term,
    `Est. (no genre)` = est_base,
    `SE (no genre)`   = se_base,
    `p (no genre)`    = p_base,
    `Est. (+ genre)`  = est_genre,
    `SE (+ genre)`    = se_genre,
    `p (+ genre)`     = p_genre
  )

genre_table |>
  xtable(digits = c(0, 0, 3, 3, 3, 3, 3, 3)) |>
  print(file = "img/summary_pca_genre.tex", type = "latex",
        include.rownames = FALSE, floating = FALSE)

## model comparison: does genre improve fit at all? ##
comp_tbl <- anova(m_base, m_genre) |>
  as.data.frame() |>
  mutate(Model = c("Party only", "Party + genre"), .before = 1) |>
  select(Model, npar, AIC, BIC, Chisq, Df, p = `Pr(>Chisq)`)

comp_tbl |>
  xtable(digits = c(0, 0, 0, 1, 1, 2, 0, 3)) |>
  print(file = "img/summary_pca_genre_lrt.tex", type = "latex",
        include.rownames = FALSE, floating = FALSE)



# ## leave one campaign out
loo_table <- function(pc) {
  map_dfr(unique(pca_out$cand_year), function(cc) {
    f <- reformulate(c("party", "(1 | cand_year)"), response = pc)
    m <- lmer(f, data = filter(pca_out, cand_year != cc))
    s <- summary(m)$coefficients["partyR", ]
    tibble(cand_year = cc, est = s[["Estimate"]], p = s[["Pr(>|t|)"]])
  })
}

loo_pc1 <- loo_table("PC1") |> rename(pc1_est = est, pc1_p = p)
loo_pc2 <- loo_table("PC2") |> rename(pc2_est = est, pc2_p = p)

loo_combined <- full_join(loo_pc1, loo_pc2, by = "cand_year") |>
  arrange(pc2_p) |>
  rename(
    `Campaign Dropped` = cand_year,
    `Dim 1 Est.`  = pc1_est,
    `Dim 1 p`     = pc1_p,
    `Dim 2 Est.`  = pc2_est,
    `Dim 2 p`     = pc2_p )

loo_combined |>
  xtable(digits = c(0, 0, 3, 3, 3, 3)) |>
  print(file = "img/summary_loo_pc.tex", type = "latex",
        include.rownames = FALSE, floating = FALSE)



cand_songs |> select(-starts_with("lyric")) |> write_csv("data/cand_songs_replication.csv") 
