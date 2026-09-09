###############################################################################
## LIWC ANALYSIS & PLOTTING FUNCTIONS
##
## Unit of analysis for the party-level functions is set by unit:
##   "song"      (default) every song counts once; n = number of songs
##               (SEs are adjusted for intra-campaign clustering via n_eff)
##   "campaign"  average within cand_year first; n = number of campaigns
##   "candidate" average within cand_full first; n = number of candidates
##               (collapses Trump 2016/2020/2024 into ONE observation)
###############################################################################

library(tidyverse)

party_colors <- c("D" = "#2a5d9c", "R" = "#a83232")
front <- c("Clinton 2016", "Biden 2020", "Harris 2024",
           "Trump 2016", "Trump 2020", "Trump 2024")

# Significance star helper
star <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    TRUE     ~ ""
  )
}

# Pivot LIWC variables into long format
liwc_long <- function(data, liwc_vars) {
  missing <- setdiff(liwc_vars, names(data))
  if (length(missing)) stop("not in data: ", paste(missing, collapse = ", "))
  data |>
    select(cand_full, party, campaign, cand_year, all_of(liwc_vars)) |>
    pivot_longer(cols = all_of(liwc_vars), names_to = "liwc", values_to = "value") |>
    filter(!is.na(value))
}

# One row per unit of analysis, carrying both id columns in every case
liwc_units <- function(long, unit) {
  switch(unit,
         song = long |>
           select(party, cand_full, cand_year, liwc, value) |>
           mutate(n_songs = 1L),
         campaign = long |>
           group_by(party, cand_year, cand_full, liwc) |>
           summarise(value = mean(value, na.rm = TRUE), n_songs = n(), .groups = "drop"),
         candidate = long |>
           group_by(party, cand_full, liwc) |>
           summarise(value = mean(value, na.rm = TRUE), n_songs = n(), .groups = "drop") |>
           mutate(cand_year = cand_full)
  )
}

# One-way random-effects ICC, design effect, and effective n per measure
liwc_icc <- function(data, liwc_vars, group = "cand_year") {
  liwc_long(data, liwc_vars) |>
    group_by(liwc) |>
    group_modify(function(d, k) {
      d[[group]] <- factor(d[[group]])
      ni <- as.numeric(table(d[[group]])); ni <- ni[ni > 0]
      if (length(ni) < 2 || sum(ni) <= length(ni)) {
        return(tibble(n_songs = sum(ni), n_groups = length(ni),
                      icc = NA_real_, design_effect = NA_real_, n_eff = NA_real_))
      }
      a   <- summary(aov(reformulate(group, "value"), data = d))[[1]]
      msb <- a[["Mean Sq"]][1]; msw <- a[["Mean Sq"]][2]
      k0  <- (sum(ni) - sum(ni^2) / sum(ni)) / (length(ni) - 1)
      icc <- max(0, (msb - msw) / (msb + (k0 - 1) * msw))
      deff <- 1 + (k0 - 1) * icc
      tibble(n_songs = sum(ni), n_groups = length(ni),
             icc = icc, design_effect = deff, n_eff = sum(ni) / deff)
    }) |>
    ungroup()
}

# Aggregated party-level means and CIs (UPDATED WITH OPTION 1: n_eff Design Effect Adjustments)
liwc_by_party <- function(data, liwc_vars, conf = 0.95,
                          unit = c("song", "campaign", "candidate")) {
  unit <- match.arg(unit)
  long <- liwc_long(data, liwc_vars)
  
  # Composite/signed measures check
  neg_flag <- long |>
    group_by(liwc) |>
    summarise(has_neg = any(value < 0, na.rm = TRUE), .groups = "drop")
  
  # Calculate ICC and design effects when running song-level units
  icc_stats <- if (unit == "song") {
    liwc_icc(data, liwc_vars, group = "cand_year") |>
      select(liwc, design_effect)
  } else {
    tibble(liwc = liwc_vars, design_effect = 1)
  }
  
  liwc_units(long, unit) |>
    group_by(party, liwc) |>
    summarise(
      n_camps = n_distinct(cand_year),
      n_cands = n_distinct(cand_full),
      n_songs = sum(n_songs),
      n_unit  = n(),                       # Raw observation count
      mean    = mean(value, na.rm = TRUE),
      sd      = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(neg_flag, by = "liwc") |>
    left_join(icc_stats, by = "liwc") |>
    mutate(
      # Adjust effective sample size (n_eff) using Design Effect if unit == "song"
      n_eff = if_else(unit == "song" & !is.na(design_effect) & design_effect > 0, 
                      n_unit / design_effect, 
                      as.numeric(n_unit)),
      
      # Standard Error adjusted for effective N
      se    = if_else(n_eff > 1, sd / sqrt(n_eff), NA_real_),
      
      # Degrees of freedom adjusted for effective N
      df    = pmax(1, n_eff - 1),
      tcrit = if_else(!is.na(df) & df > 0, qt(1 - (1 - conf) / 2, df), NA_real_),
      
      lo    = if_else(has_neg, mean - tcrit * se, pmax(0, mean - tcrit * se)),
      hi    = mean + tcrit * se,
      x_lab = paste0(party, "\n(", n_songs, " songs, ", n_camps, " campaigns)")
    )
}

# Party difference hypothesis tests (Welch t-test)
liwc_party_test <- function(data, liwc_vars,
                            unit = c("song", "campaign", "candidate"),
                            var.equal = FALSE) {
  unit <- match.arg(unit)
  long <- liwc_long(data, liwc_vars)
  liwc_units(long, unit) |>
    group_by(liwc) |>
    group_modify(function(d, k) {
      R <- d$value[d$party == "R"]; D <- d$value[d$party == "D"]
      nR <- length(R); nD <- length(D)
      if (nR < 2 || nD < 2) {
        return(tibble(n_R = nR, n_D = nD, mean_R = mean(R), mean_D = mean(D),
                      diff = mean(R) - mean(D), cohen_d = NA_real_,
                      t = NA_real_, df = NA_real_, p = NA_real_, sig = ""))
      }
      tt <- t.test(R, D, var.equal = var.equal)
      pooled <- sqrt(((nR - 1) * var(R) + (nD - 1) * var(D)) / (nR + nD - 2))
      tibble(n_R = nR, nD = nD,
             mean_R = mean(R), mean_D = mean(D),
             diff = mean(R) - mean(D),
             cohen_d = (mean(R) - mean(D)) / pooled,
             t = unname(tt$statistic), df = unname(tt$parameter),
             p = tt$p.value, sig = star(tt$p.value))
    }) |>
    ungroup() |>
    mutate(p_adj = p.adjust(p, method = "BH"), sig_adj = star(p_adj))
}

# Candidate vs rest-of-corpus, two-sample Welch, BH-adjusted within measure
liwc_by_cand <- function(data, liwc_vars, conf = 0.95) {
  long <- liwc_long(data, liwc_vars)
  long |>
    group_by(liwc) |>
    group_modify(function(dl, kl) {
      keys <- dl |>
        group_by(cand_year) |>
        summarise(cand_full = first(cand_full),
                  party     = first(party),
                  campaign  = first(campaign), .groups = "drop")
      map_dfr(seq_len(nrow(keys)), function(i) {
        k    <- keys[i, ]
        idx  <- dl$cand_year == k$cand_year
        own  <- dl$value[idx]
        rest <- dl$value[!idx]
        n_songs <- length(own)
        m  <- mean(own)
        s  <- sd(own)
        se <- if (n_songs > 1) s / sqrt(n_songs) else NA_real_
        tc <- if (n_songs > 1) qt(1 - (1 - conf) / 2, n_songs - 1) else NA_real_
        ok <- n_songs > 1 && length(rest) > 1 && !is.na(s) && s > 0
        tt <- if (ok) t.test(own, rest) else NULL
        bind_cols(k, tibble(
          n_songs = n_songs,
          mean    = m,
          sd      = s,
          ref     = mean(rest),
          ref_all = mean(dl$value),
          se      = se,
          lo      = m - tc * se,
          hi      = m + tc * se,
          diff    = m - mean(rest),
          t       = if (ok) unname(tt$statistic) else NA_real_,
          df      = if (ok) unname(tt$parameter) else NA_real_,
          p       = if (ok) tt$p.value else NA_real_
        ))
      })
    }) |>
    ungroup() |>
    group_by(liwc) |>
    mutate(
      p_adj   = p.adjust(p, method = "BH"),
      sig     = star(p),
      sig_adj = star(p_adj)
    ) |>
    ungroup()
}

# Plot party comparisons
# Plot party comparisons with in-bar estimates and SEs
plot_party <- function(data, liwc_vars,
                       unit = c("song", "campaign", "candidate"),
                       title = "LIWC by party", filename = NULL,
                       use_adj = TRUE, w = 7, h = 5) {
  unit <- match.arg(unit)
  d   <- liwc_by_party(data, liwc_vars, unit = unit)
  tst <- liwc_party_test(data, liwc_vars, unit = unit)
  star_col <- if (use_adj) "sig_adj" else "sig"
  neg      <- any(d$mean < 0, na.rm = TRUE)
  
  # Format the estimate and SE label (e.g., "3.42\n(0.15)")
  d <- d |>
    mutate(
      bar_label = if_else(
        !is.na(mean) & !is.na(se),
        sprintf("%.2f\n(%.2f)", mean, se),
        sprintf("%.2f", mean)
      ),
      # Position label at 25% height of the bar for center alignment
      label_y = mean * 0.25
    )
  
  ann <- d |>
    group_by(liwc) |>
    summarise(top  = max(hi, na.rm = TRUE),
              span = diff(range(c(lo, hi), na.rm = TRUE)),
              .groups = "drop") |>
    mutate(y = top + 0.10 * if_else(is.finite(span) & span > 0, span, abs(top))) |>
    left_join(tst, by = "liwc") |>
    mutate(sig_use = .data[[star_col]]) |>
    filter(!is.na(sig_use), sig_use != "") |>
    mutate(x = 1.5)
  
  bar_txt <- switch(unit,
                    song      = "Bars are party means across all songs (SEs adjusted via n_eff); whiskers 95% CI.",
                    campaign  = "Bars are party means across campaign means; whiskers 95% CI.",
                    candidate = "Bars are party means across candidate means; whiskers 95% CI.")
  test_txt <- switch(unit,
                     song      = "songs (nested within campaigns; see liwc_icc for effective n)",
                     campaign  = "campaigns",
                     candidate = "candidates")
  
  p <- ggplot(d, aes(x = x_lab, y = mean, fill = party)) +
    (if (neg) geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.6) else NULL) +
    geom_col(width = 0.7, colour = NA) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.9, colour = "grey40") +
    # White text inside the bars (centered at 50% bar height)
    geom_text(
      aes(y = label_y, label = bar_label),
      colour = "white",
      size = 3.5,
      fontface = "bold",
      lineheight = 0.85,
      vjust = 0.5
    ) +
    geom_text(data = ann, aes(x = x, y = y, label = sig_use),
              inherit.aes = FALSE, size = 8, vjust = 0) +
    facet_wrap(~ liwc, scales = "free_y") +
    scale_fill_manual(values = party_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = if (neg) c(0.15, 0.22) else c(0, 0.22))) +
    labs(x = NULL, y = "Mean % of words per song", title = title,
         caption = paste0(
           bar_txt, "\n",
           "* p<.05  ** p<.01  *** p<.001 (Welch t, R vs D across ", test_txt,
           if (use_adj) "; BH-adjusted)" else ")")) +
    theme(panel.grid.major.x = element_blank())
  
  if (!is.null(filename)) {
    ggsave(path = "img/", filename = filename, plot = p, width = w, height = h, dpi = 300)
  }
  print(tst, width = Inf)
  return(p)
}

# Plot individual campaigns against the corpus mean, BH-adjusted stars
plot_cand_single <- function(data, liwc_vars, title_prefix = "LIWC",
                             filename_prefix = NULL, only_front = FALSE,
                             w = 9, h = 7) {
  plots <- list()
  for (var in liwc_vars) {
    d <- liwc_by_cand(data, var)
    if (only_front) d <- d |> filter(cand_year %in% front)
    if (nrow(d) == 0) {
      warning("no rows for ", var, " after filtering; skipped")
      next
    }
    var_title    <- paste0(title_prefix, ": ", var)
    var_filename <- if (!is.null(filename_prefix)) {
      sub("(\\.[a-zA-Z0-9]+)$", paste0("_", var, "\\1"), filename_prefix)
    } else {
      NULL
    }
    overall_mean <- d$ref_all[1]
    rng     <- c(d$lo, d$hi, d$mean)
    pad_amt <- if (all(is.na(rng))) 0 else max(abs(rng), na.rm = TRUE) * 0.08
    d <- d |>
      mutate(
        star_x     = if_else(mean >= 0,
                             pmax(hi, mean, 0, na.rm = TRUE) + pad_amt,
                             pmin(lo, mean, 0, na.rm = TRUE) - pad_amt),
        star_hjust = if_else(mean >= 0, 0, 1)
      )
    p <- ggplot(d, aes(x = reorder(cand_year, mean), y = mean, fill = party)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.6) +
      geom_col(width = 0.75, colour = NA) +
      geom_errorbar(aes(ymin = pmin(lo, hi), ymax = pmax(lo, hi)),
                    width = 0.2, linewidth = 0.8, colour = "grey40") +
      geom_text(aes(y = star_x, label = sig_adj, hjust = star_hjust),
                vjust = 0.75, size = 6) +
      geom_hline(yintercept = overall_mean, linetype = "dashed",
                 color = "grey40", linewidth = 0.8) +
      coord_flip() +
      scale_fill_manual(values = party_colors, name = "Party") +
      scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
      labs(
        x = NULL,
        y = "Mean % of words per song",
        title = var_title,
        caption = paste0(
          "Dashed line = corpus mean across all R & D songs.\n",
          "* FDR p<.05  ** FDR p<.01  *** FDR p<.001 (Welch t, campaign vs rest of ",
          "corpus; Benjamini-Hochberg adjusted)")
      ) +
      theme(panel.grid.major.y = element_blank())
    if (!is.null(var_filename)) {
      ggsave(path = "img/", filename = var_filename, plot = p, width = w, height = h, dpi = 300)
    }
    plots[[var]] <- p
    print(p)
  }
  invisible(plots)
}