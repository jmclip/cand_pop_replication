###############################################################################
## 2_tech_appendix.R
## Supplemental Appendix E (text models) and F (PCA robustness)
##
## RUN 1_analysis.R FIRST, in the same session. This script uses its objects:
##   cand_songs, liwc_nrc_vars, keep_vars, num_data, pca_res, rot_load,
##   rot_scores, var_labels, party_colors, fmt_p
## It never overwrites them: every object created here is prefixed
##   bow_ / lda_ / rb_ / a2_
##
## Outputs (img/):
##   F: appendix_bow_sensitivity.png, appendix_lda_topic_k.png,
##      appendix_lda_stability_ari.png
##   G: pca_loo_similarity.png, pca_var_partygap.tex, pca_loco.tex,
##      pca_spec_compare.tex, pca_cand_boot_dim2_dim3.png
## Note: Claude Opus 5.5 used to help streamline this code and ensure all outputs from paper included
###############################################################################

## ---- switches ------------------------------------------------------------
RUN_BOW   <- TRUE    # Appendix F: bag-of-words classifier
RUN_LDA   <- TRUE    # Appendix F: LDA diagnostics (slow: searchK, ~10+ min)
RUN_PCA   <- TRUE    # Appendix G: PCA robustness + Figure A2
RUN_EXTRA <- FALSE   # exploratory output NOT reported in the paper

## ---- packages ------------------------------------------------------------
library(text2vec)     # BOW document-term matrix
library(glmnet)       # lasso classifier
library(rsample)      # grouped train/test split
library(pROC)         # AUC confidence interval
library(tm)           # LDA preprocessing
library(SnowballC)    # stemming (used by tm)
library(slam)         # sparse matrix row sums
library(topicmodels)  # LDA
library(stm)          # searchK diagnostics
library(mclust)       # adjusted Rand index
library(psych)        # Promax 
library(xtable)
library(tidyverse)    # load last so dplyr functions are not masked

need <- c("cand_songs", "liwc_nrc_vars", "keep_vars", "num_data", "pca_res",
          "rot_load", "rot_scores", "var_labels", "party_colors", "fmt_p")
miss <- need[!vapply(need, exists, logical(1))]
if (length(miss)) stop("Run 1_analysis.R first; missing: ", paste(miss, collapse = ", "))

pct <- function(x) paste0(formatC(100 * x, format = "f", digits = 1), "%")


###############################################################################
## Appendix F.1  BAG-OF-WORDS CLASSIFIER
###############################################################################
if (RUN_BOW) {
  
  NFOLDS <- 5
  
  # Same songs as the paper (cand_songs from script 1); one row per song x party
  bow_data <- cand_songs |>
    mutate(party01 = if_else(party == "D", 1, 0)) |>
    filter(!is.na(genre), nchar(lyrics_liwc) > 50) |>
    distinct(title, performer, party, .keep_all = TRUE) |>
    mutate(doc_id = as.character(row_number()))
  
  bow_iter <- function(d) itoken(d$lyrics_liwc, preprocessor = tolower,
                                 tokenizer = word_tokenizer, ids = d$doc_id,
                                 progressbar = FALSE)
  
  # Fold ids assigned per title, so all rows of a title share a fold
  bow_fold_by_title <- function(titles, k) {
    u <- unique(titles)
    f <- sample(rep_len(seq_len(k), length(u)))
    f[match(titles, u)]
  }

  # One split, vocabulary built from the TRAINING songs only (no test leakage).
  # Songs are kept together: a title appears in train or test, never both.
  bow_fit_split <- function(seed, data = bow_data, prop = 0.80) {
    set.seed(seed)
    sp    <- group_initial_split(data, group = title, prop = prop)
    train <- training(sp); test <- testing(sp)
    it_tr <- bow_iter(train); it_te <- bow_iter(test)
    vec   <- vocab_vectorizer(create_vocabulary(it_tr))
    dtm_tr <- create_dtm(it_tr, vec); dtm_te <- create_dtm(it_te, vec)
    # internal CV folds also grouped by title: a song used by both parties
    # must not sit in a training fold and its validation fold at once
    foldid <- bow_fold_by_title(train$title, NFOLDS)
    stopifnot(all(table(foldid, train$party01) > 0))           # both parties in every fold
    fit <- cv.glmnet(x = dtm_tr, y = train$party01, family = "binomial",
                     alpha = 1, type.measure = "auc", foldid = foldid,
                     control = list(thresh = 1e-3, maxit = 1e3))
    preds <- predict(fit, dtm_te, type = "response")[, 1]
    list(fit = fit, test = test, preds = preds,
         row = data.frame(seed = seed, cv_auc = max(fit$cvm),
                          test_auc = as.numeric(glmnet:::auc(test$party01, preds))))
  }
  
  ## Main split
  bow_main <- bow_fit_split(1389)
  print(bow_main$row)
  print(pROC::ci.auc(pROC::roc(bow_main$test$party01, bow_main$preds, quiet = TRUE)))
  
  ## Stability across random splits (Figure in Appendix)
  bow_seeds   <- c(1389, 13, 87, 89, 2026, 1015, 216, 410, 2012, 2016,
                   1983, 1130, 12, 16, 83, 84)
  bow_results <- map_dfr(bow_seeds, \(s) bow_fit_split(s)$row)
  
  bow_results |>
    summarise(n_splits = n(),
              cv_mean = mean(cv_auc), cv_sd = sd(cv_auc),
              test_mean = mean(test_auc), test_sd = sd(test_auc),
              test_min = min(test_auc), test_max = max(test_auc)) |>
    print()
  
  bow_results |>
    pivot_longer(c(cv_auc, test_auc), names_to = "metric", values_to = "auc") |>
    mutate(metric = recode(metric, cv_auc = "CV AUC (train)", test_auc = "Test AUC")) |>
    ggplot(aes(metric, auc)) +
    geom_jitter(width = 0.08, height = 0, alpha = .6, size = 2) +
    geom_boxplot(fill = "transparent") +
    stat_summary(fun = mean, geom = "point", size = 3, colour = "black") +
    stat_summary(fun = mean, geom = "text", label = "mean", hjust = -0.5,
                 vjust = 0.5, colour = "black", size = fig_text) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    labs(x = NULL, y = "AUC",
         title = paste("BOW classifier performance across", length(bow_seeds), "random splits"),
         subtitle = "Dashed line = chance (AUC = 0.5)") +
    theme_minimal(base_size = 15) + theme_fig()
  ggsave_fig("img/appendix_bow_sensitivity.png", width = 8, height = 6)
  
  ## ---- exploratory (not reported): stop words + n-grams, top coefficients ----
  if (RUN_EXTRA) {
    set.seed(1389)
    sp    <- group_initial_split(bow_data, group = title, prop = 0.80)
    train <- training(sp); test <- testing(sp)
    stop_words_filtered <- tidytext::stop_words |>
      filter(!word %in% c("good", "kind", "problem", "wanting", "important", "parting")) |>
      pull(word) |> c("ass", "uh", "ooh", "hoo", "ah")
    vocab_ns <- create_vocabulary(bow_iter(train), ngram = c(1L, 3L),
                                  stopwords = stop_words_filtered) |>
      prune_vocabulary(doc_count_min = 4, doc_proportion_max = 0.5)
    vec_ns <- vocab_vectorizer(vocab_ns)
    fit_ns <- cv.glmnet(create_dtm(bow_iter(train), vec_ns), train$party01,
                        family = "binomial", alpha = 1, type.measure = "auc",
                        foldid = bow_fold_by_title(train$title, NFOLDS))
    preds_ns <- predict(fit_ns, create_dtm(bow_iter(test), vec_ns), type = "response")[, 1]
    print(glmnet:::auc(test$party01, preds_ns))
    cm <- coef(fit_ns, s = "lambda.1se")
    data.frame(term = rownames(cm), estimate = as.vector(cm)) |>
      filter(term != "(Intercept)", estimate != 0) |>
      mutate(party = if_else(estimate > 0, "Democrat", "Republican")) |>
      group_by(party) |> slice_max(abs(estimate), n = 10) |> ungroup() |>
      ggplot(aes(reorder(term, estimate), estimate, fill = party)) +
      geom_col(show.legend = FALSE, width = 0.7) + coord_flip() +
      facet_wrap(~ party, scales = "free") +
      scale_fill_manual(values = c(Democrat = "#2b5c8f", Republican = "#a83232")) +
      labs(x = NULL, y = "Coefficient (log-odds, lambda.1se)") + theme_minimal()
  }
}


###############################################################################
## APPENDIX E.2  LDA TOPIC MODELS
###############################################################################
if (RUN_LDA) {
  
  lda_text <- cand_songs |>
    select(performer, title, cand_year, party, lyrics_liwc) |>
    mutate(doc_id = as.character(row_number()), text = as.character(lyrics_liwc))
  
  # Stop words: local copy so the script does not depend on a website.
  # First run downloads and saves it to data/; later runs read the local file.
  lda_stop_file <- "data/stopwords_en.txt"
  if (!file.exists(lda_stop_file)) {
    download.file("https://slcladal.github.io/resources/stopwords_en.txt", lda_stop_file)
  }
  lda_stopwords <- readLines(lda_stop_file, encoding = "UTF-8")
  
  lda_corpus <- VCorpus(DataframeSource(lda_text)) |>
    tm_map(content_transformer(tolower)) |>
    tm_map(removePunctuation, preserve_intra_word_dashes = TRUE) |>
    tm_map(removeNumbers) |>
    tm_map(removeWords, lda_stopwords) |>
    tm_map(stemDocument, language = "en") |>
    tm_map(stripWhitespace)
  
  lda_dtm <- DocumentTermMatrix(lda_corpus, control = list(
    bounds = list(global = c(5, Inf)),               # terms in >= 5 songs
    stopwords = stemDocument(lda_stopwords)))
  
  lda_remove <- c(
    "ooh", "yeah", "nigga", "whoa", "uh-huh", "hey", "woah", "ooh-ooh", "ooh-ooh-ooh",
    "the", "and", "you", "your", "for", "that", "but", "can", "just",
    "what", "all", "this", "from", "her", "when", "there", "one", "was",
    "like", "get", "got", "let", "feel", "know", "time", "man", "will", "with",
    "dont", "wanna", "aint", "gonna", "gotta", "cant", "youre", "youll",
    "ive", "ill", "thats", "its", "hes", "shes", "theyre", "were",
    "cause", "cuz", "bout", "til", "em")
  lda_dtm <- lda_dtm[, !(colnames(lda_dtm) %in% unique(c(lda_remove, stemDocument(lda_remove))))]
  
  lda_keep <- slam::row_sums(lda_dtm) > 0           # LDA cannot take empty documents
  lda_dtm  <- lda_dtm[lda_keep, ]
  lda_text <- lda_text[lda_keep, ]
  dim(lda_dtm)
  
  ## ---- choosing K: stm::searchK diagnostics (Appendix F figure) ----
  ## NOTE for the caption: these diagnostics come from STM (structural topic
  ## model) fits via stm::searchK, used as a guide to K for the LDA models.
  lda_stm_corpus <- readCorpus(lda_dtm, type = "slam")
  lda_searchK <- searchK(documents = lda_stm_corpus$documents,
                         vocab     = lda_stm_corpus$vocab,
                         K         = seq(2, 25, by = 2),
                         seed      = 1989)
  
  lda_searchK$results |>
    select(K, heldout, residual, semcoh, exclus) |>
    unnest(cols = c(K, heldout, residual, semcoh, exclus)) |>
    pivot_longer(-K, names_to = "metric", values_to = "value") |>
    mutate(metric = recode(metric,
                           heldout  = "Held-out likelihood\n(higher = better)",
                           residual = "Residuals\n(lower = better)",
                           semcoh   = "Semantic coherence\n(higher = better)",
                           exclus   = "Exclusivity\n(higher = better)")) |>
    ggplot(aes(K, value)) +
    geom_line(colour = "grey40", linewidth = 0.7) +
    geom_point(colour = "grey20", size = 2.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    labs(x = "Number of Topics (K)", y = NULL,
         title = "Diagnostic Values by Number of Topics (stm::searchK)") +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) +
    theme_fig()
  ggsave_fig("img/appendix_lda_topic_k.png", width = 8, height = 6)
  
  ## ---- seed stability: does the LDA solution replicate across seeds? ----
  lda_seeds   <- c(13, 87, 89, 8789, 2026)
  lda_K_check <- c(5, 10, 13, 16, 19, 22, 25)
  
  lda_primary <- function(fit) apply(posterior(fit)$topics, 1, which.max)
  
  lda_stability <- function(k, dtm = lda_dtm, seeds = lda_seeds, iter = 500) {
    fits <- lapply(seeds, \(s) LDA(dtm, k = k, method = "Gibbs",
                                   control = list(seed = s, iter = iter)))
    asg  <- lapply(fits, lda_primary)
    pr   <- combn(length(seeds), 2)
    data.frame(k = k,
               pair = combn(seeds, 2, paste, collapse = " vs "),
               ari = apply(pr, 2, \(i) mclust::adjustedRandIndex(asg[[i[1]]], asg[[i[2]]])))
  }
  
  lda_stab <- map_dfr(lda_K_check, lda_stability)
  lda_stab_sum <- lda_stab |>
    group_by(k) |>
    summarise(mean_ari = mean(ari), min_ari = min(ari), max_ari = max(ari))
  print(lda_stab_sum)
  
  ggplot(lda_stab_sum, aes(k, mean_ari)) +
    geom_ribbon(aes(ymin = min_ari, ymax = max_ari), alpha = .2) +
    geom_hline(yintercept = 0.3, linetype = "dashed", colour = "grey40") +
    geom_line() + geom_point() +
    labs(x = "Number of topics (k)",
         y = "Pairwise ARI across seeds (higher = more stable)",
         title = "LDA seed stability by K",
         subtitle = "Dashed line = conventional 'fair agreement' threshold (~0.3)") +
    theme_minimal(base_size = 15) + theme_fig()
  ggsave_fig("img/appendix_lda_stability_ari.png", width = 8, height = 6)
  
  ## ---- exploratory (not reported): K = 5 model, top terms, word cloud ----
  if (RUN_EXTRA) {
    library(wordcloud); library(RColorBrewer)
    lda_k5 <- LDA(lda_dtm, 5, method = "Gibbs", control = list(seed = 1989, iter = 500))
    print(terms(lda_k5, 10))
    top40 <- sort(posterior(lda_k5)$terms[5, ], decreasing = TRUE)[1:40]
    wordcloud(names(top40), top40, random.order = FALSE, colors = brewer.pal(8, "Dark2"))
  }
}


###############################################################################
## Appendix G  PCA ROBUSTNESS (rotated solution, as in the paper)
###############################################################################
if (RUN_PCA) {
  
  rb_vars <- setdiff(keep_vars, "id")                         # the 16 PCA inputs
  rb_ref  <- rot_load[rb_vars, ]                              # main rotated loadings
  
  ## ---- helpers ----
  # Tucker congruence between loading columns (common variables only)
  rb_congr <- function(A, B) {
    v <- intersect(rownames(A), rownames(B)); A <- A[v, , drop = FALSE]; B <- B[v, , drop = FALSE]
    crossprod(A, B) / sqrt(outer(colSums(A^2), colSums(B^2)))
  }
  
  # Fit PCA + rotation on a subset; align columns and signs to the main solution
  # by congruence (works even when affiliation or joy is the variable dropped)
  rb_fit <- function(songs, vars, rotation = c("varimax", "none", "promax")) {
    rotation <- match.arg(rotation)
    d  <- songs |> select(id, all_of(vars)) |> filter(if_all(all_of(vars), is.finite))
    stopifnot(nrow(d) == nrow(songs))                         # same songs in every refit
    p  <- prcomp(d |> select(-id), center = TRUE, scale. = TRUE)
    L  <- p$rotation[, 1:2] %*% diag(p$sdev[1:2])
    S  <- p$x[, 1:2] %*% diag(1 / p$sdev[1:2])
    phi <- NA_real_
    if (rotation == "varimax") { R <- varimax(L)$rotmat; L <- L %*% R; S <- S %*% R }
    if (rotation == "promax")  { pm <- psych::Promax(L); U <- pm$rotmat
    L <- L %*% U; S <- S %*% t(solve(U)); phi <- abs(pm$Phi[1, 2]) }
    C <- rb_congr(L, rb_ref)
    if (abs(C[1, 2]) + abs(C[2, 1]) > abs(C[1, 1]) + abs(C[2, 2])) {
      L <- L[, 2:1]; S <- S[, 2:1]; C <- C[2:1, ]
    }
    for (j in 1:2) if (C[j, j] < 0) { L[, j] <- -L[, j]; S[, j] <- -S[, j] }
    colnames(L) <- colnames(S) <- c("PC1", "PC2")
    list(L = L, S = S, id = d$id, phi = phi,
         ve_tot = sum(p$sdev[1:2]^2) / sum(p$sdev^2),
         congr = diag(rb_congr(L, rb_ref)))
  }
  
  # Party tests on campaign means of the scores
  rb_party <- function(f, songs) {
    cm <- data.frame(id = f$id, d1 = f$S[, 1], d2 = f$S[, 2]) |>
      inner_join(songs |> select(id, cand_full, campaign, party), by = "id") |>
      mutate(cy = paste(word(as.character(cand_full), -1), str_extract(campaign, "\\d{4}"))) |>
      group_by(cy, party) |> summarise(c1 = mean(d1), c2 = mean(d2), .groups = "drop")
    t1 <- t.test(c1 ~ party, data = cm); t2 <- t.test(c2 ~ party, data = cm)
    tibble(n_R = sum(cm$party == "R"), n_D = sum(cm$party == "D"),
           dim1_gap = unname(diff(t1$estimate)), dim1_p = t1$p.value,   # R minus D
           dim2_gap = unname(diff(t2$estimate)), dim2_p = t2$p.value)
  }
  
  rb_main <- rb_fit(cand_songs, rb_vars)
  stopifnot(all(rb_main$congr > 0.999))                       # reproduces script 1
  rb_main_party <- rb_party(rb_main, cand_songs)
  print(rb_main_party)
  
  ## ---- G.1 Leave one input variable out ----
  rb_lovo <- map_dfr(rb_vars, function(v) {
    f <- rb_fit(cand_songs, setdiff(rb_vars, v))
    bind_cols(tibble(dropped = v, congr_dim1 = f$congr[1], congr_dim2 = f$congr[2]),
              rb_party(f, cand_songs))
  })
  print(rb_lovo, n = Inf)
  
  rb_lovo |>
    mutate(dropped = coalesce(var_labels[dropped], dropped)) |>
    pivot_longer(c(congr_dim1, congr_dim2), names_to = "dim", values_to = "congr") |>
    mutate(dim = recode(dim, congr_dim1 = "Dimension 1", congr_dim2 = "Dimension 2"),
           dropped = fct_reorder(dropped, congr, .fun = min)) |>
    ggplot(aes(congr, dropped, colour = dim, shape = dim)) +
    geom_vline(xintercept = c(0.85, 0.95), linetype = "dashed", colour = "grey60") +
    geom_point(size = 3) +
    scale_shape_manual(values = c("Dimension 1" = 16, "Dimension 2" = 17), name = NULL) +
    scale_colour_manual(values = c("Dimension 1" = "grey15", "Dimension 2" = "grey60"),
                        name = NULL) +
    labs(x = "Tucker congruence with main solution", y = "Variable omitted",
         caption = "Dashed lines: .85 (fair similarity) and .95 (equal) (Lorenzo-Seva & ten Berge, 2006).") +
    theme_minimal(base_size = 15) + theme(legend.position = "bottom") + theme_fig()
  ggsave_fig("img/pca_loo_similarity.png", width = 8, height = 7)
  
  rb_lovo |>
    transmute(`Variable omitted` = coalesce(var_labels[dropped], dropped),
              `Congr. Dim 1` = congr_dim1, `Congr. Dim 2` = congr_dim2,
              `Dim 1 gap (R$-$D)` = dim1_gap, `Dim 1 $p$` = fmt_p(dim1_p),
              `Dim 2 gap (R$-$D)` = dim2_gap, `Dim 2 $p$` = fmt_p(dim2_p)) |>
    xtable(digits = c(0, 0, 2, 2, 2, 0, 2, 0)) |>
    print(file = "img/pca_var_partygap.tex", type = "latex", include.rownames = FALSE,
          floating = FALSE, sanitize.colnames.function = identity)
  
  ## ---- G.2 Leave one campaign / one candidate out (PCA re-estimated each time) ----
  rb_songs <- cand_songs |>
    mutate(cy   = paste(word(as.character(cand_full), -1), str_extract(campaign, "\\d{4}")),
           cand = word(as.character(cand_full), -1))
  # campaign ids match script 1, and every campaign meets the six-song rule
  stopifnot(all(rb_songs$cy == as.character(rb_songs$cand_year)),
            min(dplyr::count(rb_songs, cy)$n) >= 6)

  rb_multi <- rb_songs |> distinct(cand, cy) |> dplyr::count(cand) |> filter(n > 1) |> pull(cand)
  rb_drop  <- bind_rows(tibble(type = "Campaign",  unit = sort(unique(rb_songs$cy))),
                        tibble(type = "Candidate", unit = rb_multi))   # e.g. all three Trump campaigns
  
  rb_loco <- pmap_dfr(rb_drop, function(type, unit) {
    s <- if (type == "Campaign") filter(rb_songs, cy != unit) else filter(rb_songs, cand != unit)
    f <- rb_fit(s, rb_vars)
    bind_cols(tibble(type = type, dropped = unit,
                     congr_dim1 = f$congr[1], congr_dim2 = f$congr[2]),
              rb_party(f, s))
  })
  print(rb_loco, n = Inf)
  
  rb_loco |>
    transmute(Omitted = paste0(dropped, if_else(type == "Candidate", " (all campaigns)", "")),
              `Congr. Dim 1` = congr_dim1, `Congr. Dim 2` = congr_dim2,
              `Dim 1 $p$` = fmt_p(dim1_p),
              `Dim 2 gap (R$-$D)` = dim2_gap, `Dim 2 $p$` = fmt_p(dim2_p)) |>
    xtable(digits = c(0, 0, 2, 2, 0, 2, 0)) |>
    print(file = "img/pca_loco.tex", type = "latex", include.rownames = FALSE,
          floating = FALSE, sanitize.colnames.function = identity)
  
  ## Summary lines for the Section 4.5 / Appendix G text
  cat("\nLeave-one-variable-out: min congruence Dim 1 =", round(min(rb_lovo$congr_dim1), 2),
      "| Dim 2 =", round(min(rb_lovo$congr_dim2), 2),
      "| Dim 2 p range", fmt_p(min(rb_lovo$dim2_p)), "-", fmt_p(max(rb_lovo$dim2_p)),
      "| Dim 1 significant in", sum(rb_lovo$dim1_p < .05), "of", nrow(rb_lovo), "\n")
  cat("Leave-one-out (campaign/candidate): min congruence Dim 1 =", round(min(rb_loco$congr_dim1), 2),
      "| Dim 2 =", round(min(rb_loco$congr_dim2), 2),
      "| Dim 2 p max", fmt_p(max(rb_loco$dim2_p)),
      "| Dim 1 significant in", sum(rb_loco$dim1_p < .05), "of", nrow(rb_loco), "\n")
  
  ## ---- G.3 Specification comparison: I - we x rotation (Table A18) ----
  rb_spec_grid <- expand_grid(iwe = c("with", "without"),
                              rotation = c("none", "varimax", "promax"))
  rb_spec <- pmap_dfr(rb_spec_grid, function(iwe, rotation) {
    v <- if (iwe == "with") rb_vars else setdiff(rb_vars, "i_minus_we")
    f <- rb_fit(cand_songs, v, rotation)
    bind_cols(tibble(iwe = iwe, rotation = rotation, var_expl = f$ve_tot, promax_phi = f$phi,
                     congr_dim1 = f$congr[1], congr_dim2 = f$congr[2]),
              rb_party(f, cand_songs))
  })
  print(rb_spec, width = Inf)
  
  rb_spec |>
    transmute(`I $-$ we` = iwe, Rotation = rotation,
              `Var. (\\%)` = 100 * var_expl, `Promax $r$` = promax_phi,
              `Congr. Dim 1` = congr_dim1, `Congr. Dim 2` = congr_dim2,
              `Dim 1 $p$` = fmt_p(dim1_p), `Dim 2 $p$` = fmt_p(dim2_p)) |>
    xtable(digits = c(0, 0, 0, 1, 2, 2, 2, 0, 0)) |>
    print(file = "img/pca_spec_compare.tex", type = "latex", include.rownames = FALSE,
          floating = FALSE, sanitize.colnames.function = identity, NA.string = "")
  
  ## Correlations among pronoun-based inputs (Section 4.5 text)
  print(round(cor(num_data[c("Clout", "Authentic", "i_minus_we", "affiliation")]), 2))
  
  ## How well the 2-D map represents each input (supplement sentence)
  rb_R   <- cor(num_data[rb_vars])
  rb_off <- (rb_R - rb_ref %*% t(rb_ref))[upper.tri(rb_R)]
  cat("\nCommunalities (rotated 2-D solution):\n"); print(round(sort(rowSums(rb_ref^2)), 2))
  for (k in 2:4) {
    Lk <- pca_res$rotation[, 1:k] %*% diag(pca_res$sdev[1:k])
    o  <- (rb_R - Lk %*% t(Lk))[upper.tri(rb_R)]
    cat(k, "components: residuals >= |.10| =", round(mean(abs(o) >= .10), 2),
        "  RMSR =", round(sqrt(mean(o^2)), 3), "\n")
  }
  
  ## ---- G.4 Figure A2: Dimension 3 (x) by rotated Dimension 2 (y) ----
  ## Dim 3 is the third unrotated component, standardized like the rotated scores.
  ## It is orthogonal to the rotated Dim 1-2 plane.
  a2_load3 <- pca_res$rotation[, 3] * pca_res$sdev[3]
  print(round(sort(a2_load3), 2))      # use these to label Dim 3 (set a2_xlab below)
  
  a2_xlab <- "Dim 3"                    # e.g. "Dim 3: <low pole> to <high pole>"
  a2_ylab <- "Dim 2: Forward-looking to Golden Era"
  
  a2_songs <- data.frame(id = num_data$id,
                         d3 = pca_res$x[, 3] / pca_res$sdev[3],
                         d2 = rot_scores[, 2]) |>
    inner_join(cand_songs |> select(id, cand_full, campaign, party), by = "id") |>
    mutate(party = as.character(party),
           cand_year = paste(word(as.character(cand_full), -1), str_extract(campaign, "\\d{4}")))
  # caption says six or more songs per campaign: hold the code to it
  stopifnot(min(dplyr::count(a2_songs, cand_year)$n) >= 6)

  a2_B <- 2000
  a2_cent <- a2_songs |>
    group_by(cand_year, party) |>
    summarise(n_songs = n(), c1 = mean(d3), c2 = mean(d2), .groups = "drop")
  
  set.seed(1389)
  a2_boot <- map_dfr(seq_len(a2_B), function(b) {
    a2_songs |> group_by(cand_year) |> slice_sample(prop = 1, replace = TRUE) |>
      summarise(m1 = mean(d3), m2 = mean(d2), .groups = "drop")
  })
  a2_stab <- a2_boot |>
    group_by(cand_year) |>
    summarise(lo1 = quantile(m1, .025), hi1 = quantile(m1, .975),
              lo2 = quantile(m2, .025), hi2 = quantile(m2, .975), .groups = "drop") |>
    left_join(a2_cent, by = "cand_year")
  
  a2_by_camp <- split(a2_songs, a2_songs$cand_year)
  a2_party_cent <- a2_cent |> group_by(party) |>
    summarise(n_camp = n(), c1 = mean(c1), c2 = mean(c2), .groups = "drop")
  
  set.seed(8987)
  a2_party_boot <- map_dfr(seq_len(a2_B), function(b) {
    map_dfr(unique(a2_cent$party), function(pp) {
      draw <- sample(a2_cent$cand_year[a2_cent$party == pp], replace = TRUE)
      mm <- vapply(draw, function(cc) {
        s <- a2_by_camp[[cc]]; i <- sample.int(nrow(s), nrow(s), replace = TRUE)
        c(mean(s$d3[i]), mean(s$d2[i]))
      }, numeric(2))
      tibble(party = pp, p1 = mean(mm[1, ]), p2 = mean(mm[2, ]))
    })
  })
  a2_party_stab <- a2_party_boot |>
    group_by(party) |>
    summarise(lo1 = quantile(p1, .025), hi1 = quantile(p1, .975),
              lo2 = quantile(p2, .025), hi2 = quantile(p2, .975), .groups = "drop") |>
    left_join(a2_party_cent, by = "party") |>
    mutate(lab = paste0(party, " (", n_camp, " campaigns)"))
  
  # Descriptive only: songs are not independent draws, so no ANOVA p-values.
  # omega2 corrects eta2 for chance (eta2 averages ~df_b / df_tot, ~7%, under no effect)
  a2_ve <- map_dfr(c("d3", "d2"), function(v) {
    a <- anova(lm(reformulate("cand_year", v), data = a2_songs))
    ss <- a$`Sum Sq`; ms_e <- a$`Mean Sq`[2]
    tibble(eta2   = ss[1] / sum(ss),
           omega2 = max(0, (ss[1] - a$Df[1] * ms_e) / (sum(ss) + ms_e)))
  })
  print(a2_ve)
  print(t.test(c1 ~ party, data = a2_cent))    # party difference on Dim 3
  
  a2_cap <- paste0(
    "Campaigns with six or more songs. Small points are campaign means ",
    "(crossbars: 95% bootstrap CI over songs).\n",
    "Large circles (D) and triangles (R) are party means of campaign means (crossbars: nested bootstrap over ",
    "campaigns and songs).\n",
    "Beyond chance, campaign membership accounts for ", pct(a2_ve$omega2[1]), " and ",
    pct(a2_ve$omega2[2]), " of song-level variation in Dimension 3 and Dimension 2 scores ",
    "(omega squared, descriptive).")
  
  # party labels sit just off their markers, joined by a leader line
  a2_party_lab <- a2_party_stab |>
    mutate(lx = if_else(party == "R", c1 - 0.25, c1 - 0.30),
           ly = if_else(party == "R", c2 + 0.70, c2 - 0.65))

  ggplot(a2_stab, aes(c1, c2, colour = party)) +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_hline(yintercept = 0, colour = "grey88") +
    # campaigns: faded crossbars
    geom_segment(aes(x = lo1, xend = hi1, y = c2, yend = c2, colour = I(party_outline[party])),
                 alpha = 0.25, linewidth = 0.5) +
    geom_segment(aes(x = c1, xend = c1, y = lo2, yend = hi2, colour = I(party_outline[party])),
                 alpha = 0.25, linewidth = 0.5) +
    # parties: heavy crossbars, large circle (D) / triangle (R)
    geom_segment(data = a2_party_stab, linewidth = 1.1, alpha = 0.8, aes(x = lo1, xend = hi1, y = c2, yend = c2)) +
    geom_segment(data = a2_party_stab, linewidth = 1.1, alpha = 0.8, aes(x = c1, xend = c1, y = lo2, yend = hi2)) +
    geom_point(data = a2_party_stab, aes(c1, c2, shape = party), size = 7) +
    # leader lines drawn under the campaign names (halos hide crossings)
    geom_segment(data = a2_party_lab, aes(x = c1, y = c2, xend = lx, yend = ly), linewidth = 0.4) +
    # campaigns: clear points; labels with a white halo, no box
    geom_point(aes(shape = party), size = 3.2) +
    ggrepel::geom_text_repel(aes(label = cand_year), size = fig_text, seed = 13,
                             bg.color = "white", bg.r = 0.15,
                             min.segment.length = 0.2, box.padding = 0.3, show.legend = FALSE) +
    geom_label(data = a2_party_lab, aes(lx, ly, label = lab), fontface = "bold", size = fig_text,
               fill = "white", linewidth = 0, show.legend = FALSE) +
    scale_colour_manual(values = party_outline, name = "Party") +  # darker shades: legible in B&W
    scale_shape_manual(values = c(D = 16, R = 17), name = "Party") +
    labs(x = a2_xlab, y = a2_ylab,
         title = "Campaign and party centroids with bootstrap uncertainty (Dimensions 3 and 2)",
         caption = a2_cap) +
    coord_cartesian(xlim = c(-0.6, 0.85), ylim = c(-1.05, 1.25)) +  # zoomed to campaign means
    theme_minimal(base_size = 15) +
    theme(plot.caption = element_text(hjust = 0, colour = "grey30"),
          legend.position = "bottom") +
    theme_fig()
  ggsave_fig("img/pca_cand_boot_dim2_dim3.png", width = 8, height = 9.5)
}

