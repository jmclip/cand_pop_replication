###### TECHNICAL APPENDIX
# Showcase 3 approaches: 
# 1) BOW
# 2) LDA 
# 3) Embeddings (separate file: .py scripts: embeddings and bertopic)

# set options
options(stringsAsFactors = F)         # no automatic data transformation
options("scipen" = 100, "digits" = 4) # suppress math annotation
# load packages
library(knitr) 
library(kableExtra) 
library(DT)
library(tm)
library(topicmodels)
library(reshape2)
library(ggplot2)
library(wordcloud)
library(pals)
library(SnowballC)
library(stm)
# library(lda)
# library(ldatuning)
#library(flextable)
library(tidyverse)
library(text2vec)
library(glmnet)
library(rsample)
library(mclust)

## read data
cand_songs_appendix <- read_csv("data/en_songs_liwc.csv") |> 
  mutate(cand_year = fct_reorder(cand_year, year)) |>
  mutate(across(starts_with("NRC_"), ~ 100 * .x / WC), 
         across(starts_with("avg_"), ~ 100 * .x / WC),
         i_less_we = i-we, fem_male = female-male) |>  #deal with NRC wordcount issues for scaling
  # Order by 'num_plays' (or swap with 'WC' or another metric)
  group_by(cand_year) |> filter(song_ct>5) |> ungroup() 

# ####### BOW and EMBEDDDINGS #########
## BOW approach: https://cran.r-project.org/web/packages/text2vec/vignettes/text-vectorization.html

## PREP
# Convert cand_songs_appendix to a data.table and set 'id' as the key
# 1. Clean and assign train_id to the FULL dataset ONCE
set.seed(1389)
NFOLDS = 5

cand_songs_embed <- cand_songs_appendix |>
  mutate(party01 = if_else(party == "D", 1, 0)) |>
  filter(!is.na(genre) & nchar(lyrics_liwc) > 50) |>
  distinct(title, performer, party, .keep_all = TRUE) |>
  mutate(train_id = as.character(row_number())) # Explicit character ID

# 2. Split into train/test objects -- KEEP SONGS TOGETHER!
split <- group_initial_split(cand_songs_embed, group = title, prop = 0.80)

# Use these directly
train <- training(split)
test  <- testing(split)

# 3. Create iterator directly from 'train'
prep_fun = tolower
tok_fun  = word_tokenizer

it_train <- itoken(
  train$lyrics_liwc,
  preprocessor = tolower,
  tokenizer = word_tokenizer,
  ids = train$train_id,
  progressbar = FALSE
)

it_test <- itoken(
  test$lyrics_liwc,
  preprocessor = tolower,
  tokenizer = word_tokenizer,
  ids = test$train_id,
  progressbar = FALSE
)


####### fitting the first model #####
vocab = create_vocabulary(it_train)
vocab

vectorizer = vocab_vectorizer(vocab)
t1 = Sys.time()
dtm_train = create_dtm(it_train, vectorizer)
print(difftime(Sys.time(), t1, units = 'sec'))

# check:
dim(dtm_train)
train$train_id <- as.character(train$train_id)

identical(rownames(dtm_train), as.character(train$train_id))

# maybe a problem re: as character
identical(rownames(dtm_train), train$train_id)

glmnet_classifier = cv.glmnet(x = dtm_train, y = train[['party01']],
                              family = 'binomial',
                              # L1 penalty
                              alpha = 1,
                              # interested in the area under ROC curve
                              type.measure = "auc",
                              # 5-fold cross-validation
                              nfolds = NFOLDS,
                              # high value is less accurate, but has faster training
                              thresh = 1e-3,
                              # again lower number of iterations for faster training
                              maxit = 1e3)

plot(glmnet_classifier)

print(paste("max AUC =", round(max(glmnet_classifier$cvm), 4)))

## how we are doing overall
dtm_test = create_dtm(it_test, vectorizer)
preds = predict(glmnet_classifier, dtm_test, type = 'response')[,1]
glmnet:::auc(test$party01, preds)
pROC::ci.auc(pROC::roc(test$party01, preds))


## stability of model:
seed_runs <- c(1389, 13, 87, 89, 2026, 1015,216,410,2012,2016,1983,1130,12,16,83,84)
auc_split_stability <- function(seed, data = cand_songs_embed, vectorizer, prop = 0.80) {
  set.seed(seed)
  split <- group_initial_split(data, group = title, prop = prop)
  train <- training(split); test <- testing(split)
  
  it_train <- itoken(train$lyrics_liwc, preprocessor = tolower,
                     tokenizer = word_tokenizer, ids = train$train_id, progressbar = FALSE)
  it_test  <- itoken(test$lyrics_liwc, preprocessor = tolower,
                     tokenizer = word_tokenizer, ids = test$train_id, progressbar = FALSE)
  
  dtm_train <- create_dtm(it_train, vectorizer)
  dtm_test  <- create_dtm(it_test, vectorizer)
  
  fit   <- cv.glmnet(x = dtm_train, y = train[['party01']], family = 'binomial',
                     alpha = 1, type.measure = "auc", nfolds = NFOLDS,
                     thresh = 1e-3, maxit = 1e3)
  preds <- predict(fit, dtm_test, type = 'response')[, 1]
  
  data.frame(seed = seed, cv_auc = max(fit$cvm),
             test_auc = as.numeric(glmnet:::auc(test$party01, preds)))
}

split_results <- do.call(rbind, lapply(seed_runs, auc_split_stability,
                                       vectorizer = vectorizer))
split_results

split_summary_row <- split_results |>
  summarise(
    seed     = "Mean (SD)",
    cv_auc   = sprintf("%.3f (%.3f)", mean(cv_auc), sd(cv_auc)),
    test_auc = sprintf("%.3f (%.3f)", mean(test_auc), sd(test_auc))
  )

split_range_row <- split_results |>
  summarise(
    seed     = "Range",
    cv_auc   = sprintf("%.3f – %.3f", min(cv_auc), max(cv_auc)),
    test_auc = sprintf("%.3f – %.3f", min(test_auc), max(test_auc))
  )

split_results_fmt <- split_results |>
  mutate(across(c(cv_auc, test_auc), ~ sprintf("%.3f", .x)),
         seed = as.character(seed))

table_data <- bind_rows(split_results_fmt, split_summary_row, split_range_row)

table_data |>
  kbl(col.names = c("Seed", "CV AUC (train)", "Test AUC"),
      align = "c",
      caption = "BOW classifier AUC across random train/test splits") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) |>
  row_spec(nrow(split_results_fmt), extra_css = "border-top: 2px solid black;")

split_long <- split_results |>
  pivot_longer(c(cv_auc, test_auc), names_to = "metric", values_to = "auc") |>
  mutate(metric = recode(metric, cv_auc = "CV AUC (train)", test_auc = "Test AUC"))

ggplot(split_long, aes(x = metric, y = auc)) +
  geom_jitter(width = 0.08, height = 0, alpha = .6, size = 2) +
  geom_boxplot(fill="transparent")+
  stat_summary(fun = mean,  geom = "point", size = 3,color = "firebrick") +
  stat_summary(fun = mean, geom = "text", label = "mean", hjust = -0.5,vjust = 0.5,   
    color = "firebrick") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(x = NULL, y = "AUC",
       title = "BOW classifier performance across 15 random splits",
       subtitle = "Dashed line = chance (AUC = 0.5)") +
  theme_minimal()
ggsave(filename = "img/appendix_bow_sensitivity.png", width = 8)

#### VERSION TWO refining the BOW model ####
## refined model: stop words and ngrams
filter_stop <- c("good", "kind", "problem", "wanting", "important", "parting")
stop_add <- c("ass", "uh", "ooh", "hoo", "ah")
stop_words_filtered <- tidytext::stop_words |> filter(! word %in% filter_stop) |>
  pull(word) |>  c(stop_add)

vocab_nostop <- create_vocabulary(
  it_train,
  ngram = c(ngram_min = 1L, ngram_max = 3L),
  stopwords = stop_words_filtered)

pruned_vocab = prune_vocabulary(
  vocab_nostop,
  doc_count_min = 4,        # Drop rare phrases
  doc_proportion_max = 0.5   # Drop overly common phrases (phrases in over 50% of docs)
)

# Build DTM and rerun
vectorizer_nostop = vocab_vectorizer(pruned_vocab)
dtm_train_nostop = create_dtm(it_train, vectorizer_nostop)
glmnet_classifier_nostop = cv.glmnet(x = dtm_train_nostop, y = train[['party01']],
                              family = 'binomial',
                              # L1 penalty
                              alpha = 1,
                              # interested in the area under ROC curve
                              type.measure = "auc",
                              # 5-fold cross-validation
                              nfolds = NFOLDS)

plot(glmnet_classifier_nostop)

print(paste("max AUC =", round(max(glmnet_classifier_nostop$cvm), 4)))

## apply to test dataset:
# apply vectorizer
dtm_test_nostop = create_dtm(it_test, vectorizer_nostop)
preds = predict(glmnet_classifier_nostop, dtm_test_nostop, type = 'response')[,1]
glmnet:::auc(test$party01, preds)
pROC::ci.auc(pROC::roc(test$party01, preds))


######## plot
## plot of relevant words

# 1. Extract non-zero coefficients at lambda.1se
coef_matrix <- coef(glmnet_classifier_nostop, s = "lambda.1se")

# 2. Convert sparse matrix to a data frame with exact party labels
coef_df <- data.frame(
  term = rownames(coef_matrix),
  estimate = as.vector(coef_matrix)
) |>
  filter(term != "(Intercept)", estimate != 0) |>
  mutate(party_label = if_else(estimate > 0, "Democrat (1)", "Republican (0)"))

# 3. Get top 10 terms for each party
top_terms <- coef_df |>
  group_by(party_label) |>
  slice_max(order_by = abs(estimate), n = 10) |>
  ungroup()

# 4. Plot top predictive features
ggplot(top_terms, aes(x = reorder(term, estimate), y = estimate, fill = party_label)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  coord_flip() +
  facet_wrap(~ party_label, scales = "free") +
  scale_fill_manual(values = c(
    "Democrat (1)"           = "#2b5c8f", # Blue
    "Republican (0)" = "#a83232"  # Red
  )) +
  labs(
    x = NULL,
    y = "Coefficient Estimate (Log-Odds Ratio at lambda.1se)",
    title = "Top Predictive Language Features by Party",
    subtitle = "Positive = Democrat | Negative = Republican / Other",
    caption = "Extracted from Lasso Logistic Regression on Campaign Songs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 11))


######## LDA ############
# tutorial: https://slcladal.github.io/topicmodels.html
# https://tm4ss.github.io/docs/Tutorial_6_Topic_Models.html


# load data
textdata <- cand_songs_appendix %>% select(performer, title, cand_year, party, lyrics_liwc) %>%
  mutate(doc_id = as.character(row_number()), text = as.character(lyrics_liwc)) 


# load stopwords
english_stopwords <- readLines("https://slcladal.github.io/resources/stopwords_en.txt", encoding = "UTF-8")

# corpus object
corpus <- VCorpus(DataframeSource(textdata))

processedCorpus <- tm_map(corpus, content_transformer(tolower))
processedCorpus <- tm_map(processedCorpus, removePunctuation, preserve_intra_word_dashes = TRUE)
processedCorpus <- tm_map(processedCorpus, removeNumbers)
processedCorpus <- tm_map(processedCorpus, removeWords, english_stopwords)
processedCorpus <- tm_map(processedCorpus, stemDocument, language = "en")
processedCorpus <- tm_map(processedCorpus, stripWhitespace)

# compute document term matrix with terms >= minimumFrequency
minimumFrequency <- 5
DTM <- DocumentTermMatrix(processedCorpus, control = list(
  bounds = list(global = c(minimumFrequency, Inf)),
  stopwords = stemDocument(english_stopwords)  # filter stemmed forms at DTM stage
))
# have a look at the number of documents and terms in the matrix


remove_terms <- c(
  # original
  "ooh", "yeah", "nigga", "whoa", "uh-huh", "hey", "woah", "ooh-ooh", "ooh-ooh-ooh",
  # common words slipping through
  "the", "and", "you", "your", "for", "that", "but", "can", "just",
  "what", "all", "this", "from", "her", "when", "there", "one", "was",
  "like", "get", "got", "let", "feel", "know", "time", "man", "will", "with",
  # informal/lyric contractions standard stopword lists miss
  "dont", "wanna", "aint", "gonna", "gotta", "cant", "youre", "youll",
  "ive", "ill", "thats", "its", "hes", "shes", "theyre", "were",
  "cause", "cuz", "bout", "til", "em"
)
remove_terms_stemmed <- unique(c(remove_terms, stemDocument(remove_terms)))

DTM <- DTM[, !(colnames(DTM) %in% remove_terms_stemmed)]

dim(DTM)

# due to vocabulary pruning, we have empty rows in our DTM
# LDA does not like this. So we remove those docs from the
# DTM and the metadata

sel_idx <- slam::row_sums(DTM) > 0
DTM <- DTM[sel_idx, ]
textdata <- textdata[sel_idx, ]

## used claude to convert bc ldatuning no longer available
cao_juan <- function(model) {
  phi <- posterior(model)$terms
  pairs <- combn(nrow(phi), 2)
  mean(apply(pairs, 2, function(i) sum(phi[i[1],] * phi[i[2],]) /
               (sqrt(sum(phi[i[1],]^2)) * sqrt(sum(phi[i[2],]^2)))))
}

deveaud <- function(model) {
  phi <- posterior(model)$terms
  pairs <- combn(nrow(phi), 2)
  js <- function(p, q) { m <- 0.5*(p+q); 0.5*sum(p*log(p/m), na.rm=TRUE) + 0.5*sum(q*log(q/m), na.rm=TRUE) }
  mean(apply(pairs, 2, function(i) js(phi[i[1],], phi[i[2],])))
}

fit_lda <- function(k) LDA(DTM, k = k, method = "Gibbs", control = list(seed = 13))

models <- lapply(2:25, fit_lda)

result <- data.frame(
  topics      = 2:25,
  CaoJuan2009 = sapply(models, cao_juan),
  Deveaud2014 = sapply(models, deveaud)
)

# Look for max of D and min of C
result |>
  mutate(across(-topics, scale)) |>
  pivot_longer(-topics) |>
  ggplot(aes(topics, value)) +
  geom_line() + geom_point() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  labs(x = "Number of Topics", y = "Value (scaled)")


#### select k: apply insight ####
# number of topics
K <- 5
set.seed(1989)
topicModel <- LDA(DTM, K, method = "Gibbs", control = list(iter = 500, verbose = 25))

tmResult <- posterior(topicModel)
attributes(tmResult)

ncol(DTM)        # lengthOfVocab  (replaces nTerms)
beta <- tmResult$terms
dim(beta)
rowSums(beta)

theta <- tmResult$topics
dim(theta)
nrow(DTM)        # size of collection  (replaces nDocs)

terms(topicModel, 15)
exampleTermData <- terms(topicModel, 15)

top5termsPerTopic <- terms(topicModel, 5)
topicNames <- apply(top5termsPerTopic, 2, paste, collapse = " ")

topicToViz <- 5 # 
top40terms <- sort(tmResult$terms[topicToViz, ], decreasing = TRUE)[1:40]
words <- names(top40terms)
probabilities <- top40terms  # already sorted, no need to sort twice

mycolors <- brewer.pal(8, "Dark2")
wordcloud(words, probabilities, random.order = FALSE, color = mycolors)


## ADJUST alpha for different distribution within topics:
# see alpha from previous model
attr(topicModel, "alpha") 

topicModel2 <- LDA(DTM, K, method="Gibbs", control=list(iter = 500, verbose = 25, alpha = 0.2))

# top N keywords per topic as a readable table
N <- 10
keywords <- as.data.frame(t(terms(topicModel2, N)))
keywords$topic <- 1:K
keywords <- keywords[, c("topic", paste0("V", 1:N))]

# ----------------------------------------------------------------------
# FREX TERM RE-RANKING (stm::calcfrex)
# ----------------------------------------------------------------------
# 1. Compute FREX terms directly using stm's vector-based calculation
# (Works seamlessly on topicmodels LDA objects using posterior beta)
log_beta <- log(posterior(topicModel2)$terms)
beta <- exp(log_beta)

# Calculate FREX matrix (Topics x Terms)
w <- 0.5 # Weight parameter (0.5 balances frequency & exclusivity)
freq_rank <- apply(beta, 1, function(x) ecdf(x)(x))
excl_ratio <- beta / colSums(beta)[col(beta)]
excl_rank <- apply(excl_ratio, 1, function(x) ecdf(x)(x))
frex_matrix <- 1 / ((w / freq_rank) + ((1 - w) / excl_rank))

# 1. Extract top N FREX terms per topic (Yields a matrix of 10 rows x K columns)
terms_all <- colnames(beta)
top_frex_terms <- apply(frex_matrix, 1, function(row) {
  terms_all[order(row, decreasing = TRUE)[1:10]]
})

# 2. Transpose (t) so rows = K topics and columns = top 10 terms
frex_keywords <- as.data.frame(t(top_frex_terms))

# 3. Rename columns and attach topic numbers
colnames(frex_keywords) <- paste0("V", 1:10)
frex_keywords$topic <- 1:ncol(top_frex_terms)

# Reorder columns so 'topic' is first
frex_keywords <- frex_keywords[, c("topic", paste0("V", 1:10))]

print("--- TOP FREX TERMS PER TOPIC ---")
print(frex_keywords)


# Update topicNames using top 4 FREX terms for plot legends
topicNames <- apply(top_frex_terms[1:4, ], 2, paste, collapse = " ")
# ----------------------------------------------------------------------

# or tidy version, easier to read
terms(topicModel2, N) |>
  as.data.frame() |>
  setNames(1:K) |>
  pivot_longer(everything(), names_to = "topic", values_to = "term") |>
  mutate(topic = as.integer(topic)) |>
  group_by(topic) |>
  summarise(keywords = paste(term, collapse = ", ")) |>
  print(n = K)

tmResult <- posterior(topicModel2)
theta <- tmResult$topics
beta <- tmResult$terms
topicNames <- apply(terms(topicModel2, 5), 2, paste, collapse = " ")  # reset topicnames


# re-rank top topic terms for topic names
topicNames <- apply(terms(topicModel, 4), 2, paste, collapse = " ")

# mean topic proportions over all documents
topicProportions <- colSums(theta) / nrow(DTM)
names(topicProportions) <- topicNames

soP <- sort(topicProportions, decreasing = TRUE)
paste(round(soP, 3), ":", names(soP))

## distinctive topics:
# primary by cand:
primary_topic <- apply(theta, 1, which.max)
textdata$primary_topic <- topicNames[primary_topic]

table(textdata$cand_year, textdata$primary_topic)

## more info:
corpus_avg <- colMeans(theta)

topic_proportion_per_cand <- aggregate(theta, by = list(cand_year = textdata$cand_year), mean)

lift <- sweep(
  as.matrix(topic_proportion_per_cand[, -1]),
  MARGIN = 2,
  STATS  = corpus_avg,
  FUN    = "/"
)
rownames(lift) <- topic_proportion_per_cand$cand_year

# most distinctive topic per candidate
apply(lift, 1, function(x) topicNames[which.max(x)])

##
lift_df <- as.data.frame(lift)
colnames(lift_df) <- topicNames

lift_df |>
  tibble::rownames_to_column("cand_year") |>
  pivot_longer(-cand_year, names_to = "topic", values_to = "lift") |>
  group_by(cand_year) |>
  slice_max(lift, n = 3) |>
  arrange(cand_year)


### assessment:
## ---- k-fold perplexity: how well does the model predict held-out documents? ----
set.seed(1989)
K_grid  <- seq(2, 20, by = 2)
n_folds <- 5
fold_id <- sample(rep(1:n_folds, length.out = nrow(DTM)))

perplexity_cv <- function(k) {
  sapply(1:n_folds, function(f) {
    train_dtm <- DTM[fold_id != f, ]
    test_dtm  <- DTM[fold_id == f, ]
    # perplexity() requires newdata's terms to be a subset of the fitted
    # model's vocabulary -- restrict, then drop any doc left empty by that
    test_dtm <- test_dtm[, colnames(test_dtm) %in% colnames(train_dtm)]
    test_dtm <- test_dtm[slam::row_sums(test_dtm) > 0, ]
    
    fit <- LDA(train_dtm, k = k, method = "Gibbs",
               control = list(seed = 13, iter = 500))
    perplexity(fit, newdata = test_dtm)
  })
}

perp_results <- data.frame(
  k    = rep(K_grid, each = n_folds),
  fold = rep(1:n_folds, times = length(K_grid)),
  perplexity = unlist(lapply(K_grid, perplexity_cv))
)

perp_summary <- perp_results |>
  group_by(k) |>
  summarise(mean_perplexity = mean(perplexity), se = sd(perplexity) / sqrt(n()))

ggplot(perp_summary, aes(k, mean_perplexity)) +
  geom_ribbon(aes(ymin = mean_perplexity - se, ymax = mean_perplexity + se), alpha = .2) +
  geom_line() + geom_point() +
  labs(x = "Number of topics (k)", y = "Held-out perplexity (lower = better fit)",
       title = "5-fold cross-validated perplexity")

stm_corpus <- readCorpus(DTM, type = "slam")   # converts your existing tm DTM directly

storage <- searchK(
  documents = stm_corpus$documents,
  vocab     = stm_corpus$vocab,
  K         = seq(2, 25, by = 2),
  seed      = 1989)

k_metrics <- storage$results %>%
  select(K, heldout, residual, semcoh, exclus) %>%
  # 1. Unnest the list-columns into standard numeric vectors
  unnest(cols = c(K, heldout, residual, semcoh, exclus)) %>%
  pivot_longer(cols = -K, names_to = "metric", values_to = "value") %>%
  mutate(metric = case_when(
    metric == "heldout"  ~ "Held-out Likelihood (Higher = Better)",
    metric == "residual" ~ "Residuals (Lower = Better)",
    metric == "semcoh"   ~ "Semantic Coherence (Higher = Better)",
    metric == "exclus"   ~ "Exclusivity (Higher = Better)"
  ))

# Verify that K and value are now numeric (<dbl>) rather than lists (<list>)
head(k_metrics)

p_k <- ggplot(k_metrics, aes(x = K, y = value)) +
  geom_line(color = "grey40", linewidth = 0.7) +
  geom_point(color = "#2b5c8f", size = 2.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs( x = "Number of Topics (K)", y = NULL, title = "Diagnostic Values by Number of Topics (stm::searchK)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),strip.text = element_text(face = "bold"))

ggsave("img/appendix_lda_topic_k.png", plot = p_k, width = 8, height = 6, dpi = 300)

## ---- seed stability: does the LDA solution replicate across random seeds? ----
seeds <- c(13, 87, 89, 8789, 2026)

primary_topic <- function(fit) apply(posterior(fit)$topics, 1, which.max)

lda_stability <- function(k, dtm = DTM, seeds, iter = 500) {
  fits <- lapply(seeds, function(s)
    LDA(dtm, k = k, method = "Gibbs", control = list(seed = s, iter = iter)))
  
  assignments <- lapply(fits, primary_topic)
  pairs <- combn(length(seeds), 2)
  
  data.frame(
    k    = k,
    pair = combn(seeds, 2, paste, collapse = " vs "),
    ari  = round(apply(pairs, 2, function(i)
      mclust::adjustedRandIndex(assignments[[i[1]]], assignments[[i[2]]])), 3)
  )
}

K_check <- c(5, 10, 13, 16, 19, 22, 25)   # coherence check k

stability_results <- do.call(rbind, lapply(K_check, lda_stability, dtm = DTM, seeds = seeds))

stability_results

stability_summary <- stability_results |>
  group_by(k) |>
  summarise(mean_ari = mean(ari), min_ari = min(ari), max_ari = max(ari))

ggplot(stability_summary, aes(k, mean_ari)) +
  geom_ribbon(aes(ymin = min_ari, ymax = max_ari), alpha = .2) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey40") +
  geom_line() + geom_point() +
  labs(x = "Number of topics (k)", y = "Pairwise ARI across seeds (higher = more stable)",
       title = "LDA seed stability by K",
       subtitle = "Dashed line = conventional 'fair agreement' threshold (~0.3)") +
  theme_minimal()
ggsave(filename = "img/appendix_lda_stability_ari.png", width = 8)
