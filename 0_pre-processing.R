###############################################################################
## Campaign rally music — preprocessing
###############################################################################

library(xtable)
library(tidytext)
library(syuzhet)
library(psych)
library(ggrepel)
library(factoextra)
library(tm)
library(cld3)
library(tidyverse)
library(ggpattern)
library(textdata)


############################################
############ pre-processing ################
############################################
set.seed(13)

## load data
cand_songs_raw1 <- read_csv("data/DISTINCT_FULL_COMPLETED_CLEAN_DATASET_2026_08.csv") |>
  filter(event_type == "campaign rally") |>
  mutate(
    year      = as.integer(str_extract(campaign, "\\d{4}")),
    cand_year = str_c(cand_l, " ", year),
    cand_year = fct_reorder(cand_year, year) ) |>
  select(-c(date, event_title)) |>
  mutate(performer = recode(performer, "Brooks and Dunn" = "Brooks & Dunn"),
         genre = case_when(
           title == "United States of America's national anthem" ~ "patriotic",
           title == "This Land Is Your Land" ~ "folk",
           title == "Higher Love" ~ "pop/rock",    
           TRUE ~ genre)) |>
  # one canonical lyric per song: prefer a Genius scrape, otherwise the first row
  group_by(title, performer) |>
  mutate(canon_row = { g <- which(str_detect(url_used, "genius\\.com")); if (length(g)) g[1] else 1L },
         lyrics    = lyrics[canon_row],
         url_used  = url_used[canon_row]) |>
  ungroup() |>
  select(-canon_row)

###### OTHER CLEANING
# strip structural tags once; both downstream texts start from this
strip_tags <- function(x) {
  x |>
    str_remove_all("\\d*Embed\\s*$") |>
    str_remove_all("\\[[^\\]]*\\]") |>
    str_remove_all("\\{[^\\}]*\\}") |>
    str_replace_all("[\u2019\u2018\u02bc]", "'") |>      
    str_replace_all("[\u201c\u201d]", "\"")  
}

# for LANGUAGE DETECTION: repeats kept, so cld3 has enough text
flatten_lyrics <- function(x) {
  str_squish(str_replace_all(strip_tags(x), "[\\r\\n]+", " "))
}


clean_lyrics <- function(x) strip_tags(x)


norm_txt <- function(x) {
  x |>
    str_to_lower() |>
    str_replace_all("[‘’“”]", "'") |>   # curly quotes
    str_replace_all("[–—]", "-") |>               # en/em dashes
    str_remove_all("[^a-z0-9 ]") |>
    str_squish()
}


########## TEXT NORMALIZATION ##########

## a) invisible + curly characters (keeps real accents)
tidy_text <- function(x) {
  x |>
    str_replace_all(" ", " ") |>                 # non-breaking space
    str_replace_all("[‘’ʼ]", "'") |>   # curly apostrophes
    str_replace_all("[“”]", '"') |>         # curly quotes
    str_replace_all("[–—]", "-") |>         # en / em dash
    str_squish()
}

cand_songs_norm <- cand_songs_raw1 |>
  mutate(across(c(title, performer), tidy_text))

## b) the one mojibake iconv couldn't repair (verify against source)
cand_songs_norm <- cand_songs_norm |>
  mutate(
    title     = if_else(str_detect(title, "^ç‹‚æ°"), "狂気沈殿", title),
    performer = if_else(str_detect(performer, "å¦–ç²¾"),
                        "YUI, Tatsuya Tachibana, 妖精帝國", performer)
  )

## c) canonical titles -- re-run AFTER (a), since tidying apostrophes
##    can create new variants to collapse
strip_version <- function(x) {
  str_remove(x, regex(
    "\\s*[\\(\\[-]\\s*(\\d{4}\\s*)?(remaster(ed)?|remix|radio edit|single version|album version|live|mono|stereo|deluxe)\\b.*$",
    ignore_case = TRUE)) |> str_squish()
}

#manual fix for duplicates
title_map <- c(
  "I Got A Feeling"     = "I Gotta Feeling",
  "Hold On, I'm Comin'" = "Hold On, I'm Coming"
)

# ok get right title with all cleaned up (remove all the extra stuff)
canon <- cand_songs_norm |>
  mutate(t = coalesce(title_map[title], title) |> strip_version()) |>
  count(t, name = "n") |>
  mutate(key = norm_txt(t)) |>
  group_by(key) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(key, title_canon = t)

# align titles to all match each other when they are the same
cand_songs_canon <- cand_songs_norm |>
  mutate(key = norm_txt(coalesce(title_map[title], title) |> strip_version())) |>
  left_join(canon, by = "key") |>
  mutate(title = coalesce(title_canon, title)) |>
  select(-key, -title_canon)



## Work on lyrics
cand_songs_cleaned <- cand_songs_canon |>
  mutate(
    is_instrumental = str_detect(str_squish(str_to_upper(lyrics)),
                                 "^\\[?INSTRUMENTAL\\]?$") %in% TRUE,
    
    lyrics_detect = flatten_lyrics(lyrics),        # must come before case_when
    n_words       = str_count(lyrics_detect, "\\S+"),
    lyrics_clean  = clean_lyrics(lyrics),
    
    lang_song = case_when(
      is_instrumental      ~ "instrumental",
      is.na(lyrics_detect) ~ NA_character_,
      n_words < 20         ~ "too_short",
      TRUE                 ~ detect_language(lyrics_detect)
    )
  )

count(cand_songs_cleaned, lang_song, sort = TRUE)

# any weird things going on here re missingness, etc.[all clear!]
bad_keys <- cand_songs_cleaned |>
  filter(!is.na(lyrics_clean), nchar(lyrics_clean) > 100) |>
  mutate(lyr_key = norm_txt(lyrics_clean)) |>
  group_by(lyr_key) |>
  filter(n_distinct(title) > 3) |>
  pull(lyr_key) |> unique()

cand_songs_cleanedf <- cand_songs_cleaned |>
  mutate(
    bad_scrape    = !is.na(lyrics_clean) & norm_txt(lyrics_clean) %in% bad_keys, #nothing here
    lyrics        = if_else(bad_scrape, NA_character_, lyrics),
    lyrics_clean  = if_else(bad_scrape, NA_character_, lyrics_clean),
    lyrics_detect = if_else(bad_scrape, NA_character_, lyrics_detect),
    lang_song     = if_else(bad_scrape, NA_character_, lang_song)
  )

# removed bad scrape (missing lyrics)
cat(sum(cand_songs_cleanedf$bad_scrape), "rows blanked as bad scrapes\n")


########## TYPO PROCESSING
## clean typos
cand_songs_cleanedf <- cand_songs_cleanedf |>
  mutate(genre = str_replace(genre,
                             regex("^hip[- ]?hop/rap$", ignore_case = TRUE),
                             "hip-hop/rap"))
#now combine duplicates
## FIXING


## CHECK ITS FIXED TEN TO FIX!!
## 1. NON-ASCII CHARACTERS in titles and performers -------------------
## Curly quotes, em dashes, non-breaking spaces, zero-width joiners.
cand_songs_cleanedf |>
  distinct(title, performer) |>
  filter(str_detect(title, "[^\\x20-\\x7E]") |
           str_detect(performer, "[^\\x20-\\x7E]")) |>
  mutate(
    odd_chars = map_chr(paste(title, performer), ~
                          paste(unique(unlist(str_extract_all(.x, "[^\\x20-\\x7E]"))), collapse = " "))
  ) |>
  print(n = Inf)

## 2. TITLES THAT NORMALIZE TO THE SAME STRING ------------------------[all fixed!]
## Catches "Don't Stop" vs "Don’t Stop", "Hey, Ya!" vs "Hey Ya".
cand_songs_cleanedf |>
  distinct(title, performer) |>
  mutate(key = norm_txt(title)) |>
  group_by(key) |>
  filter(n_distinct(title) > 1) |>
  summarise(variants = paste(unique(title), collapse = "  ||  "),
            performers = paste(unique(performer), collapse = "  ||  "),
            .groups = "drop") |>
  print(n = Inf)

## 3. SAME LYRICS, DIFFERENT TITLES -----------------------------------
## The strongest test: if the lyrics are identical it's the same song,
## whatever the title says.
cand_songs_cleanedf |>
  filter(!is.na(lyrics_clean), nchar(lyrics_clean) > 100) |>
  mutate(lyr_key = norm_txt(lyrics_clean)) |>
  group_by(lyr_key) |>
  filter(n_distinct(title) > 1) |>
  summarise(
    n_rows     = n(),
    titles     = paste(sort(unique(title)),     collapse = "  ||  "),
    performers = paste(sort(unique(performer)), collapse = "  ||  "),
    snippet    = str_trunc(first(lyrics_clean), 60),
    .groups    = "drop"
  ) |>
  select(-lyr_key) |>
  arrange(desc(n_rows)) |>
  print(n = Inf, width = Inf)


### LANGUAGE
## fix coding for songs not recognized (claude)
## NOTE: joined on title alone. Generic titles ("Free", "America", "99.9",
## "Do You Love Me") could match more than one song — check 7b below and add a
## performer column here if any turn out to be ambiguous in your data.
lang_fix <- tribble(
  ~title,                            ~lang_manual,
  # --- English ---------------------------------------------------
  "Whole Lotta Love",                "en",
  "Please, Please, Please",          "en",
  "Free",                            "en", 
  "RITMO",                           "en"
)

## ---- title key ------------------------------------------------------------
## (title_key isn't consumed anywhere yet — kept for downstream matching)
VW <- paste(c(
  "instrumental", "remix", "live", "remaster(?:ed)?", "radio edit",
  "single version", "album version", "extended", "edit", "mix", "acoustic",
  "demo", "mono", "stereo", "sessions?", "karaoke", "clean", "explicit",
  "deluxe", "bonus track", "re-?recorded", "feat\\.?", "featuring", "version"
), collapse = "|")

clean_title <- function(x) {
  x |>
    str_to_lower() |>
    str_replace_all("[‘’‛]", "'") |>
    str_remove_all(regex(paste0("\\s*[\\(\\[][^\\)\\]]*\\b(?:", VW, ")\\b[^\\)\\]]*[\\)\\]]"),
                         ignore_case = TRUE)) |>
    str_remove(regex(paste0("\\s*-\\s*(?:", VW, ").*$"), ignore_case = TRUE)) |>
    str_squish()
}

## ---- genre backfill -------------------------------------------------------

## fill with the modal non-NA genre for that song
fill_genre <- function(g) {
  if (all(is.na(g))) return(g)
  coalesce(g, names(sort(table(g), decreasing = TRUE))[1])
}

cand_songs_cleanedf <- cand_songs_cleanedf |>
  group_by(title, performer) |>
  mutate(genre = fill_genre(genre)) |>
  ungroup()

genre_fill <- tribble(
  ~title,                                                ~performer,                    ~genre_new,
  "Ain't Nothing Like The Real Thing", "Marvin Gaye & Tammi Terrell", "r&b",
  "America", "Los Tigres Del Norte", "latin",
  "Awake My Soul", "Mumford & Sons", "folk",
  "Back in My Body", "Maggie Rogers", "pop/rock",
  "Bodak Yellow", "Cardi B, Kodak Black", "hip-hop/rap",
  "Different Kind of World", "Maggie Rogers", "pop/rock",
  "Disco Man", "Remi Wolf", "pop/rock",
  "Don't Forget Me", "Maggie Rogers", "pop/rock",
  "Fake Empire", "The National", "pop/rock",
  "Fightin' Words", "Trace Adkins", "country",
  "Free Now (Apple Music Live)", "Gracie Abrams", "pop",
  "I Love You, I'm Sorry", "Gracie Abrams", "pop/rock",
  "I Need My Girl", "The National", "pop/rock",
  "I Will Wait", "Mumford & Sons", "folk",
  "Light On", "Maggie Rogers", "pop/rock",
  "Little Lion Man (City Sessions / Amazon Music Live)", "Mumford & Sons", "folk",
  "Love You For A Long Time", "Maggie Rogers", "pop/rock",
  "Only in America", "Brooks & Dunn", "country",
  "Piece of My Heart", "Beverley Knight", "r*b",
  "Risk", "Gracie Abrams", "pop",
  "Rumor", "Mike Parr, Andrew Shubin", "pop/rock",
  "September", "Earth, Wind & Fire", "r&b",
  "Something More", "Sugarland", "country",
  "The Bloodbuzz Ohio", "The National", "pop/rock",
  "The Cave (Apple Music Sessions)", "Mumford & Sons", "folk",
  "The Promised Land", "Bruce Springsteen", "pop/rock",
  "The Star-Spangled Banner - National Anthem USA", "United States of America", "patriotic",
  "Turn Down for What", "DJ Snake & Lil Jon", "electronic",
  "You Should Be Here", "Cole Swindell", "country",
)

## ---- apply ----------------------------------------------------------------
cand_songs <- cand_songs_cleanedf |>
  mutate(title_key = clean_title(title)) |>
  left_join(genre_fill, by = c("title", "performer")) |>
  left_join(lang_fix,   by = "title") |>
  mutate(
    genre     = coalesce(genre, genre_new),
    genre     = if_else(genre == "rap", "hip-hop/rap", genre),
    # CHANGED: lang_manual first, so a manual code overrides a confident-but-
    # wrong cld3 guess. The old order only filled NA / "too_short".
    lang_song = coalesce(lang_manual, na_if(lang_song, "too_short"))
  ) |>
  select(-genre_new, -lang_manual)

## ---- checks ---------------------------------------------------------------
cand_songs |> filter(is.na(genre))     |> distinct(title, performer)
cand_songs |> filter(is.na(lang_song)) |> distinct(title, performer)
count(cand_songs, genre,     sort = TRUE)
count(cand_songs, lang_song, sort = TRUE)

## 7b. did any manual language code actually override cld3?
cand_songs_cleanedf |>
  select(title, performer, lang_auto = lang_song) |>
  left_join(lang_fix, by = "title") |>
  filter(!is.na(lang_manual), !is.na(lang_auto),
         lang_auto != "too_short", lang_auto != lang_manual) |>
  distinct() |>
  print(n = Inf)

## 7c. are any lang_fix titles ambiguous (same title, >1 performer)?
cand_songs |>
  filter(title %in% lang_fix$title) |>
  distinct(title, performer) |>
  count(title) |>
  filter(n > 1) |>
  print(n = Inf)

## how many got dumped:
cand_songs |> 
  # play counts, on the full play-level data
  add_count(cand_year, title, performer, name = "num_plays") |>
  # campaign size = distinct songs
  group_by(cand_year) |> filter(n_words >= 20 ) |>
  mutate(song_ct = n_distinct(title, performer)) |>
  ungroup() |>
  filter(song_ct > 5) |>
  # collapse to one row per song  <- without this, slice_max picks 15 rows of
  # the same song, since num_plays is constant within cand_year x title
  distinct(cand_year, title, performer, .keep_all = TRUE) |>
  group_by(cand_year) |>
  slice_max(order_by = num_plays, n = 15, with_ties = FALSE) |>
  ungroup() |> count(lang_song)

########################### Dataset setup ##############
## MOVED HERE from the bottom of the script.
## English first, THEN top 15 — so every campaign contributes up to 15 English
## songs. Selecting first and filtering after would let playlist composition
## (correlated with party) decide N per campaign, focus on most played English songs.

cand_songs_sel <- cand_songs |>
  filter(lang_song == "en") |>
  # play counts, on the full play-level data using english-language songs
  add_count(cand_year, title, performer, name = "num_plays") |>
  # campaign size = distinct songs
  group_by(cand_year) |>
  mutate(song_ct = n_distinct(title, performer)) |>
  ungroup() |>
  filter(song_ct > 5) |>
  # collapse to one row per song 
  # the same song, since num_plays is constant within cand_year x title
  distinct(cand_year, title, performer, .keep_all = TRUE) |>
  group_by(cand_year) |>
  slice_sample(prop = 1) |> # Shuffles row order randomly
  # slice up to 15 songs per cand
  slice_max(order_by = num_plays, n = 15, with_ties = FALSE)|>
  ungroup() |>
  mutate(id = row_number())

## ---- verify the sample ----------------------------------------------------
cand_songs_sel |> count(cand_year) |> print(n = Inf)          # <= 15 each
stopifnot(!any(duplicated(cand_songs_sel$id)))
stopifnot(nrow(count(cand_songs_sel, cand_year, title, performer) |>
                 filter(n > 1)) == 0)

## what did English-only cost, and did any campaign drop out?
count(cand_songs, lang_song, useNA = "ifany") |> print(n = Inf)
setdiff(levels(droplevels(cand_songs$cand_year)),
        levels(droplevels(cand_songs_sel$cand_year)))


## Work with stop words, tokenization
# edit stop_words:
stop_words_filtered <- tidytext::stop_words

#
# ## lyrics: use get_tokens from syuzhet package
# cand_songs_tokenized <- cand_songs_sel |>
#   select(id, lyrics_clean) |>
#   filter(!is.na(lyrics_clean)) |>
#   #distinct(row_id, .keep_all = TRUE) |>
#   mutate(word = map(lyrics_clean, get_tokens)) |>
#   unnest(word) |>
#   anti_join(stop_words_filtered, by = "word") |>
#   group_by(id) |>
#   summarise(lyrics_tok = paste(word, collapse = " "), .groups = "drop")


## bring in NRC data # word: get sentiments, also syuzhet package
## produce full dataset w/NRC counts
# VERIFY XXX
lyr <- cand_songs_sel |>
  select(id, lyrics_clean) |>
  filter(!is.na(lyrics_clean))       # id is unique by construction now

sent <- get_nrc_sentiment(lyr$lyrics_clean)
names(sent) <- paste0("NRC_", names(sent))
sent$id <- lyr$id                      # safe: same object, same order

cand_songs_nrc <- left_join(cand_songs_sel, sent, by = "id")
stopifnot(nrow(cand_songs_nrc) == nrow(cand_songs_sel))

## NOTE: NRC counts are raw and scale with lyric length — rescale as
## 100 * x / WC once the LIWC word count comes back (see bottom).

## get NRC sD:
# Tokenize text and map it to VAD scores
# 1. normalize the lexicon's column names first
#    textdata::lexicon_nrc_vad() returns word / Valence / Arousal / Dominance
#    (v1, 20,007 words, scores on 0-1)
## NRC VAD v2.1 (Mohammad 2025): 55k terms, scores on -1..1
vad_file <- list.files("data", pattern = "NRC-VAD-Lexicon-v2.*\\.txt$",
                       recursive = TRUE, full.names = TRUE)[1]
vad_file                                   # check it found the right file

nrc_vad <- read_tsv(vad_file, show_col_types = FALSE)
names(nrc_vad)                             # expect term / valence / arousal / dominance
range(nrc_vad[[4]])                        # expect about -1 to 1

nrc_vad_std <- nrc_vad |>
  rename_with(tolower) |>
  rename(word = 1) |>                      # first column = the term
  select(word, valence, arousal, dominance) |>
  mutate(word = tolower(word)) |>
  filter(!str_detect(word, "\\s")) |>      # single words only (tokens are single words)
  mutate(across(c(valence, arousal, dominance), ~ (.x + 1) / 2))  # rescale to 0-1


# 2. tokenize the data frame, keeping id
vad_tokens <- cand_songs_nrc |>
  select(id, lyrics_clean) |>
  filter(!is.na(lyrics_clean)) |>
  unnest_tokens(word, lyrics_clean) |>
  left_join(nrc_vad_std, by = "word")     # left_join keeps zero-match songs visible

# 3. average per document, and keep coverage
summary_scores <- vad_tokens |>
  group_by(id) |>
  summarise(
    n_tokens      = n(),
    n_vad         = sum(!is.na(valence)),
    vad_coverage  = n_vad / n_tokens,
    avg_valence   = mean(valence,   na.rm = TRUE),
    avg_arousal   = mean(arousal,   na.rm = TRUE),
    avg_dominance = mean(dominance, na.rm = TRUE),
    n_vad_unique  = n_distinct(word[!is.na(valence)]),
    .groups = "drop"
  )

# low coverage = unreliable mean
summary_scores |> arrange(vad_coverage) |> print(n = 25)

cand_songs_final <- left_join(cand_songs_nrc, summary_scores, by = "id")


######## export, do liwc, re-import

## ONE text column, newlines collapsed, so the wrong column can't be picked in
## the LIWC import dialog and the CSV parses cleanly. Rejoin full lyrics on id.
cand_songs_final |>
  mutate(lyrics_liwc = str_squish(str_replace_all(lyrics_clean, "[\\r\\n]+", " "))) |>
  select(-lyrics, -lyrics_clean, -lyrics_detect) |>
  group_by(title, performer) |>
  mutate(performer = recode(performer, "Brooks and Dunn" = "Brooks & Dunn"),
         genre = case_when(
           title == "United States of America's national anthem" ~ "patriotic",
           title == "This Land Is Your Land" ~ "folk",
           title == "Higher Love" ~ "pop/rock",    
           TRUE ~ genre)) |>
    mutate(canon_row = { g <- which(str_detect(url_used, "genius\\.com")); if (length(g)) g[1] else 1L },
         lyrics_liwc    = lyrics_liwc[canon_row],
         url_used  = url_used[canon_row]) |>
  write_csv("data/en_songs.csv")

## >>> run LIWC-22 on data/en_songs.csv, text column = lyrics_liwc,
## >>> save as data/en_songs_liwc.csv, then continue.

en <- read_csv("data/en_songs_liwc.csv")

stopifnot(nrow(en) == nrow(cand_songs_final))
stopifnot(!any(duplicated(en$id)))
stopifnot(setequal(en$id, cand_songs_final$id))


# ---------------------------------------------------------------------------
# SPANISH BRANCH — parked (English-only decision). Everything below is intact;
# uncomment the block plus the esp write/read to bring Spanish back in.
# ---------------------------------------------------------------------------
# ## a litle convoluted but getting translation for spanish songs ###
# # write out spanish only songs
# cand_songs_final |> filter(lang_song!= "instrumental" & lang_song == "es") |> write_csv("data/esp_songs.csv")
#
# # need to use the above export to get the esp songs
# esp <- read_csv("data/esp_songs_liwc.csv")
#
# #claude used
# # prep for merge:
# meta <- c("title","performer","genre","sentiment","party","cand_full","cand_f","cand_l",
#           "event_title","campaign","event_type","lyrics_exists","lyrics","url_used",
#           "id","is_instrumental","lyrics_clean","lang_song")
#
# # new_name = old_name
# esp_map <- c(
#   Segment = "Segment_1", WC = "WC_1", WPS = "WPS_1", Dic = "Dic_1",
#   AllPunc = "AllPunc_1", Period = "Period_1", Comma = "Comma_1", QMark = "QMark_1",
#   Exclam = "Exclam_1", Apostro = "Apostro_1", OtherP = "OtherP_1", Emoji = "Emoji_1",
#   BigWords = "Sixltr", `function` = "Funct", pronoun = "TotPron", ppron = "PronPer",
#   i = "Yo", we = "Nosotro", you = "TuUtd", shehe = "ElElla", they = "Ellos",
#   ipron = "PronImp", article = "Articulo", verb = "Verbos", auxverb = "VerbAux",
#   adverb = "Adverb", prep = "Prepos", conj = "Conjunc", negate = "Negacio",
#   quantity = "Cuantif", number = "Numeros", swear = "Maldec",
#   focuspast = "Pasado", focuspresent = "Present", focusfuture = "Futuro",
#   Affect = "Afect", tone_pos = "EmoPos", tone_neg = "EmoNeg",
#   emo_anx = "Ansiedad", emo_anger = "Enfado", emo_sad = "Triste",
#   cogproc = "MecCog", insight = "Insight", cause = "Causa", discrep = "Discrep",
#   tentat = "Tentat", certitude = "Certeza",
#   Social = "Social_1", family = "Familia", friend = "Amigos",
#   Perception = "Percept", visual = "Ver", auditory = "Oir", feeling = "Sentir",
#   health = "Salud", sexual = "Sexual", food = "Ingerir",
#   motion = "Movim", space = "Espacio", time = "Tiempo",
#   work = "Trabajo", achieve = "Logro", leisure = "Placer", home = "Hogar",
#   money = "Dinero", relig = "Relig", death = "Muerte",
#   assent = "Asentir", nonflu = "NoFluen", filler = "Relleno"
# )
#
# # Spanish-only categories with no LIWC-22 counterpart — keep, prefixed
# esp_only <- c("Humanos","Biolog","Cuerpo","Relativ","Incl","Excl","Inhib",
#               "Subjuntiv","formal","informal","VosUtds",
#               "verbYO","verbTU","verbNOS","verbosEL","verbELLOS","verbVos")
#
# # BEFORE USING: confirm the _1 columns are the SPANISH run, not the stale
# # English one carried in from cand_songs. Spanish-dictionary coverage on
# # Spanish text should be much higher than English-dictionary coverage:
# #   esp |> summarise(dic_plain = mean(Dic, na.rm=TRUE),
# #                    dic_sfx   = mean(Dic_1, na.rm=TRUE),
# #                    i_plain   = mean(i,  na.rm=TRUE),
# #                    yo_es     = mean(Yo, na.rm=TRUE))
# # And see what else the Spanish run produced: grep("_1$", names(esp), value=TRUE)
# #
# # all_of() not any_of() for required columns — any_of() drops typos silently.
# map_missing <- esp_map[!esp_map %in% names(esp)]
# if (length(map_missing)) message("Unmapped (dropped): ",
#   paste(names(map_missing), "<-", map_missing, collapse = "; "))
# esp_map_ok <- esp_map[esp_map %in% names(esp)]
#
# esp_clean <- esp |>
#   select(all_of(meta), all_of(esp_map_ok), any_of(esp_only)) |>
#   rename_with(~ paste0("es_", .x), any_of(esp_only)) |>
#   mutate(liwc_dict = "es")
#
# # pull in new spanish songs
# # en_clean <- en |> mutate(liwc_dict = "en")   # <- en needs the tag too
# # cand_songs_full <- bind_rows(en_clean, esp_clean)


########################### Analysis variables ##############

front = c( "Clinton 2016", "Biden 2020",
           "Harris 2024", "Trump 2016", "Trump 2020", "Trump 2024")
women = c("Elizabeth Warren", "Hillary Clinton", "Kamala Harris",
          "Kirsten Gillibrand", "Nikki Haley", "Amy Klobuchar")

cand_songs_final_export <- en |>
  # NOTE: no second add_count() here — it used to overwrite num_plays, which
  # is the variable that drove the top-15 selection. It survives from above.
  mutate(
    cand_year     = fct_reorder(cand_year, year),   # CSV round-trip drops order
    woman         = if_else(cand_full %in% women, 1L, 0L),
    lastname_year = paste(cand_l, year),            # cand_l / year already exist
    front_runner  = if_else(lastname_year %in% front, 1L, 0L),
  ) |>
  mutate(performer = recode(performer, "Brooks and Dunn" = "Brooks & Dunn"),
         genre = case_when(
           title == "United States of America's National Anthem" ~ "patriotic",
           title == "This Land Is Your Land" ~ "folk",
           title == "Higher Love" ~ "pop/rock",    
           TRUE ~ genre))

#

# cand_songs_final_export <- cand_songs_final |> group_by(cand_full, campaign) |>
#   filter(!is_instrumental, (lang_song %in% "en" | lang_song %in% "es")) |>
#   mutate(song_ct = n()) |>
#   ungroup() |> filter(song_ct > 10)   # date no longer works for these events bc compressed


## ---- final sanity ---------------------------------------------------------
count(cand_songs_final_export, cand_year, front_runner) |> print(n = Inf)
count(cand_songs_final_export, woman)
cand_songs_final_export |>
  summarise(across(c(num_plays, WC, vad_coverage),
                   list(min = min, med = median, max = max), na.rm = TRUE)) |>
  glimpse()

write_csv(cand_songs_final_export, "data/cand_songs_bycampaign.csv")
