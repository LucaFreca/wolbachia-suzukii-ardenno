# =====================================================================
# PARTE 1 - PREVALENZA DI WOLBACHIA IN DROSOPHILA SUZUKII
#           + PREPARAZIONE DEI DATI METEO E CORRELAZIONI ESPLORATIVE
# =====================================================================
# Cosa fa questo script, in ordine:
#   1) carica pacchetti, stile grafico e funzioni di supporto
#   2) importa e pulisce i dati PCR e meteo
#   3) calcola la prevalenza (con IC 95%) e disegna i grafici A-H
#   4) esegue i test di Fisher (sesso, luogo)
#   5) costruisce le variabili meteo mensili e le unisce ai dati PCR
#   6) studia le correlazioni meteo / concentrazione DNA
#   7) salva tutto in "analisi_esplorativa.RData" (usato dalla Parte 2)
#
# NB: i nomi degli oggetti salvati (data_completo, meteo_mensile, ...)
#     NON vanno cambiati: la Parte 2 li usa così come sono.
# =====================================================================


# =====================================================================
# 1. PACCHETTI, CARTELLA DI LAVORO E STILE
# =====================================================================

library(tidyverse)   # dplyr, ggplot2, readr, tidyr, purrr, stringr
library(lubridate)   # gestione delle date
library(GGally)      # matrice di scatterplot (ggpairs)
library(corrplot)    # matrice di correlazione con clustering

# Cartella che contiene i file di dati (cambiarla solo qui)
cartella_lavoro <- "."   # cartella che contiene i file di dati (modificare se serve)
setwd(cartella_lavoro)

# Stile comune a TUTTI i grafici ggplot (non serve ripeterlo nei grafici)
theme_set(theme_classic(base_size = 13, base_family = "Arial"))

# Palette fissa: ogni gruppo/variabile ha sempre lo stesso colore
col_prevalenza  <- "#264653"  # blu petrolio  -> prevalenza
col_female      <- "#D55E00"  # vermiglio     -> Female
col_male        <- "#0072B2"  # blu           -> Male
col_bosco       <- "#009E73"  # verde         -> Bosco
col_tavola      <- "#E69F00"  # arancio       -> Tavola
col_temperatura <- "#D55E00"  # vermiglio     -> temperatura
col_umidita     <- "#0072B2"  # blu           -> umidità
col_vento       <- "#009E73"  # verde         -> vento
col_pressione   <- "#7B2CBF"  # viola         -> pressione

colori_sesso <- c("Female" = col_female, "Male" = col_male)
colori_luogo <- c("Bosco" = col_bosco, "Tavola" = col_tavola)


# =====================================================================
# 2. FUNZIONI DI SUPPORTO
# =====================================================================

# --- Statistiche "sicure": restituiscono NA se non ci sono dati validi
#     (evitano gli errori/warning di min/max/mean su vettori tutti NA)
min_sicuro   <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
max_sicuro   <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
media_sicura <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
mediana_sicura <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
somma_sicura <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
range_sicuro <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
sd_sicura    <- function(x) if (sum(!is.na(x)) < 2) NA_real_ else sd(x, na.rm = TRUE)

# --- Numero di osservazioni di un modello (GLM, glmmTMB, ecc.)
#     Usata anche nella Parte 2.
numero_osservazioni <- function(modello) {
  if (inherits(modello, "glmmTMB")) {
    n <- nrow(modello$frame)
    if (!is.null(n) && length(n) == 1L) return(n)
  }
  n <- tryCatch(stats::nobs(modello), error = function(e) NA_integer_)
  if (!is.na(n) && length(n) == 1L) return(as.integer(n))
  mf <- tryCatch(stats::model.frame(modello), error = function(e) NULL)
  if (!is.null(mf)) return(nrow(mf))
  stop("Impossibile determinare N per il modello: ",
       paste(class(modello), collapse = "/"))
}

# --- Prevalenza (%) e intervallo di confidenza 95% esatto (Clopper-Pearson)
#     gruppi = colonne per cui raggruppare (vuoto = prevalenza totale)
calcola_ic <- function(dati, gruppi = character()) {
  if (length(gruppi) > 0) {
    dati <- dati %>% group_by(across(all_of(gruppi)))
  }
  dati %>%
    summarise(
      total_n    = n(),
      positive_n = sum(wolb_presence == "yes"),
      prevalence = positive_n / total_n * 100,
      IC_lower   = binom.test(positive_n, total_n)$conf.int[1] * 100,
      IC_upper   = binom.test(positive_n, total_n)$conf.int[2] * 100,
      .groups = "drop"
    )
}

# --- Aggiunge data ed etichetta del mese, ordinate cronologicamente
#     (serve perché i mesi nei grafici non siano in ordine alfabetico)
aggiungi_mese_ordinato <- function(tabella) {
  tabella %>%
    mutate(
      month_date  = lubridate::my(month),
      month_label = format(month_date, "%b %Y")
    ) %>%
    arrange(month_date) %>%
    mutate(month_label = factor(month_label, levels = unique(month_label)))
}

# --- Pezzi di grafico riutilizzati
asse_prevalenza <- function() {
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))
}
etichette_x_inclinate <- function() {
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# --- Grafico di una serie meteo nel tempo (oraria o mensile)
grafico_serie_meteo <- function(dati, x, y, colore, titolo, etichetta_y,
                                con_punti = FALSE) {
  g <- ggplot(dati, aes(x = .data[[x]], y = .data[[y]])) +
    geom_line(color = colore, linewidth = if (con_punti) 0.8 else 0.7)
  if (con_punti) g <- g + geom_point(color = colore, size = 2.5)
  g +
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
    labs(title = titolo, x = "Mese", y = etichetta_y) +
    theme(axis.text.x = element_text(size = 9, hjust = 0.5))
}

# --- TRUE se una colonna numerica è utilizzabile in una correlazione:
#     almeno 2 valori validi, non costante, non tutta a zero.
variabile_correlabile <- function(x) {
  x_validi <- x[!is.na(x)]
  length(x_validi) >= 2 &&
    dplyr::n_distinct(x_validi) > 1 &&
    !all(x_validi == 0)
}


# =====================================================================
# 3. IMPORTAZIONE E PULIZIA DEI DATI PCR
# =====================================================================

data  <- read_csv("dsuz_prevalence_data.xlsx - PCR.csv")
meteo <- read_csv("talamona_hourly.csv")

data <- data %>%
  select(-((ncol(data) - 5):ncol(data))) %>%       # elimina le ultime 6 colonne
  filter(count_check == TRUE) %>%                  # solo righe controllate
  mutate(
    wolb_presence = tolower(trimws(wolb_presence)),  # "Yes " -> "yes"
    GADPH         = tolower(trimws(GADPH)),
    month         = paste(month, year, sep = "-")    # es. "May-2025"
  ) %>%
  filter(wolb_presence %in% c("yes", "no")) %>%    # solo PCR interpretabili
  mutate(
    month_date  = lubridate::my(month),             # mese come data (1° giorno)
    month_label = format(month_date, "%b %Y")
  )

# Data meteo in formato Date
meteo <- meteo %>%
  mutate(date = as.Date(date, format = "%d/%m/%Y"))


# =====================================================================
# 4. TABELLE DI PREVALENZA CON IC 95%
# =====================================================================

IC_totale         <- calcola_ic(data)                                  # A
IC_sesso          <- calcola_ic(data, "sex")                           # B
IC_mese           <- calcola_ic(data, "month")                         # C
IC_trappola       <- calcola_ic(data, "place")                         # D
IC_sesso_trappola <- calcola_ic(data, c("sex", "place"))               # E
IC_sesso_mese     <- calcola_ic(data, c("sex", "month"))               # F
IC_mese_trappola  <- calcola_ic(data, c("month", "place"))             # G
IC_completo       <- calcola_ic(data, c("sex", "month", "place"))      # H

# Ordine cronologico dei mesi per le tabelle temporali
IC_mese          <- aggiungi_mese_ordinato(IC_mese)
IC_sesso_mese    <- aggiungi_mese_ordinato(IC_sesso_mese)
IC_mese_trappola <- aggiungi_mese_ordinato(IC_mese_trappola)
IC_completo      <- aggiungi_mese_ordinato(IC_completo)


# =====================================================================
# 5. GRAFICI DI PREVALENZA (A-H)
# =====================================================================

# A) Prevalenza totale
grafico_A <- ggplot(IC_totale, aes(x = "Total", y = prevalence)) +
  geom_col(fill = col_prevalenza, width = 0.5) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.8) +
  asse_prevalenza() +
  labs(title = "Wolbachia prevalence", x = NULL, y = "Prevalence (%)")
print(grafico_A)

# B) Per sesso
grafico_B <- ggplot(IC_sesso, aes(x = sex, y = prevalence, fill = sex)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.8) +
  asse_prevalenza() +
  scale_fill_manual(values = colori_sesso) +
  labs(title = "Wolbachia prevalence by sex", x = "Sex", y = "Prevalence (%)") +
  theme(legend.position = "none")
print(grafico_B)

# C) Per mese
grafico_C <- ggplot(IC_mese, aes(x = month_label, y = prevalence, group = 1)) +
  geom_line(color = col_prevalenza, linewidth = 1) +
  geom_point(color = col_prevalenza, size = 3) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.7) +
  asse_prevalenza() +
  labs(title = "Wolbachia prevalence by month", x = "Month", y = "Prevalence (%)") +
  etichette_x_inclinate()
print(grafico_C)

# D) Per luogo
grafico_D <- ggplot(IC_trappola, aes(x = place, y = prevalence, fill = place)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.8) +
  asse_prevalenza() +
  scale_fill_manual(values = colori_luogo) +
  labs(title = "Wolbachia prevalence by sampling location",
       x = "Sampling location", y = "Prevalence (%)") +
  theme(legend.position = "none")
print(grafico_D)

# E) Per sesso e luogo
grafico_E <- ggplot(IC_sesso_trappola, aes(x = place, y = prevalence, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper),
                position = position_dodge(width = 0.8), width = 0.15, linewidth = 0.7) +
  asse_prevalenza() +
  scale_fill_manual(values = colori_sesso) +
  labs(title = "Wolbachia prevalence by sex and sampling location",
       x = "Sampling location", y = "Prevalence (%)", fill = "Sex")
print(grafico_E)

# F) Per sesso e mese
grafico_F <- ggplot(IC_sesso_mese,
                    aes(x = month_label, y = prevalence, color = sex, group = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.6) +
  asse_prevalenza() +
  scale_color_manual(values = colori_sesso) +
  labs(title = "Wolbachia prevalence by sex and month",
       x = "Month", y = "Prevalence (%)", color = "Sex") +
  etichette_x_inclinate()
print(grafico_F)

# G) Per luogo e mese
grafico_G <- ggplot(IC_mese_trappola,
                    aes(x = month_label, y = prevalence, color = place, group = place)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.6) +
  asse_prevalenza() +
  scale_color_manual(values = colori_luogo) +
  labs(title = "Wolbachia prevalence by sampling location and month",
       x = "Month", y = "Prevalence (%)", color = "Sampling location") +
  etichette_x_inclinate()
print(grafico_G)

# H) Per sesso, mese e luogo (un pannello per luogo)
grafico_H <- ggplot(IC_completo,
                    aes(x = month_label, y = prevalence, color = sex, group = sex)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper), width = 0.15, linewidth = 0.6) +
  facet_wrap(~ place) +
  asse_prevalenza() +
  scale_color_manual(values = colori_sesso) +
  labs(title = "Wolbachia prevalence by sex, month and sampling location",
       x = "Month", y = "Prevalence (%)", color = "Sex") +
  etichette_x_inclinate()
print(grafico_H)


# =====================================================================
# 6. TEST ESATTO DI FISHER
# =====================================================================
# Verifica se la prevalenza differisce tra i gruppi.
# Per ogni test: conteggi, percentuali di riga e risultato.

# A) Female vs Male
tab_sesso <- table(Sesso = data$sex, Risultato_PCR = data$wolb_presence)
print(tab_sesso)
percentuali_sesso <- prop.table(tab_sesso, margin = 1) * 100
print(round(percentuali_sesso, 2))
fisher_sesso <- fisher.test(tab_sesso)
print(fisher_sesso)

# B) Bosco vs Tavola
tab_luogo <- table(Luogo = data$place, Risultato_PCR = data$wolb_presence)
print(tab_luogo)
percentuali_luogo <- prop.table(tab_luogo, margin = 1) * 100
print(round(percentuali_luogo, 2))
fisher_luogo <- fisher.test(tab_luogo)
print(fisher_luogo)


# =====================================================================
# 7. SERIE TEMPORALI METEO (DATI ORARI)
# =====================================================================

grafico_temperatura_oraria <- grafico_serie_meteo(
  meteo, "date", "temp_C", col_temperatura,
  "Temperatura nel tempo - Talamona", "Temperatura (°C)")
print(grafico_temperatura_oraria)

grafico_umidita_oraria <- grafico_serie_meteo(
  meteo, "date", "rh_pct", col_umidita,
  "Umidità relativa nel tempo - Talamona", "Umidità relativa (%)")
print(grafico_umidita_oraria)

grafico_vento_orario <- grafico_serie_meteo(
  meteo, "date", "wind_kmh", col_vento,
  "Velocità del vento nel tempo - Talamona", "Velocità del vento (km/h)")
print(grafico_vento_orario)

grafico_pressione_oraria <- grafico_serie_meteo(
  meteo, "date", "press_hPa", col_pressione,
  "Pressione atmosferica nel tempo - Talamona", "Pressione atmosferica (hPa)")
print(grafico_pressione_oraria)


# =====================================================================
# 8. INDICI METEO DERIVATI DAI DATI GIORNALIERI
# =====================================================================
# Climatologia = media e SD delle medie giornaliere per mese dell'anno
# (tutti i gennaio insieme, tutti i febbraio insieme, ...).
# Con pochi anni di dati è da considerare provvisoria.

# Mese di riferimento (1° giorno del mese) e numero del mese
meteo <- meteo %>%
  mutate(
    anno_mese   = floor_date(date, unit = "month"),
    mese_numero = month(date)
  )

# A) Pioggia oraria: rain_col39 è un CUMULATO mensile, quindi si calcola
#    la differenza tra righe consecutive. Se il cumulato cala/si azzera
#    (reset), il valore corrente è già l'incremento.
meteo_precip_incrementi <- meteo %>%
  arrange(date) %>%
  group_by(anno_mese) %>%
  mutate(
    rain_increment_mm = case_when(
      is.na(rain_col39)                ~ NA_real_,
      is.na(lag(rain_col39))           ~ rain_col39,
      rain_col39 >= lag(rain_col39)    ~ rain_col39 - lag(rain_col39),
      TRUE                             ~ rain_col39
    )
  ) %>%
  ungroup()

# B) Valori giornalieri
meteo_giornaliero <- meteo_precip_incrementi %>%
  group_by(date, anno_mese, mese_numero) %>%
  summarise(
    temp_media_giornaliera        = media_sicura(temp_C),
    temp_range_giornaliero        = range_sicuro(temp_C),
    rh_media_giornaliera          = media_sicura(rh_pct),
    rh_range_giornaliero          = range_sicuro(rh_pct),
    precipitazione_giornaliera_mm = somma_sicura(rain_increment_mm),
    .groups = "drop"
  )

# C) Climatologia per mese dell'anno
climatologia_mensile <- meteo_giornaliero %>%
  group_by(mese_numero) %>%
  summarise(
    temp_climatologica_media = media_sicura(temp_media_giornaliera),
    temp_climatologica_sd    = sd_sicura(temp_media_giornaliera),
    rh_climatologica_media   = media_sicura(rh_media_giornaliera),
    rh_climatologica_sd      = sd_sicura(rh_media_giornaliera),
    .groups = "drop"
  )

# D) Soglie (media ± 1 SD) e giorni fuori soglia
#    (un confronto con NA conta come FALSE, cioè "non fuori soglia")
sopra <- function(x, soglia) !is.na(x) & !is.na(soglia) & x > soglia
sotto <- function(x, soglia) !is.na(x) & !is.na(soglia) & x < soglia

meteo_giornaliero_con_soglie <- meteo_giornaliero %>%
  left_join(climatologia_mensile, by = "mese_numero") %>%
  mutate(
    temp_soglia_superiore = temp_climatologica_media + temp_climatologica_sd,
    temp_soglia_inferiore = temp_climatologica_media - temp_climatologica_sd,
    rh_soglia_superiore   = rh_climatologica_media + rh_climatologica_sd,
    rh_soglia_inferiore   = rh_climatologica_media - rh_climatologica_sd,
    temp_sopra_soglia     = sopra(temp_media_giornaliera, temp_soglia_superiore),
    temp_sotto_soglia     = sotto(temp_media_giornaliera, temp_soglia_inferiore),
    rh_sopra_soglia       = sopra(rh_media_giornaliera, rh_soglia_superiore),
    rh_sotto_soglia       = sotto(rh_media_giornaliera, rh_soglia_inferiore),
    giorno_pioggia_gt_1mm = !is.na(precipitazione_giornaliera_mm) &
                            precipitazione_giornaliera_mm > 1
  )

# E) Riepilogo mensile degli indici
indici_meteo_mensili <- meteo_giornaliero_con_soglie %>%
  group_by(anno_mese) %>%
  summarise(
    # Temperatura
    range_giornaliero_medio_temp    = media_sicura(temp_range_giornaliero),
    sd_medie_giornaliere_temp       = sd_sicura(temp_media_giornaliera),
    giorni_temp_sopra_clim_plus_sd  = sum(temp_sopra_soglia),
    giorni_temp_sotto_clim_minus_sd = sum(temp_sotto_soglia),
    # Umidità relativa
    range_giornaliero_medio_rh      = media_sicura(rh_range_giornaliero),
    sd_medie_giornaliere_rh         = sd_sicura(rh_media_giornaliera),
    giorni_rh_sopra_clim_plus_sd    = sum(rh_sopra_soglia),
    giorni_rh_sotto_clim_minus_sd   = sum(rh_sotto_soglia),
    # Precipitazioni
    giorni_pioggia_gt_1mm           = sum(giorno_pioggia_gt_1mm),
    .groups = "drop"
  )

print(indici_meteo_mensili)   # controllare prima di usarli


# =====================================================================
# 9. STATISTICHE METEO MENSILI (MIN / MAX / MEDIA) E GRAFICI
# =====================================================================
# Maggio 2025 e novembre 2026 sono esclusi SOLO dai grafici
# (copertura dati insufficiente); restano in tabelle e analisi.
mesi_esclusi_grafici <- as.Date(c("2025-05-01", "2026-11-01"))

# Temperatura: min e max dalle colonne temp_min_C / temp_max_C, media da temp_C
temperatura_mensile <- meteo %>%
  group_by(anno_mese) %>%
  summarise(
    temperatura_min   = min_sicuro(temp_min_C),
    temperatura_max   = max_sicuro(temp_max_C),
    temperatura_media = media_sicura(temp_C),
    .groups = "drop"
  )

umidita_mensile <- meteo %>%
  group_by(anno_mese) %>%
  summarise(
    umidita_min   = min_sicuro(rh_pct),
    umidita_max   = max_sicuro(rh_pct),
    umidita_media = media_sicura(rh_pct),
    .groups = "drop"
  )

vento_mensile <- meteo %>%
  group_by(anno_mese) %>%
  summarise(
    vento_min   = min_sicuro(wind_kmh),
    vento_max   = max_sicuro(wind_kmh),
    vento_media = media_sicura(wind_kmh),
    .groups = "drop"
  )

pressione_mensile <- meteo %>%
  group_by(anno_mese) %>%
  summarise(
    pressione_min   = min_sicuro(press_hPa),
    pressione_max   = max_sicuro(press_hPa),
    pressione_media = media_sicura(press_hPa),
    .groups = "drop"
  )

# Precipitazioni:
# - rain_col39 = cumulato mensile (mm) che si azzera a inizio mese:
#   il totale del mese è il suo MASSIMO (sommarlo sovrastimerebbe).
# - rain_rate_mmh = intensità di pioggia (mm/h): min, max e media.
precipitazioni_mensili <- meteo %>%
  group_by(anno_mese) %>%
  summarise(
    precipitazioni_totali     = max_sicuro(rain_col39),
    precipitazioni_min_rate   = min_sicuro(rain_rate_mmh),
    precipitazioni_max_rate   = max_sicuro(rain_rate_mmh),
    precipitazioni_media_rate = media_sicura(rain_rate_mmh),
    .groups = "drop"
  )

# Tabelle filtrate per i grafici
temperatura_mensile_grafico <- filter(temperatura_mensile, !anno_mese %in% mesi_esclusi_grafici)
umidita_mensile_grafico     <- filter(umidita_mensile,     !anno_mese %in% mesi_esclusi_grafici)
vento_mensile_grafico       <- filter(vento_mensile,       !anno_mese %in% mesi_esclusi_grafici)
pressione_mensile_grafico   <- filter(pressione_mensile,   !anno_mese %in% mesi_esclusi_grafici)

print(grafico_serie_meteo(
  temperatura_mensile_grafico, "anno_mese", "temperatura_media", col_temperatura,
  "Temperatura media mensile a Talamona", "Temperatura media (°C)", con_punti = TRUE))

print(grafico_serie_meteo(
  umidita_mensile_grafico, "anno_mese", "umidita_media", col_umidita,
  "Umidità media mensile a Talamona", "Umidità relativa media (%)", con_punti = TRUE))

print(grafico_serie_meteo(
  vento_mensile_grafico, "anno_mese", "vento_media", col_vento,
  "Velocità media mensile del vento a Talamona", "Velocità media del vento (km/h)",
  con_punti = TRUE))

print(grafico_serie_meteo(
  pressione_mensile_grafico, "anno_mese", "pressione_media", col_pressione,
  "Pressione atmosferica media mensile a Talamona", "Pressione media (hPa)",
  con_punti = TRUE))


# =====================================================================
# 10. TABELLA METEO MENSILE UNICA E UNIONE CON I DATI PCR
# =====================================================================

meteo_mensile <- temperatura_mensile %>%
  full_join(umidita_mensile,        by = "anno_mese") %>%
  full_join(vento_mensile,          by = "anno_mese") %>%
  full_join(pressione_mensile,      by = "anno_mese") %>%
  full_join(precipitazioni_mensili, by = "anno_mese") %>%
  full_join(indici_meteo_mensili,   by = "anno_mese")

# Ogni campione PCR riceve le variabili meteo del proprio mese
data_completo <- data %>%
  left_join(meteo_mensile, by = c("month_date" = "anno_mese"))

# Controllo: una riga per mese con le principali medie meteo
print(
  data_completo %>%
    distinct(month_date, temperatura_media, umidita_media, vento_media,
             pressione_media, precipitazioni_totali, precipitazioni_media_rate) %>%
    arrange(month_date)
)


# =====================================================================
# 10B. GRAFICO C2 - PREVALENZA MENSILE + GIORNI CALDI (DOPPIO ASSE Y)
# =====================================================================
# Come il grafico C, con in più la linea dei giorni del mese con
# temperatura sopra la media climatologica + 1 SD
# (giorni_temp_sopra_clim_plus_sd, la variabile meteo del modello M11).
#   - asse Y sinistro: prevalenza di Wolbachia (%), linea continua + IC 95%
#   - asse Y destro:   numero di giorni caldi, linea tratteggiata
# ggplot ha un solo asse "vero": i giorni vengono moltiplicati per un
# fattore per stare sulla scala 0-100, e l'asse destro li riporta
# ai valori originali (dividendo per lo stesso fattore).

# Prevalenza mensile + giorni caldi dello stesso mese
dati_grafico_C2 <- IC_mese %>%
  left_join(
    indici_meteo_mensili %>% select(anno_mese, giorni_temp_sopra_clim_plus_sd),
    by = c("month_date" = "anno_mese")
  ) %>%
  # Agosto escluso SOLO da questo grafico (resta in tabelle e analisi).
  # Per togliere solo un anno specifico usare ad esempio:
  # filter(month_date != as.Date("2025-08-01"))
  filter(month(month_date) != 8) %>%
  droplevels()

# Fattore di scala: il massimo dei giorni corrisponde a 100 sull'asse sinistro
max_giorni_C2 <- max(1, ceiling(max(dati_grafico_C2$giorni_temp_sopra_clim_plus_sd,
                                    na.rm = TRUE)))
fattore_C2 <- 100 / max_giorni_C2

grafico_C2 <- ggplot(dati_grafico_C2, aes(x = month_label)) +
  # Linea tratteggiata e punti: giorni caldi (riscalati)
  geom_line(
    aes(y = giorni_temp_sopra_clim_plus_sd * fattore_C2,
        color = "Days above mean + 1 SD", group = 1),
    linewidth = 1, linetype = "dashed", na.rm = TRUE
  ) +
  geom_point(
    aes(y = giorni_temp_sopra_clim_plus_sd * fattore_C2,
        color = "Days above mean + 1 SD"),
    size = 2.5, shape = 17, na.rm = TRUE
  ) +
  # Linea, punti e IC 95%: prevalenza (come nel grafico C)
  geom_line(aes(y = prevalence, color = "Wolbachia prevalence", group = 1),
            linewidth = 1) +
  geom_point(aes(y = prevalence, color = "Wolbachia prevalence"), size = 3) +
  geom_errorbar(aes(ymin = IC_lower, ymax = IC_upper),
                color = col_prevalenza, width = 0.15, linewidth = 0.7) +
  scale_y_continuous(
    name     = "Prevalence (%)",
    limits   = c(0, 100),
    breaks   = seq(0, 100, 20),
    sec.axis = sec_axis(~ . / fattore_C2,
                        name   = "Days above climatological mean + 1 SD",
                        breaks = scales::breaks_pretty())
  ) +
  scale_color_manual(
    values = c("Wolbachia prevalence"   = col_prevalenza,
               "Days above mean + 1 SD" = col_temperatura),
    breaks = c("Wolbachia prevalence", "Days above mean + 1 SD"),
    name   = NULL
  ) +
  labs(title = "Wolbachia prevalence and warm days by month", x = "Month") +
  etichette_x_inclinate() +
  theme(
    legend.position    = "top",
    axis.title.y.left  = element_text(color = col_prevalenza),
    axis.text.y.left   = element_text(color = col_prevalenza),
    axis.title.y.right = element_text(color = col_temperatura),
    axis.text.y.right  = element_text(color = col_temperatura)
  )
print(grafico_C2)

# Togliere il # per salvarlo in alta risoluzione:
# ggsave("grafico_C2_prevalenza_giorni_caldi.png", grafico_C2,
#        width = 9, height = 5.5, dpi = 300)


# =====================================================================
# 11. CORRELAZIONI METEO / CONCENTRAZIONE DNA (ESPLORATIVE)
# =====================================================================
# Il meteo è mensile, quindi anche la concentrazione DNA viene
# riassunta per mese (evita la pseudoreplicazione).
# Le correlazioni sono esplorative: non dimostrano causalità.

# A) Concentrazione in formato numerico (gestisce la virgola decimale)
data_concentrazione <- data %>%
  mutate(
    concentration_num = suppressWarnings(
      as.numeric(gsub(",", ".", as.character(concentration), fixed = TRUE))
    )
  )

cat("\nCampioni con concentration non convertibile in numero:",
    sum(is.na(data_concentrazione$concentration_num) &
          !is.na(data_concentrazione$concentration)), "\n")

# B) Concentrazione media/mediana per mese
concentrazione_mensile <- data_concentrazione %>%
  group_by(month_date) %>%
  summarise(
    concentration_media       = media_sicura(concentration_num),
    concentration_mediana     = mediana_sicura(concentration_num),
    n_campioni_concentrazione = sum(!is.na(concentration_num)),
    .groups = "drop"
  )

# C) Tabella mensile: DNA + tutte le variabili meteo
correlazioni_mensili <- concentrazione_mensile %>%
  left_join(meteo_mensile, by = c("month_date" = "anno_mese")) %>%
  arrange(month_date)

print(correlazioni_mensili)
if (interactive()) View(correlazioni_mensili)

# D) Variabili numeriche utilizzabili (escluse: costanti, tutte 0, tutte NA)
dati_numerici_correlazioni <- correlazioni_mensili %>% select(where(is.numeric))

variabili_correlabili <- names(dati_numerici_correlazioni)[
  vapply(dati_numerici_correlazioni, variabile_correlabile, logical(1))
]

dati_matrice_correlazioni <- dati_numerici_correlazioni %>%
  select(all_of(variabili_correlabili))

cat("\nVariabili incluse nella matrice di correlazione:\n")
print(names(dati_matrice_correlazioni))
cat("\nVariabili escluse (costanti, tutte zero o tutte NA):\n")
print(setdiff(names(dati_numerici_correlazioni), variabili_correlabili))

# E) Matrice di correlazione di Spearman
matrice_correlazione <- cor(dati_matrice_correlazioni,
                            method = "spearman",
                            use = "pairwise.complete.obs")
print(round(matrice_correlazione, 3))

# F) Corrplot con clustering gerarchico
#    Raggruppa visivamente le variabili simili (possibile collinearità).
#    Se la matrice contiene NA il clustering non è possibile:
#    in quel caso si usa l'ordine originale.
ordine_corrplot <- if (anyNA(matrice_correlazione)) "original" else "hclust"
if (ordine_corrplot == "original") {
  message("La matrice contiene NA: corrplot disegnato senza clustering.")
}

disegna_corrplot <- function() {
  corrplot::corrplot(
    matrice_correlazione,
    method = "color", type = "upper", order = ordine_corrplot,
    addCoef.col = "black", number.cex = 0.45,
    tl.col = "black", tl.srt = 45, tl.cex = 0.65,
    diag = FALSE,
    title = "Correlazioni di Spearman - clustering gerarchico",
    mar = c(0, 0, 3, 0)
  )
}

# Versione PDF grande (più leggibile) ...
pdf("correlazioni_spearman_clustering_V12.pdf", width = 16, height = 14)
disegna_corrplot()
dev.off()

# ... e versione nel pannello Plots di RStudio
disegna_corrplot()

# G) Heatmap ggplot con clustering gerarchico
#    Distanza = 1 - |rho|: vicine le variabili molto correlate
#    (sia positivamente sia negativamente). Gli NA valgono 0 solo
#    per il calcolo del clustering.
matrice_per_clustering <- matrice_correlazione
matrice_per_clustering[is.na(matrice_per_clustering)] <- 0

clustering_correlazioni <- hclust(as.dist(1 - abs(matrice_per_clustering)),
                                  method = "complete")
ordine_variabili <- colnames(matrice_correlazione)[clustering_correlazioni$order]

matrice_long <- as.data.frame(as.table(matrice_correlazione)) %>%
  rename(Variabile_1 = Var1, Variabile_2 = Var2, rho = Freq) %>%
  mutate(
    Variabile_1 = factor(as.character(Variabile_1), levels = ordine_variabili),
    # sull'asse y il primo livello sta in basso: si invertono i livelli
    Variabile_2 = factor(as.character(Variabile_2), levels = rev(ordine_variabili))
  )

grafico_heatmap_correlazioni <- ggplot(
  matrice_long, aes(x = Variabile_1, y = Variabile_2, fill = rho)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(rho), "", sprintf("%.2f", rho))), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  labs(title = "Correlazioni tra variabili meteorologiche e biologiche",
       subtitle = "Spearman su dati aggregati mensilmente - clustering gerarchico",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))
print(grafico_heatmap_correlazioni)

# H) Test di Spearman: ogni variabile meteo vs concentrazione media DNA
dati_meteo_numerici <- meteo_mensile %>% select(where(is.numeric))

variabili_meteo_test <- names(dati_meteo_numerici)[
  vapply(dati_meteo_numerici, variabile_correlabile, logical(1))
]

# Restituisce rho e p-value; NA se ci sono meno di 3 mesi o valori costanti
calcola_correlazione_dna <- function(dati, variabile_meteo) {
  coppie <- dati %>%
    select(all_of(variabile_meteo), concentration_media) %>%
    drop_na()
  n_mesi <- nrow(coppie)

  dati_insufficienti <- n_mesi < 3 ||
    n_distinct(coppie[[variabile_meteo]]) < 2 ||
    n_distinct(coppie$concentration_media) < 2

  if (dati_insufficienti) {
    return(tibble(variabile_meteo = variabile_meteo, n_mesi = n_mesi,
                  rho_spearman = NA_real_, p_value = NA_real_))
  }

  test <- cor.test(coppie[[variabile_meteo]], coppie$concentration_media,
                   method = "spearman", exact = FALSE)

  tibble(variabile_meteo = variabile_meteo, n_mesi = n_mesi,
         rho_spearman = unname(test$estimate), p_value = test$p.value)
}

risultati_correlazione_DNA <- map_dfr(
  variabili_meteo_test,
  ~ calcola_correlazione_dna(correlazioni_mensili, .x)
) %>%
  mutate(rho_spearman = round(rho_spearman, 3),
         p_value      = signif(p_value, 3))

print(risultati_correlazione_DNA, n = Inf)

# I) Scatterplot: ogni variabile meteo vs concentrazione DNA (con LOESS)
grafici_meteo_DNA <- list()

for (variabile in variabili_meteo_test) {

  dati_grafico <- correlazioni_mensili %>%
    select(month_date, concentration_media, all_of(variabile)) %>%
    drop_na()

  if (nrow(dati_grafico) < 3 ||
      n_distinct(dati_grafico[[variabile]]) < 2 ||
      n_distinct(dati_grafico$concentration_media) < 2) {
    message("Grafico saltato per ", variabile, ": dati insufficienti o costanti.")
    next
  }

  grafico <- ggplot(dati_grafico, aes(x = .data[[variabile]], y = concentration_media)) +
    geom_point(color = col_prevalenza, size = 2.5, alpha = 0.8) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                color = "#D55E00", fill = "#D55E00", alpha = 0.15) +
    labs(title = paste("Concentrazione DNA vs", variabile),
         x = variabile, y = "Concentrazione media DNA (ng/µL)") +
    theme_classic(base_size = 12, base_family = "Arial")

  grafici_meteo_DNA[[variabile]] <- grafico
  print(grafico)
}


# =====================================================================
# 12. MATRICE GGPAIRS (TUTTE LE VARIABILI NUMERICHE VALIDE)
# =====================================================================
# Diagonale: distribuzione; sotto: scatterplot; sopra: rho di Spearman.
# Si lavora su una COPIA con nomi abbreviati, così la tabella
# originale dati_matrice_correlazioni resta intatta.

dati_ggpairs <- dati_matrice_correlazioni
names(dati_ggpairs) <- names(dati_ggpairs) %>%
  str_replace_all("_", " ") %>%
  str_replace_all("temperatura", "Temp.") %>%
  str_replace_all("umidita", "Umid.") %>%
  str_replace_all("vento", "Vento") %>%
  str_replace_all("pressione", "Press.") %>%
  str_replace_all("precipitazioni", "Precip.") %>%
  str_replace_all("concentration", "DNA") %>%
  str_replace_all("giorni", "gg") %>%
  str_replace_all("media", "med.") %>%
  str_replace_all("minimo", "min.") %>%
  str_replace_all("massimo", "max.")

matrice_scatter_correlazioni <- GGally::ggpairs(
  dati_ggpairs,
  upper = list(continuous = GGally::wrap("cor", method = "spearman",
                                         stars = TRUE, digits = 2)),
  lower = list(continuous = GGally::wrap("points", alpha = 0.65, size = 1.2)),
  diag  = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)),
  title = "Matrice di correlazione: tutte le variabili meteo e DNA"
) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(strip.text = element_text(size = 6),
        axis.text  = element_text(size = 5),
        plot.title = element_text(size = 13, face = "bold"))

# Con molte variabili il grafico è molto pesante: meglio salvarlo su file.
# Togliere il # per esportarlo:
# ggsave("matrice_correlazioni_tutte_variabili.png",
#        plot = matrice_scatter_correlazioni, width = 20, height = 20, dpi = 300)


# =====================================================================
# 13. SALVATAGGIO PER LA PARTE 2
# =====================================================================

save.image(file.path(cartella_lavoro, "analisi_esplorativa.RData"))
cat("\nParte 1 completata: salvato analisi_esplorativa.RData\n")
