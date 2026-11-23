# =====================================================================
# PARTE 2 - WOLBACHIA E VARIABILI METEO: MODELLI BINOMIALI
# =====================================================================
# Richiede il file "analisi_esplorativa.RData" creato dalla Parte 1.
#
# Cosa fa questo script, in ordine:
#   1) correlazioni tra variabili meteo (multicollinearità)
#   2) screening delle variabili meteo con GLM (una alla volta)
#   3) confronto AIC/BIC dello screening
#   4) costruzione dei 16 modelli (M1-M16)
#   5) confronto AIC/BIC dei 16 modelli
#   6) diagnostica dei residui DHARMa per M3 e M11
#   7) temperatura, umidità e pressione aggiunte a 4 strutture (28 modelli)
#   8) MODELLO FINALE M11: summary e odds ratio con IC 95%
#
# Tutti i file prodotti (CSV, PNG, TXT, RData) vengono salvati
# nella cartella di lavoro.
# =====================================================================


# =====================================================================
# 0. PACCHETTI, CARTELLA E DATI
# =====================================================================

library(tidyverse)   # dplyr, ggplot2, tidyr, purrr, ...
library(lubridate)
library(glmmTMB)     # GLM e GLMM binomiali

if (!requireNamespace("DHARMa", quietly = TRUE)) install.packages("DHARMa")
library(DHARMa)      # diagnostica dei residui simulati

cartella_lavoro <- "."   # cartella che contiene i file di dati (modificare se serve)
setwd(cartella_lavoro)

# Carica tutti gli oggetti della Parte 1 (data_completo, meteo_mensile, ...)
load("analisi_esplorativa.RData")


# =====================================================================
# FUNZIONI DI SUPPORTO
# =====================================================================

# TRUE se una colonna è utilizzabile in una correlazione
# (almeno 2 valori validi, non costante, non tutta a zero)
variabile_correlabile <- function(x) {
  x_validi <- x[!is.na(x)]
  length(x_validi) >= 2 &&
    dplyr::n_distinct(x_validi) > 1 &&
    !all(x_validi == 0)
}

# Interrompe lo script se mancano colonne in data_completo
controlla_variabili <- function(variabili, step) {
  mancanti <- setdiff(variabili, names(data_completo))
  if (length(mancanti) > 0) {
    stop(step, " interrotto. Variabili mancanti: ", paste(mancanti, collapse = ", "))
  }
}

# Prepara il campione per i modelli:
# - tiene solo le colonne richieste
# - risposta come fattore no/yes ("no" = riferimento -> si modella P(yes))
# - sex, place e month come fattori
# - tiene solo le righe senza NA (tutti i modelli usano lo stesso campione)
prepara_dati_modello <- function(dati, variabili) {
  dati_out <- dati %>%
    select(all_of(variabili)) %>%
    mutate(
      wolb_presence = factor(tolower(trimws(as.character(wolb_presence))),
                             levels = c("no", "yes")),
      across(any_of(c("sex", "place", "month")), ~ factor(.x))
    ) %>%
    filter(complete.cases(.)) %>%
    droplevels()

  if (nrow(dati_out) == 0) {
    stop("Nessuna osservazione completa per: ", paste(variabili, collapse = ", "))
  }
  if (nlevels(dati_out$wolb_presence) != 2) {
    stop("La risposta wolb_presence deve avere entrambi i livelli no/yes.")
  }
  if ("month" %in% names(dati_out) && nlevels(dati_out$month) < 2) {
    stop("La variabile month deve contenere almeno due gruppi.")
  }
  dati_out
}

# TRUE se un modello glmmTMB è andato a convergenza
# (ottimizzatore OK e matrice Hessiana definita positiva)
convergenza_ok <- function(m) {
  isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess)
}

# Formula del modello come testo su una riga
formula_testo <- function(m) paste(deparse(formula(m)), collapse = "")

# Tabella AIC/BIC: baseline + modelli estesi, con differenze dal baseline.
# ΔAIC/ΔBIC negativi = modello esteso migliore del baseline.
crea_tabella_ic <- function(modello_base, modelli_estesi) {
  modelli <- c(list(Baseline = modello_base), modelli_estesi)
  numerosita <- vapply(modelli, function(m) as.numeric(nobs(m)), numeric(1))

  if (length(unique(numerosita)) != 1) {
    stop("Confronto interrotto: i modelli non usano lo stesso numero di osservazioni.")
  }

  tibble(
    modello = names(modelli),
    formula = vapply(modelli, formula_testo, character(1)),
    N       = numerosita,
    AIC     = vapply(modelli, function(m) as.numeric(AIC(m)), numeric(1)),
    BIC     = vapply(modelli, function(m) as.numeric(BIC(m)), numeric(1))
  ) %>%
    mutate(delta_AIC = AIC - AIC[1],
           delta_BIC = BIC - BIC[1])
}

# Salva un grafico base-R in PNG (il codice del grafico va tra { })
salva_png <- function(file, codice_grafico, larghezza = 1200, altezza = 900) {
  png(file, width = larghezza, height = altezza, res = 180)
  on.exit(dev.off())
  force(codice_grafico)
  invisible(file)
}


# =====================================================================
# 1. CORRELAZIONI TRA VARIABILI METEOROLOGICHE
# =====================================================================
# Scopo: trovare coppie di predittori meteo molto correlati
# (rischio di multicollinearità). Soglia esplorativa: |rho| >= 0.70.
# La scelta finale considera anche la rilevanza biologica.

# 1.1 Variabili meteo numeriche utilizzabili
dati_meteo_step1 <- meteo_mensile %>% select(where(is.numeric))
variabili_meteo_step1 <- names(dati_meteo_step1)[
  vapply(dati_meteo_step1, variabile_correlabile, logical(1))
]
dati_meteo_step1 <- dati_meteo_step1 %>% select(all_of(variabili_meteo_step1))

cat("\nVariabili meteorologiche incluse nello Step 1:\n")
print(names(dati_meteo_step1))

# 1.2 Matrice di Spearman
matrice_correlazione_meteo <- cor(dati_meteo_step1,
                                  method = "spearman",
                                  use = "pairwise.complete.obs")
write.csv(matrice_correlazione_meteo, "STEP1_matrice_correlazione_meteo.csv",
          row.names = TRUE)

# 1.3 Coppie con |rho| sopra soglia (solo triangolo superiore = no doppioni)
soglia_correlazione_step1 <- 0.70

indici_coppie_correlate <- which(
  abs(matrice_correlazione_meteo) >= soglia_correlazione_step1 &
    upper.tri(matrice_correlazione_meteo),
  arr.ind = TRUE
)

coppie_correlate_step1 <- tibble(
  variabile_1  = rownames(matrice_correlazione_meteo)[indici_coppie_correlate[, 1]],
  variabile_2  = colnames(matrice_correlazione_meteo)[indici_coppie_correlate[, 2]],
  rho_spearman = matrice_correlazione_meteo[indici_coppie_correlate]
) %>%
  mutate(
    rho_assoluto = abs(rho_spearman),
    intensita = case_when(
      rho_assoluto >= 0.90 ~ "Molto elevata",
      rho_assoluto >= 0.80 ~ "Elevata",
      rho_assoluto >= 0.70 ~ "Moderata-alta",
      TRUE                 ~ "Inferiore alla soglia"
    )
  ) %>%
  arrange(desc(rho_assoluto))

cat("\nCoppie meteo con |Spearman rho| >=", soglia_correlazione_step1, ":\n")
print(coppie_correlate_step1, n = Inf)
write.csv(coppie_correlate_step1, "STEP1_coppie_meteo_correlate.csv", row.names = FALSE)

# 1.4 Quante coppie "forti" ha ogni variabile
riepilogo_correlazioni_step1 <- tibble(
  variabile = c(coppie_correlate_step1$variabile_1, coppie_correlate_step1$variabile_2)
) %>%
  count(variabile, name = "n_coppie_correlate") %>%
  arrange(desc(n_coppie_correlate))

cat("\nVariabili coinvolte in correlazioni elevate:\n")
print(riepilogo_correlazioni_step1, n = Inf)
write.csv(riepilogo_correlazioni_step1, "STEP1_riepilogo_variabili_correlate.csv",
          row.names = FALSE)

cat("\nSTEP 1 completato: salvati 3 file CSV.\n")


# =====================================================================
# 2. SCREENING DELLE VARIABILI METEO CON GLM BINOMIALI
# =====================================================================
# Tre serie di modelli: A) sex   B) place   C) sex + place
# In ogni serie si aggiunge UN predittore meteo alla volta e si
# confronta con il baseline della stessa serie (stesso campione).
# Analisi esplorativa, non causale.

# 2.1 Predittori meteo candidati
variabili_step2 <- c(
  "temperatura_media",
  "precipitazioni_media_rate",
  "umidita_media",
  "range_giornaliero_medio_rh",
  "pressione_min",
  "pressione_max",
  "pressione_media",
  "vento_min",
  "vento_max",
  "vento_media",
  "giorni_temp_sopra_clim_plus_sd",
  "giorni_temp_sotto_clim_minus_sd",
  "giorni_rh_sopra_clim_plus_sd",
  "giorni_rh_sotto_clim_minus_sd"
)

controlla_variabili(c("wolb_presence", "sex", "place", variabili_step2), "STEP 2")

# 2.2 Funzione che stima baseline + un modello per ogni variabile meteo
stima_serie_glm <- function(dati, formula_base) {
  base <- glm(as.formula(paste("wolb_presence ~", formula_base)),
              data = dati, family = binomial(link = "logit"))

  estesi <- setNames(
    lapply(variabili_step2, function(v) {
      glm(as.formula(paste("wolb_presence ~", formula_base, "+", v)),
          data = dati, family = binomial(link = "logit"))
    }),
    variabili_step2
  )

  list(base = base, estesi = estesi)
}

# 2.3 SERIE A - baseline: sex
dati_step2_sex <- prepara_dati_modello(
  data_completo, c("wolb_presence", "sex", variabili_step2))

cat("\n================ STEP 2A: BASELINE SEX ================\n")
cat("N usato:", nrow(dati_step2_sex), "\n")
print(table(dati_step2_sex$wolb_presence))
print(table(dati_step2_sex$sex))

serie_sex          <- stima_serie_glm(dati_step2_sex, "sex")
modello_base_sex   <- serie_sex$base
modelli_step2_sex  <- serie_sex$estesi

# 2.4 SERIE B - baseline: place
dati_step2_place <- prepara_dati_modello(
  data_completo, c("wolb_presence", "place", variabili_step2))

cat("\n================ STEP 2B: BASELINE PLACE ================\n")
cat("N usato:", nrow(dati_step2_place), "\n")
print(table(dati_step2_place$wolb_presence))
print(table(dati_step2_place$place))

serie_place          <- stima_serie_glm(dati_step2_place, "place")
modello_base_place   <- serie_place$base
modelli_step2_place  <- serie_place$estesi

# 2.5 SERIE C - baseline: sex + place
dati_step2_sex_place <- prepara_dati_modello(
  data_completo, c("wolb_presence", "sex", "place", variabili_step2))

cat("\n================ STEP 2C: BASELINE SEX + PLACE ================\n")
cat("N usato:", nrow(dati_step2_sex_place), "\n")
print(table(dati_step2_sex_place$wolb_presence))
print(table(dati_step2_sex_place$sex))
print(table(dati_step2_sex_place$place))

serie_sex_place          <- stima_serie_glm(dati_step2_sex_place, "sex + place")
modello_base_sex_place   <- serie_sex_place$base
modelli_step2_sex_place  <- serie_sex_place$estesi


# =====================================================================
# 3. CONFRONTO AIC / BIC DELLO SCREENING
# =====================================================================

risultati_ic_sex <- crea_tabella_ic(modello_base_sex, modelli_step2_sex)
cat("\n================ CONFRONTO AIC/BIC - SEX ================\n")
print(risultati_ic_sex, n = Inf)
write.csv(risultati_ic_sex, "STEP2A_confronto_AIC_BIC_sex.csv", row.names = FALSE)

risultati_ic_place <- crea_tabella_ic(modello_base_place, modelli_step2_place)
cat("\n================ CONFRONTO AIC/BIC - PLACE ================\n")
print(risultati_ic_place, n = Inf)
write.csv(risultati_ic_place, "STEP2B_confronto_AIC_BIC_place.csv", row.names = FALSE)

risultati_ic_sex_place <- crea_tabella_ic(modello_base_sex_place, modelli_step2_sex_place)
cat("\n================ CONFRONTO AIC/BIC - SEX + PLACE ================\n")
print(risultati_ic_sex_place, n = Inf)
write.csv(risultati_ic_sex_place, "STEP2C_confronto_AIC_BIC_sex_place.csv", row.names = FALSE)


# =====================================================================
# 4. COSTRUZIONE DEI 16 MODELLI (M1-M16)
# =====================================================================
# 4 strutture x 4 configurazioni meteo = 16 modelli:
#
#                          nessun meteo  vento  giorni_temp  entrambe
#   sex + place                 M1        M2        M3          M4
#   sex * place                 M5        M6        M7          M8
#   sex + place + (1|month)     M9        M10       M11         M12
#   sex * place + (1|month)     M13       M14       M15         M16
#
# (1|month) = intercetta casuale per mese di campionamento.
# Tutti i modelli usano lo stesso campione.

# 4.1 Campione comune
variabili_step4 <- c("wolb_presence", "sex", "place", "month",
                     "vento_media", "giorni_temp_sopra_clim_plus_sd")
controlla_variabili(variabili_step4, "STEP 4")

dati_step4 <- prepara_dati_modello(data_completo, variabili_step4)

cat("\n================ STEP 4: CAMPIONE COMUNE ================\n")
cat("Numero di osservazioni:", nrow(dati_step4), "\n")
cat("Numero di mesi:", nlevels(dati_step4$month), "\n")
print(table(dati_step4$wolb_presence))
print(table(dati_step4$sex))
print(table(dati_step4$place))
print(table(dati_step4$month))

# 4.2 Strutture e configurazioni meteo
strutture_modello <- c(
  "sex + place",
  "sex * place",
  "sex + place + (1 | month)",
  "sex * place + (1 | month)"
)

configurazioni_meteo <- list(
  nessun_meteo       = character(0),
  vento              = "vento_media",
  giorni_temperatura = "giorni_temp_sopra_clim_plus_sd",
  entrambe           = c("vento_media", "giorni_temp_sopra_clim_plus_sd")
)

# 4.3 Stima (l'ordine dei cicli dà la numerazione della tabella sopra)
modelli_step4 <- list()
indice_modello <- 1

for (struttura in strutture_modello) {
  for (meteo_vars in configurazioni_meteo) {

    formula_modello <- as.formula(
      paste("wolb_presence ~", paste(c(struttura, meteo_vars), collapse = " + "))
    )
    nome_modello <- paste0("M", indice_modello)

    modelli_step4[[nome_modello]] <- glmmTMB(
      formula_modello, data = dati_step4, family = binomial(link = "logit")
    )

    cat("\nStimato", nome_modello, ":", deparse(formula_modello), "\n")
    indice_modello <- indice_modello + 1
  }
}

cat("\nSTEP 4 completato: stimati", length(modelli_step4), "modelli.\n")


# =====================================================================
# 5. CONFRONTO DEI 16 MODELLI CON AIC E BIC
# =====================================================================
# Differenze calcolate:
#   - rispetto al baseline della propria struttura (M1, M5, M9, M13)
#   - rispetto a M1
# Valori negativi = modello migliore del riferimento.

numerosita_modelli <- vapply(modelli_step4, function(m) as.numeric(nobs(m)), numeric(1))
if (length(unique(numerosita_modelli)) != 1) {
  stop("Confronto interrotto: i 16 modelli non usano lo stesso numero di osservazioni.")
}

risultati_ic_step4 <- tibble(
  modello     = names(modelli_step4),
  formula     = vapply(modelli_step4, formula_testo, character(1)),
  N           = numerosita_modelli,
  AIC         = vapply(modelli_step4, function(m) as.numeric(AIC(m)), numeric(1)),
  BIC         = vapply(modelli_step4, function(m) as.numeric(BIC(m)), numeric(1)),
  convergenza = vapply(modelli_step4, convergenza_ok, logical(1))
) %>%
  mutate(
    baseline           = rep(c("M1", "M5", "M9", "M13"), each = 4),
    delta_AIC_baseline = AIC - AIC[match(baseline, modello)],
    delta_BIC_baseline = BIC - BIC[match(baseline, modello)],
    delta_AIC_M1       = AIC - AIC[modello == "M1"],
    delta_BIC_M1       = BIC - BIC[modello == "M1"]
  ) %>%
  arrange(AIC)

cat("\n================ CONFRONTO DEI 16 MODELLI ================\n")
print(risultati_ic_step4, n = Inf)
write.csv(risultati_ic_step4, "STEP4_confronto_16_modelli_AIC_BIC.csv", row.names = FALSE)

# Controllo convergenza
modelli_non_convergenti <- filter(risultati_ic_step4, !convergenza)
if (nrow(modelli_non_convergenti) > 0) {
  cat("\nATTENZIONE: verificare i seguenti modelli:\n")
  print(modelli_non_convergenti)
} else {
  cat("\nTutti i 16 modelli hanno superato i controlli di convergenza.\n")
}

cat("\nSTEP 5 completato.\n")


# =====================================================================
# 6. DIAGNOSTICA DEI RESIDUI CON DHARMa (M3 e M11)
# =====================================================================
# Controlli: uniformità, dispersione, outlier, e pattern dei residui
# rispetto a sex, place, giorni_temp e mese.
# Il controllo per mese è importante soprattutto per M11 (1|month).

dir.create("DHARMa_grafici", showWarnings = FALSE)

modelli_diagnostica <- c("M3", "M11")

# 6.1 Simulazione dei residui (seed fisso = risultati riproducibili)
set.seed(123)
residui_dharma <- lapply(
  setNames(modelli_diagnostica, modelli_diagnostica),
  function(nome) simulateResiduals(fittedModel = modelli_step4[[nome]], n = 1000)
)
residui_M3  <- residui_dharma$M3
residui_M11 <- residui_dharma$M11

# Variabili rispetto a cui controllare i residui (nome nel file = colonna)
predittori_residui <- c(
  sex                = "sex",
  place              = "place",
  giorni_temperatura = "giorni_temp_sopra_clim_plus_sd",
  month              = "month"
)

# 6.2 Test e grafici per ciascun modello
for (nome in modelli_diagnostica) {

  residui <- residui_dharma[[nome]]

  cat("\n================ DIAGNOSTICA MODELLO", nome, "================\n")
  print(testUniformity(residui))
  print(testDispersion(residui))
  print(testOutliers(residui))

  # Grafico generale (QQ-plot + residui vs predetti)
  salva_png(file.path("DHARMa_grafici", paste0(nome, "_diagnostica_generale.png")),
            plot(residui), larghezza = 1600, altezza = 1200)

  # Residui rispetto a ciascun predittore
  for (etichetta in names(predittori_residui)) {
    colonna <- predittori_residui[[etichetta]]
    salva_png(file.path("DHARMa_grafici", paste0(nome, "_residui_", etichetta, ".png")),
              plotResiduals(residui, form = dati_step4[[colonna]]))
  }
}

cat("\nDiagnostica DHARMa completata: grafici nella cartella DHARMa_grafici.\n")


# =====================================================================
# 7. TEMPERATURA, UMIDITÀ E PRESSIONE: INFORMAZIONE AGGIUNTIVA?
# =====================================================================
# Le 7 combinazioni di temperatura_media, umidita_media e
# pressione_media vengono aggiunte a 4 strutture di partenza:
#   A_BASE: sex + place
#   B_BASE: sex + place + (1 | month)
#   C_M3  : sex + place + giorni_temp_sopra_clim_plus_sd
#   D_M11 : sex + place + giorni_temp_sopra_clim_plus_sd + (1 | month)
# = 4 modelli base + 28 estesi, tutti sullo stesso campione.
# TUTTI i modelli (3 basi + 28 estesi = 31) sono confrontati con lo
# stesso riferimento A_BASE (sex + place) tramite AIC, BIC e test del
# rapporto di verosimiglianza (LRT).

# 7.1 Campione comune
variabili_step7 <- c("wolb_presence", "sex", "place", "month",
                     "giorni_temp_sopra_clim_plus_sd",
                     "temperatura_media", "umidita_media", "pressione_media")
controlla_variabili(variabili_step7, "STEP 7")

dati_step7 <- prepara_dati_modello(data_completo, variabili_step7)

cat("\n================ STEP 7: CAMPIONE COMUNE ================\n")
cat("Numero di osservazioni:", nrow(dati_step7), "\n")
cat("Numero di mesi:", nlevels(dati_step7$month), "\n")
print(table(dati_step7$wolb_presence))
print(table(dati_step7$sex))
print(table(dati_step7$place))

# 7.2 Correlazioni tra le tre variabili (possibile multicollinearità;
#     da sola la correlazione non prova un'interferenza nel modello)
matrice_correlazione_step7 <- cor(
  dati_step7 %>% select(temperatura_media, umidita_media, pressione_media),
  method = "spearman", use = "complete.obs"
)
cat("\nCorrelazioni di Spearman tra temperatura, umidità e pressione:\n")
print(round(matrice_correlazione_step7, 3))
write.csv(matrice_correlazione_step7,
          "STEP7_correlazioni_temperatura_umidita_pressione.csv", row.names = TRUE)

# 7.3 Le 7 combinazioni meteo
combinazioni_meteo_step7 <- list(
  TEMP          = "temperatura_media",
  RH            = "umidita_media",
  PRESS         = "pressione_media",
  TEMP_RH       = c("temperatura_media", "umidita_media"),
  TEMP_PRESS    = c("temperatura_media", "pressione_media"),
  RH_PRESS      = c("umidita_media", "pressione_media"),
  TEMP_RH_PRESS = c("temperatura_media", "umidita_media", "pressione_media")
)

# 7.4 Le 4 strutture di partenza, ristimate sul campione dello STEP 7
formule_base_step7 <- list(
  A_BASE = wolb_presence ~ sex + place,
  B_BASE = wolb_presence ~ sex + place + (1 | month),
  C_M3   = wolb_presence ~ sex + place + giorni_temp_sopra_clim_plus_sd,
  D_M11  = wolb_presence ~ sex + place + giorni_temp_sopra_clim_plus_sd + (1 | month)
)

modelli_base_step7 <- lapply(formule_base_step7, function(f) {
  glmmTMB(f, data = dati_step7, family = binomial(link = "logit"))
})

# 7.5 Stima dei 28 modelli estesi
#     Per ogni modello si memorizza anche la struttura di partenza.
modelli_estesi_step7 <- list()
formule_estese_step7 <- list()
struttura_estesi_step7 <- character()

for (nome_base in names(formule_base_step7)) {
  for (nome_combo in names(combinazioni_meteo_step7)) {

    formula_estesa <- update(
      formule_base_step7[[nome_base]],
      paste(". ~ . +", paste(combinazioni_meteo_step7[[nome_combo]], collapse = " + "))
    )
    nome_modello <- paste(nome_base, nome_combo, sep = "_")

    modelli_estesi_step7[[nome_modello]] <- glmmTMB(
      formula_estesa, data = dati_step7, family = binomial(link = "logit")
    )
    formule_estese_step7[[nome_modello]]   <- formula_estesa
    struttura_estesi_step7[[nome_modello]] <- nome_base

    cat("\nStimato", nome_modello, ":", deparse(formula_estesa), "\n")
  }
}

cat("\nSTEP 7: stimati", length(modelli_base_step7), "modelli base e",
    length(modelli_estesi_step7), "modelli estesi.\n")

# 7.6 Controllo: stesso campione per tutti
tutti_modelli_step7 <- c(modelli_base_step7, modelli_estesi_step7)
numerosita_step7 <- vapply(tutti_modelli_step7, function(m) as.numeric(nobs(m)), numeric(1))
if (length(unique(numerosita_step7)) != 1) {
  stop("Confronto STEP 7 interrotto: i modelli non usano lo stesso numero di osservazioni.")
}

# 7.7 Confronto di OGNI modello con lo STESSO riferimento: A_BASE (sex + place)
#     - Tutti i ΔAIC/ΔBIC sono quindi direttamente confrontabili tra loro.
#     - Tutti i modelli contengono sex + place, quindi sono annidati in
#       A_BASE e il test LRT è valido.
#     - Per i modelli con (1 | month) il LRT è conservativo (la varianza
#       casuale è testata al limite, cioè a 0): il p-value reale è
#       un po' più piccolo di quello riportato.
nome_riferimento_step7 <- "A_BASE"
modello_riferimento_step7 <- modelli_base_step7[[nome_riferimento_step7]]

confronta_con_base <- function(nome_modello, modello_esteso) {

  modello_base <- modello_riferimento_step7

  # Test del rapporto di verosimiglianza (modelli annidati)
  tabella_anova <- as.data.frame(anova(modello_base, modello_esteso))
  colonna_p <- grep("Pr\\(>Chisq\\)", names(tabella_anova), value = TRUE)
  p_lrt <- if (length(colonna_p) == 1) tabella_anova[[colonna_p]][2] else NA_real_

  tibble(
    struttura   = if (nome_modello %in% names(struttura_estesi_step7)) {
                    struttura_estesi_step7[[nome_modello]]
                  } else {
                    nome_modello
                  },
    confronto_con = nome_riferimento_step7,
    modello     = nome_modello,
    formula     = formula_testo(modello_esteso),
    N           = as.numeric(nobs(modello_esteso)),
    AIC_base    = as.numeric(AIC(modello_base)),
    AIC_esteso  = as.numeric(AIC(modello_esteso)),
    delta_AIC   = AIC_esteso - AIC_base,
    BIC_base    = as.numeric(BIC(modello_base)),
    BIC_esteso  = as.numeric(BIC(modello_esteso)),
    delta_BIC   = BIC_esteso - BIC_base,
    p_LRT       = as.numeric(p_lrt),
    convergenza = convergenza_ok(modello_esteso)
  )
}

# Si confrontano con A_BASE:
#   - le altre 3 strutture di partenza (B_BASE, C_M3, D_M11)
#   - i 28 modelli estesi
modelli_da_confrontare_step7 <- c(
  modelli_base_step7[names(modelli_base_step7) != nome_riferimento_step7],
  modelli_estesi_step7
)

risultati_step7 <- imap_dfr(
  modelli_da_confrontare_step7,
  function(modello, nome) confronta_con_base(nome, modello)
) %>%
  mutate(
    interpretazione_AIC = case_when(
      delta_AIC <= -2 ~ "AIC inferiore: possibile miglioramento",
      delta_AIC >=  2 ~ "AIC superiore: peggioramento",
      TRUE            ~ "Differenza AIC ridotta"
    ),
    interpretazione_BIC = case_when(
      delta_BIC <= -2 ~ "BIC inferiore: possibile miglioramento",
      delta_BIC >=  2 ~ "BIC superiore: peggioramento",
      TRUE            ~ "Differenza BIC ridotta"
    )
  ) %>%
  arrange(AIC_esteso)   # dal modello migliore (AIC più basso) al peggiore

cat("\n================ RISULTATI STEP 7 ================\n")
print(risultati_step7, n = Inf)
write.csv(risultati_step7, "STEP7_confronto_vs_sex_place_AIC_BIC_LRT.csv", row.names = FALSE)

# 7.8 Controllo convergenza
modelli_step7_problematici <- filter(risultati_step7, !convergenza)
if (nrow(modelli_step7_problematici) > 0) {
  cat("\nATTENZIONE: alcuni modelli richiedono verifiche:\n")
  print(modelli_step7_problematici, n = Inf)
} else {
  cat("\nTutti i modelli confrontati hanno superato i controlli di convergenza.\n")
}

# 7.9 Salvataggio
save(modelli_base_step7, modelli_estesi_step7, risultati_step7,
     matrice_correlazione_step7, file = "STEP7_modelli_e_risultati.RData")
cat("\nSTEP 7 completato.\n")


# =====================================================================
# 8. MODELLO FINALE M11: SUMMARY E ODDS RATIO
# =====================================================================
# M11 = sex + place + (1 | month) + giorni_temp_sopra_clim_plus_sd
# Modello scelto nel confronto AIC/BIC (STEP 5), stimato su dati_step4.
#
# Come leggere gli odds ratio (OR):
#   OR = exp(coefficiente logit)
#   OR > 1  -> odds di infezione da Wolbachia più alte
#   OR < 1  -> odds più basse
#   OR = 1  -> nessun effetto
#   Se l'IC 95% contiene 1, l'effetto non è significativo (p > 0.05).
# - sex e place: OR rispetto al livello di riferimento (stampato sotto)
# - giorni_temp_sopra_clim_plus_sd: OR per OGNI giorno in più al mese
#   con temperatura sopra la media climatologica + 1 SD
# - Intercetta: NON è un OR, sono gli odds di infezione per il gruppo
#   di riferimento con 0 giorni caldi
# IC 95% di Wald: exp(stima ± 1.96 × errore standard).

modello_M11 <- modelli_step4$M11

# 8.1 Livelli di riferimento (per interpretare gli OR)
cat("\n================ MODELLO FINALE M11 ================\n")
cat("Formula:", formula_testo(modello_M11), "\n")
cat("N osservazioni:", nobs(modello_M11), " - Mesi:", nlevels(dati_step4$month), "\n")
cat("Convergenza OK:", convergenza_ok(modello_M11), "\n")
cat("Riferimento sex:  ", levels(dati_step4$sex)[1], "\n")
cat("Riferimento place:", levels(dati_step4$place)[1], "\n")

# 8.2 Summary completo
summary_M11 <- summary(modello_M11)
print(summary_M11)

# 8.3 Varianza dell'intercetta casuale per mese
#     (quanto la prevalenza varia da un mese all'altro oltre ai predittori)
cat("\nEffetto casuale (1 | month):\n")
print(VarCorr(modello_M11))

# 8.4 Tabella degli odds ratio con IC 95%
coefficienti_M11 <- summary_M11$coefficients$cond
z_95 <- qnorm(0.975)

odds_ratio_M11 <- tibble(
  termine         = rownames(coefficienti_M11),
  stima_logit     = coefficienti_M11[, "Estimate"],
  errore_standard = coefficienti_M11[, "Std. Error"],
  z_value         = coefficienti_M11[, "z value"],
  p_value         = coefficienti_M11[, "Pr(>|z|)"],
  OR              = exp(stima_logit),
  IC95_inferiore  = exp(stima_logit - z_95 * errore_standard),
  IC95_superiore  = exp(stima_logit + z_95 * errore_standard),
  variazione_odds_percentuale = (OR - 1) * 100
) %>%
  mutate(
    tipo = if_else(termine == "(Intercept)",
                   "Odds di riferimento (intercetta, non è un OR)",
                   "Odds ratio"),
    significativo = if_else(p_value < 0.05, "Sì (p < 0.05)", "No")
  )

cat("\nOdds ratio del modello M11 (IC 95% di Wald):\n")
print(
  odds_ratio_M11 %>%
    mutate(across(c(stima_logit, errore_standard, z_value,
                    OR, IC95_inferiore, IC95_superiore,
                    variazione_odds_percentuale), ~ round(.x, 3)),
           p_value = signif(p_value, 3)),
  n = Inf, width = Inf
)

# 8.5 Salvataggio
write.csv(odds_ratio_M11, "STEP8_M11_odds_ratio.csv", row.names = FALSE)

writeLines(
  c(capture.output(print(summary_M11)),
    "",
    "Effetto casuale (1 | month):",
    capture.output(print(VarCorr(modello_M11)))),
  "STEP8_M11_summary.txt"
)

save(modello_M11, summary_M11, odds_ratio_M11, dati_step4,
     file = "STEP8_modello_M11.RData")

cat("\nSTEP 8 completato: salvati STEP8_M11_odds_ratio.csv,",
    "STEP8_M11_summary.txt e STEP8_modello_M11.RData\n")
