## =============================================================================
## PREPARAÇÃO DOS DADOS
## -----------------------------------------------------------------------------
## Lê a planilha de observações focais (BORIS), limpa valores ausentes nas
## colunas comportamentais, deriva os eixos comportamentais utilizados no
## modelo (ex.: Repouso, Atividade, Agitação, Vigilância) e salva os objetos
## `dados` e `dados_abertos` em .rds para uso em `modelo_unificado.R`.
##
## Este script reproduz apenas a etapa de preparação de dados da análise
## descritiva completa do TCC; os painéis gráficos e tabelas exploratórias
## (histogramas, mapas de correlação, pré-triagem de famílias de distribuição
## etc.) não fazem parte da reprodução dos testes estatísticos e por isso
## não foram incluídos aqui.
## =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

## -----------------------------------------------------------------------------
## 1) Leitura da planilha de dados
## -----------------------------------------------------------------------------
## Ajuste o caminho abaixo para o local da planilha `dados-observaçao-onca-pintada.xlsx`
## deste repositório (pasta `dados/`).
arquivo <- "CAMINHO/PARA/dados-observaçao-onca-pintada.xlsx"
if (!file.exists(arquivo)) stop("Arquivo não encontrado. Ajuste o caminho da variável 'arquivo'.")

dados <- read_excel(arquivo) %>%
  mutate(
    Data = as.Date(Data),
    Zoologico_aberto = factor(Zoologico_aberto, levels = c("Sim", "Nao"), labels = c("Aberto", "Fechado"))
  )

## -----------------------------------------------------------------------------
## 2) Definição dos grupos de colunas do etograma
## -----------------------------------------------------------------------------
DURACAO_SESSAO <- 600

cod_posturas <- c("REP", "SON", "LOC", "PAC", "OCU", "NEC", "AMBpar")
cod_sobre    <- c("AMB", "EXP", "HIG", "ENR", "VOC")
cod_visit    <- c("GRI", "BAT", "CHA", "AGL", "ASS", "PAL", "OBS", "Outro_visit", "CNV")

dur_posturas  <- paste0(cod_posturas, "_s")
freq_posturas <- paste0(cod_posturas, "_n")
dur_sobre     <- paste0(cod_sobre, "_s")
freq_sobre    <- paste0(cod_sobre, "_n")
dur_visit     <- paste0(cod_visit, "_s")
freq_visit    <- c(paste0(cod_visit, "_n"), "ATI_n", "NEN_n")

colunas_comportamento <- c(dur_posturas, freq_posturas, dur_sobre, freq_sobre, "NEM_n", dur_visit, freq_visit)

## -----------------------------------------------------------------------------
## 3) Ausentes nas colunas comportamentais viram zero (ausência de ocorrência)
## -----------------------------------------------------------------------------
dados <- dados %>%
  mutate(across(all_of(colunas_comportamento), ~ replace_na(.x, 0)))

## -----------------------------------------------------------------------------
## 4) Disponibilidade do número de visitantes
## -----------------------------------------------------------------------------
tem_visitantes <- any(!is.na(dados$N_visitantes))

## -----------------------------------------------------------------------------
## 5) Eixos comportamentais (variáveis de resposta do modelo) e preditores
## -----------------------------------------------------------------------------
dados <- dados %>%
  mutate(
    Repouso_s = REP_s, Sono_s = SON_s, Atividade_s = LOC_s, Agitacao_s = PAC_s,
    Vigilancia_s = AMB_s, Ocultacao_s = OCU_s, EmPe_s = AMBpar_s,
    Perturbacao_s = GRI_s + BAT_s + ASS_s + PAL_s,
    Presenca_calma_s = CNV_s + OBS_s
  ) %>%
  mutate(across(c(Repouso_s, Sono_s, Atividade_s, Agitacao_s, Vigilancia_s, Ocultacao_s, EmPe_s),
                ~ .x / DURACAO_SESSAO, .names = "{str_replace(.col, '_s$', '_prop')}"))

eixos_s    <- c("Repouso_s", "Sono_s", "Atividade_s", "Agitacao_s", "Vigilancia_s", "Ocultacao_s", "EmPe_s")
eixos_prop <- str_replace(eixos_s, "_s$", "_prop")
preditores <- c("Perturbacao_s", "Presenca_calma_s", "ATI_n", if (tem_visitantes) "N_visitantes")

## Subconjunto de sessões em dias com o zoológico aberto (usado nas análises
## que dependem da presença de visitantes)
dados_abertos <- dados %>% filter(Zoologico_aberto == "Aberto")

## -----------------------------------------------------------------------------
## 6) Salvar os objetos preparados
## -----------------------------------------------------------------------------
## Ajuste o caminho abaixo para a pasta de resultados onde `modelo_unificado.R`
## espera encontrar `dados_descritiva.rds`.
pasta_resultados <- "CAMINHO/PARA/pasta_de_resultados"
saveRDS(dados, file.path(pasta_resultados, "dados_descritiva.rds"))
saveRDS(dados_abertos, file.path(pasta_resultados, "dados_abertos_descritiva.rds"))
cat("Arquivos salvos em:", pasta_resultados, "\n")
