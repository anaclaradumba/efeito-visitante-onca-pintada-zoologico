# Modelo unificado: um modelo por desfecho comportamental (não um por pergunta),
# com todas as covariáveis de visitação entrando juntas (estimativas mutuamente
# ajustadas). Unidade de análise: sessão (10 min). Efeito aleatório: (1 | Data).
# Distribuição: Beta ordenada (ordbeta), que aceita 0 e 1 exatos. Correção de
# Holm por modelo. R² de Nakagawa (marginal/condicional). Agitação (pacing) é
# só descritiva (poucas ocorrências, todas em dias abertos).

pacotes <- c("dplyr", "glmmTMB", "DHARMa", "emmeans", "broom.mixed",
             "writexl", "ggplot2", "gridExtra", "tidyr", "scales",
             "performance")
faltando <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(faltando) > 0) install.packages(faltando)

library(dplyr)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(broom.mixed)
library(writexl)
library(ggplot2)
library(gridExtra)
library(performance)


# Seção 0A — paleta de cores e tema gráfico
paleta_tcc <- c(
  geral     = "#8c4327",
  Repouso   = "#8b5a2b",
  Sono      = "#d95a2b",
  Atividade = "#4a5d36",
  Ocultação = "#e6a135",
  Agitação  = "#c2b280"
)

tema_tcc <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x      = element_line(color = "grey40"),
    plot.title       = element_text(face = "bold"),
    plot.margin      = margin(0, 0, 0, 0)
  )

# Seção 0 — caminhos
pasta_resultados <- "CAMINHO/PARA/pasta_de_resultados"

cat("Pasta de resultados:", pasta_resultados, "\n")
cat("Arquivo dados_descritiva.rds encontrado?", file.exists(file.path(pasta_resultados, "dados_descritiva.rds")), "\n")


# Seção 1 — leitura dos dados preparados (ver 01_preparacao_dados.R)
dados <- readRDS(file.path(pasta_resultados, "dados_descritiva.rds"))


# Seção 2 — Tratamento (3 níveis: Fechado / moderado / intenso)
dados <- dados %>%
  mutate(
    Tratamento = case_when(
      Zoologico_aberto == "Fechado"                 ~ "Fechado",
      grepl("quarta",  tolower(Dia_semana))          ~ "Aberto (moderado)",
      grepl("domingo", tolower(Dia_semana))          ~ "Aberto (intenso)",
      TRUE                                           ~ NA_character_
    ),
    Tratamento = factor(Tratamento, levels = c("Fechado", "Aberto (moderado)", "Aberto (intenso)"))
  )

stopifnot(!anyNA(dados$Tratamento))

cat("\n--- Sessões por nível de Tratamento ---\n")
print(table(dados$Tratamento))


# Seção 3 — padronização (z-score) dos preditores contínuos
padronizar_aberto <- function(x, aberto) {
  valores_abertos <- x[aberto == "Aberto"]
  m <- mean(valores_abertos)
  s <- sd(valores_abertos)
  ifelse(aberto == "Aberto", (x - m) / s, 0)
}

dados <- dados %>%
  mutate(
    Perturbacao_z     = padronizar_aberto(Perturbacao_s, Zoologico_aberto),
    Presenca_calma_z  = padronizar_aberto(Presenca_calma_s, Zoologico_aberto),
    N_visitantes_z    = padronizar_aberto(N_visitantes, Zoologico_aberto)
  )


# Seção 4 — desfechos contínuos e fórmula unificada
eixos_continuos <- c("Repouso_prop", "Sono_prop", "Atividade_prop", "Ocultacao_prop")

nomes_eixos <- c(
  Repouso_prop   = "Repouso",
  Sono_prop      = "Sono",
  Atividade_prop = "Atividade",
  Ocultacao_prop = "Ocultação"
)

formula_unificada <- paste(
  "%s ~ Tratamento +",
  "Perturbacao_z + Presenca_calma_z + N_visitantes_z +",
  "(1 | Data)"
)


# Seção 5 — ajuste dos 4 desfechos (beta ordenada, um modelo por eixo)
resultados_continuos <- list()

for (eixo in eixos_continuos) {
  formula_txt <- sprintf(formula_unificada, eixo)
  f <- as.formula(formula_txt)

  modelo <- glmmTMB(f, data = dados, family = ordbeta())

  resultados_continuos[[eixo]] <- list(
    modelo       = modelo,
    nome_eixo    = nomes_eixos[[eixo]],
    var_resposta = eixo
  )

  cat(sprintf("%s: modelo ajustado (ordbeta, (1 | Data) mantido)\n", nomes_eixos[[eixo]]))
}


# Seção 6 — diagnóstico de resíduos (DHARMa)
diagnostico_ok <- function(modelo, n_sim = 1000) {
  sim   <- simulateResiduals(modelo, n = n_sim)
  teste <- testResiduals(sim, plot = FALSE)
  list(sim = sim, teste = teste)
}

for (eixo in eixos_continuos) {
  resultados_continuos[[eixo]]$diagnostico <- diagnostico_ok(resultados_continuos[[eixo]]$modelo)
}

diagnostico_residuos_tabela <- bind_rows(lapply(resultados_continuos, function(r) {
  data.frame(
    eixo                     = r$nome_eixo,
    familia                  = "Beta ordenada (ordbeta)",
    uniformidade_estatistica = unname(r$diagnostico$teste$uniformity$statistic),
    uniformidade_p           = r$diagnostico$teste$uniformity$p.value,
    dispersao_estatistica    = unname(r$diagnostico$teste$dispersion$statistic),
    dispersao_p              = r$diagnostico$teste$dispersion$p.value,
    outliers_p               = r$diagnostico$teste$outliers$p.value,
    prop_zeros               = mean(dados[[r$var_resposta]] == 0, na.rm = TRUE)
  )
}))

cat("\n--- Diagnóstico de resíduos (DHARMa), por desfecho — zero só é ressalva se dispersão/uniformidade falharem ---\n")
print(diagnostico_residuos_tabela)


# Seção 7 — extração de resultados e correção de Holm
extrair_resultado <- function(modelo, nome_eixo) {
  broom.mixed::tidy(modelo, effects = "fixed", conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      eixo   = nome_eixo,
      p_holm = p.adjust(p.value, method = "holm")
    )
}

resultados_tabela <- bind_rows(
  lapply(resultados_continuos, function(r) {
    extrair_resultado(r$modelo, r$nome_eixo)
  })
)

resultados_tabela <- resultados_tabela %>%
  left_join(
    diagnostico_residuos_tabela %>% select(eixo, dispersao_p, uniformidade_p),
    by = "eixo"
  )

cat("\n--- Resultados (desfechos contínuos), com p de Holm ---\n")
print(resultados_tabela)


# Seção 8 — salvar resultados (.rds dos modelos, tabelas em .xlsx)
saveRDS(resultados_continuos, file.path(pasta_resultados, "modelos_unificados_continuos.rds"))

write_xlsx(resultados_tabela,           file.path(pasta_resultados, "resultados_unificados_continuos.xlsx"))
write_xlsx(diagnostico_residuos_tabela, file.path(pasta_resultados, "resultados_diagnostico_residuos.xlsx"))


# Seção 9 — tabelas formatadas para leitura (e em imagem PNG)
rotulos_termo <- c(
  "TratamentoAberto (moderado)" = "Aberto, visitação moderada (vs. Fechado)",
  "TratamentoAberto (intenso)"  = "Aberto, visitação intensa (vs. Fechado)",
  "Perturbacao_z"               = "Perturbação",
  "Presenca_calma_z"            = "Presença calma",
  "N_visitantes_z"              = "Nº de visitantes"
)

formatar_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
}

tabela_continuos_exibicao <- resultados_tabela %>%
  transmute(
    `Eixo comportamental` = eixo,
    Termo              = rotulos_termo[term],
    Estimativa         = round(estimate, 3),
    `Erro padrão`      = round(std.error, 3),
    `IC 95% (inf.)`    = round(conf.low, 3),
    `IC 95% (sup.)`    = round(conf.high, 3),
    `p (Wald, individual)` = formatar_p(p.value),
    `p (Holm)`         = formatar_p(p_holm),
    `p dispersão (DHARMa)` = round(dispersao_p, 3)
  )

tabela_diagnostico_exibicao <- diagnostico_residuos_tabela %>%
  transmute(
    `Eixo comportamental` = eixo,
    `Família`            = familia,
    `% de zeros`         = round(prop_zeros * 100, 1),
    `Uniformidade (p)`   = round(uniformidade_p, 3),
    `Dispersão (p)`      = round(dispersao_p, 3),
    `Outliers (p)`       = round(outliers_p, 3)
  )

write_xlsx(tabela_continuos_exibicao,   file.path(pasta_resultados, "tabela_resultados_continuos_exibicao.xlsx"))
write_xlsx(tabela_diagnostico_exibicao, file.path(pasta_resultados, "tabela_diagnostico_exibicao.xlsx"))

salvar_tabela_png <- function(tabela, caminho_arquivo, cor_cabecalho = paleta_tcc[["geral"]],
                               largura = 11, altura_por_linha = 0.35) {
  grob_tabela <- tableGrob(
    tabela,
    rows = NULL,
    theme = ttheme_minimal(
      core    = list(fg_params = list(fontsize = 9), bg_params = list(fill = c("white", "#fbf6ee")),
                     padding = unit(c(3, 2.2), "mm")),
      colhead = list(fg_params = list(fontsize = 9, fontface = "bold", col = "white"),
                     bg_params = list(fill = cor_cabecalho),
                     padding = unit(c(3, 2.2), "mm"))
    )
  )
  largura_grob <- as.numeric(grid::convertWidth(sum(grob_tabela$widths), "in"))
  altura_grob  <- as.numeric(grid::convertHeight(sum(grob_tabela$heights), "in"))
  largura_final <- min(largura_grob, largura)
  ggsave(caminho_arquivo, grob_tabela, width = largura_final, height = altura_grob,
         units = "in", dpi = 600, bg = "white", limitsize = FALSE)
}

salvar_tabela_png(tabela_continuos_exibicao,   file.path(pasta_resultados, "tabela_resultados_continuos.png"),
                   cor_cabecalho = paleta_tcc[["geral"]], largura = 13)
salvar_tabela_png(tabela_diagnostico_exibicao, file.path(pasta_resultados, "tabela_diagnostico_residuos.png"),
                   cor_cabecalho = paleta_tcc[["geral"]], largura = 10)


# Seção 10 — forest plots em razão de chances (OR)
formatar_rotulo_or <- function(or_vals, low, high, p_vals) {
  fmt_num <- function(x) ifelse(abs(x) >= 100, sprintf("%.0f", x),
                          ifelse(abs(x) >= 10,  sprintf("%.1f", x),
                                                 sprintf("%.2f", x)))
  fmt_p <- function(p) ifelse(p >= 0.999, "p > 0,99",
                        ifelse(p < 0.001, "p < 0,001",
                               sprintf("p = %.3f", p)))
  sprintf(
    "OR = %s  [%s; %s]   %s",
    fmt_num(or_vals), fmt_num(low), fmt_num(high), fmt_p(p_vals)
  )
}

arquivo_forest_or <- c(
  Repouso   = "forest_or_repouso.png",
  Sono      = "forest_or_sono.png",
  Atividade = "forest_or_atividade.png",
  Ocultação = "forest_or_ocultacao.png"
)

desenhar_forest_or <- function(tabela_eixo, nome_eixo, caminho_arquivo) {
  cor_eixo <- paleta_tcc[["geral"]]

  tabela_plot <- tabela_eixo %>%
    mutate(
      Termo = rotulos_termo[term],
      OR         = exp(estimate),
      OR_low     = exp(conf.low),
      OR_high    = exp(conf.high),
      significativo = p_holm < 0.05,
      rotulo = formatar_rotulo_or(OR, OR_low, OR_high, p_holm)
    )

  grafico <- tabela_plot %>%
    ggplot(aes(x = OR, y = Termo)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = OR_low, xmax = OR_high, alpha = significativo),
                  orientation = "y", width = 0.15, color = cor_eixo, linewidth = 1) +
    geom_point(aes(alpha = significativo), size = 3.2, color = cor_eixo) +
    geom_label(
      aes(label = rotulo, fill = significativo),
      nudge_y = 0.32, size = 3.1, label.r = unit(0.15, "lines"),
      label.padding = unit(0.22, "lines"), color = ifelse(tabela_plot$significativo, "white", "grey20"),
      show.legend = FALSE
    ) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.55), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = cor_eixo, `FALSE` = "grey92"), guide = "none") +
    scale_x_log10() +
    labs(
      x = "Razão de Chances (OR) — escala logarítmica",
      y = NULL,
      title = sprintf("Efeito de cada preditor sobre %s", nome_eixo)
    ) +
    tema_tcc

  ggsave(caminho_arquivo, grafico, width = 8.6, height = 0.95 * nrow(tabela_eixo) + 1.9,
         dpi = 600, bg = "white", limitsize = FALSE)
}

for (eixo_nome in names(arquivo_forest_or)) {
  desenhar_forest_or(
    tabela_eixo     = resultados_tabela %>% filter(eixo == eixo_nome),
    nome_eixo       = eixo_nome,
    caminho_arquivo = file.path(pasta_resultados, arquivo_forest_or[[eixo_nome]])
  )
}

arquivo_tabela_or <- c(
  Repouso   = "tabela_repouso_or.png",
  Sono      = "tabela_sono_or.png",
  Atividade = "tabela_atividade_or.png",
  Ocultação = "tabela_ocultacao_or.png"
)
for (eixo_nome in names(arquivo_tabela_or)) {
  tabela_or_eixo <- resultados_tabela %>%
    filter(eixo == eixo_nome) %>%
    transmute(
      Termo            = rotulos_termo[term],
      OR               = round(exp(estimate), 2),
      `IC 95% (OR)`    = sprintf("[%s; %s]", round(exp(conf.low), 2), round(exp(conf.high), 2)),
      `p (Wald, individual)` = formatar_p(p.value),
      `p (Holm)`       = formatar_p(p_holm)
    )
  salvar_tabela_png(
    tabela_or_eixo,
    file.path(pasta_resultados, arquivo_tabela_or[[eixo_nome]]),
    cor_cabecalho = paleta_tcc[["geral"]], largura = 9.5
  )
}


# Seção 11 — médias marginais por Tratamento (emmeans), figura complementar
arquivo_emm_tratamento <- c(
  Repouso   = "emm_tratamento_repouso.png",
  Sono      = "emm_tratamento_sono.png",
  Atividade = "emm_tratamento_atividade.png",
  Ocultação = "emm_tratamento_ocultacao.png"
)
arquivo_contrastes_tratamento <- c(
  Repouso   = "contrastes_tratamento_repouso.xlsx",
  Sono      = "contrastes_tratamento_sono.xlsx",
  Atividade = "contrastes_tratamento_atividade.xlsx",
  Ocultação = "contrastes_tratamento_ocultacao.xlsx"
)

tabela_emm_todos <- list()

for (eixo in eixos_continuos) {
  r <- resultados_continuos[[eixo]]

  emm <- emmeans(r$modelo, ~ Tratamento, type = "response")
  emm_df <- as.data.frame(emm)

  nome_coluna_media <- intersect(c("response", "prob", "emmean"), names(emm_df))[1]

  grafico_emm <- emm_df %>%
    ggplot(aes(x = Tratamento, y = .data[[nome_coluna_media]])) +
    geom_point(size = 3.2, color = paleta_tcc[["geral"]]) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                  width = 0.15, color = paleta_tcc[["geral"]], linewidth = 1) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      x = NULL, y = "Proporção de tempo estimada (IC 95%)",
      title = sprintf("Média marginal estimada de %s, por Nível de Visitação", r$nome_eixo)
    ) +
    tema_tcc

  ggsave(file.path(pasta_resultados, arquivo_emm_tratamento[[r$nome_eixo]]),
         grafico_emm, width = 6.5, height = 4.8, dpi = 600, bg = "white")

  contrastes <- emmeans(r$modelo, ~ Tratamento) %>%
    regrid(transform = "response") %>%
    pairs(adjust = "holm")
  contrastes_df <- as.data.frame(contrastes)

  contrastes_exibicao <- contrastes_df %>%
    transmute(
      Contraste          = contrast,
      `Diferença (p.p.)` = round(estimate * 100, 2),
      `Erro padrão`      = round(SE * 100, 2),
      `p (Holm)`         = formatar_p(p.value)
    )

  write_xlsx(contrastes_exibicao, file.path(pasta_resultados, arquivo_contrastes_tratamento[[r$nome_eixo]]))

  tabela_emm_todos[[eixo]] <- emm_df %>% mutate(eixo = r$nome_eixo)

  cat(sprintf("%s: médias marginais e contrastes de Tratamento salvos\n", r$nome_eixo))
}


# Seção 12 — R² marginal e condicional (Nakagawa)
metricas_desempenho <- bind_rows(lapply(resultados_continuos, function(r) {
  r2 <- performance::r2_nakagawa(r$modelo, tolerance = 1e-10)

  data.frame(
    `Eixo comportamental` = r$nome_eixo,
    `R² marginal`    = round(as.numeric(r2$R2_marginal), 3),
    `R² condicional` = round(as.numeric(r2$R2_conditional), 3),
    check.names = FALSE
  )
}))

write_xlsx(metricas_desempenho, file.path(pasta_resultados, "tabela_metricas_desempenho.xlsx"))
salvar_tabela_png(metricas_desempenho, file.path(pasta_resultados, "tabela_metricas_desempenho.png"), largura = 8)

cat("\n--- Métricas de desempenho (R² de Nakagawa: marginal e condicional) ---\n")
print(metricas_desempenho)


# Seção 13 — Agitação (pacing): só descritivo, sem modelo estatístico
dados <- dados %>%
  mutate(PAC_ocorreu = as.integer(Agitacao_s > 0))

taxa_agitacao <- dados %>%
  group_by(Tratamento) %>%
  summarise(
    N_sessoes       = n(),
    N_com_agitacao  = sum(PAC_ocorreu),
    Proporcao       = N_com_agitacao / N_sessoes,
    .groups = "drop"
  )

taxa_agitacao_exibicao <- taxa_agitacao %>%
  rename(`Nível de Visitação` = Tratamento)

write_xlsx(taxa_agitacao_exibicao, file.path(pasta_resultados, "tabela_agitacao_descritiva.xlsx"))
salvar_tabela_png(taxa_agitacao_exibicao, file.path(pasta_resultados, "tabela_agitacao_descritiva.png"),
                   cor_cabecalho = paleta_tcc[["geral"]], largura = 8)

grafico_taxa_agitacao <- taxa_agitacao %>%
  ggplot(aes(x = Tratamento, y = Proporcao)) +
  geom_col(width = 0.5, fill = paleta_tcc[["geral"]]) +
  geom_text(aes(label = sprintf("%d / %d sessões", N_com_agitacao, N_sessoes)), vjust = -0.4, size = 3.5) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(taxa_agitacao$Proporcao) * 1.3)) +
  labs(
    x = NULL, y = "Proporção de sessões com Agitação",
    title = "Taxa observada de Agitação, por Nível de Visitação",
    subtitle = "Descritivo apenas — sem teste estatístico (ver Seção 13)"
  ) +
  tema_tcc +
  theme(legend.position = "none")

ggsave(
  file.path(pasta_resultados, "taxa_agitacao_aberto_fechado.png"),
  grafico_taxa_agitacao, width = 6.5, height = 4.5, dpi = 600, bg = "white"
)


# Seção 14 — gráficos descritivos: Tratamento (dados brutos)
dados_longo <- dados %>%
  select(Sessao, Tratamento, all_of(eixos_continuos)) %>%
  tidyr::pivot_longer(
    cols      = all_of(eixos_continuos),
    names_to  = "eixo_var",
    values_to = "proporcao"
  ) %>%
  mutate(eixo_var = nomes_eixos[eixo_var])

boxplot_tratamento <- dados_longo %>%
  ggplot(aes(x = Tratamento, y = proporcao, fill = Tratamento)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.3, color = paleta_tcc[["geral"]]) +
  facet_wrap(~ eixo_var, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(
    "Fechado"           = paleta_tcc[["geral"]],
    "Aberto (moderado)" = paleta_tcc[["Sono"]],
    "Aberto (intenso)"  = paleta_tcc[["Atividade"]]
  )) +
  labs(
    x = NULL, y = "Proporção do tempo da sessão",
    title = "Comportamento da Maya, por Nível de Visitação (dados brutos)",
    fill = "Nível de Visitação"
  ) +
  tema_tcc +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom",
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(
  file.path(pasta_resultados, "boxplot_tratamento.png"),
  boxplot_tratamento, width = 9, height = 7, dpi = 600, bg = "white"
)


# Seção 15 — gráficos de efeito marginal (dado bruto + suavização loess)
efeito_repouso_perturbacao <- dados %>%
  filter(Zoologico_aberto == "Aberto") %>%
  ggplot(aes(x = Perturbacao_s, y = Repouso_prop)) +
  geom_point(alpha = 0.6, color = paleta_tcc[["Repouso"]]) +
  geom_smooth(method = "loess", se = TRUE, color = paleta_tcc[["geral"]], fill = paleta_tcc[["geral"]]) +
  labs(
    x = "Perturbação na sessão (segundos)",
    y = "Repouso (proporção do tempo da sessão)",
    title = "Repouso x Perturbação (dias abertos)"
  ) +
  tema_tcc

ggsave(
  file.path(pasta_resultados, "efeito_marginal_repouso_perturbacao.png"),
  efeito_repouso_perturbacao, width = 7, height = 5, dpi = 600, bg = "white"
)

efeito_sono_visitantes <- dados %>%
  filter(Zoologico_aberto == "Aberto") %>%
  ggplot(aes(x = N_visitantes, y = Sono_prop)) +
  geom_point(alpha = 0.6, color = paleta_tcc[["Sono"]]) +
  geom_smooth(method = "loess", se = TRUE, color = paleta_tcc[["geral"]], fill = paleta_tcc[["geral"]]) +
  labs(
    x = "Nº de visitantes na sessão",
    y = "Sono (proporção do tempo da sessão)",
    title = "Sono x Nº de visitantes (dias abertos)"
  ) +
  tema_tcc

ggsave(
  file.path(pasta_resultados, "efeito_marginal_sono_visitantes.png"),
  efeito_sono_visitantes, width = 7, height = 5, dpi = 600, bg = "white"
)


# Seção 16 — % de zeros exatos por eixo (apoio visual, descritivo)
zeros_por_eixo <- sapply(eixos_continuos, function(v) mean(dados[[v]] == 0, na.rm = TRUE))

grafico_zeros <- data.frame(
  `Eixo comportamental` = nomes_eixos[eixos_continuos],
  Prop_zeros  = zeros_por_eixo,
  check.names = FALSE
) %>%
  ggplot(aes(x = reorder(`Eixo comportamental`, Prop_zeros), y = Prop_zeros)) +
  geom_col(fill = paleta_tcc[["geral"]], width = 0.5) +
  geom_text(aes(label = scales::percent(Prop_zeros, accuracy = 1)), hjust = -0.15, size = 3.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL, y = "% de sessões com valor exato 0",
    title = "Acúmulo de zeros por eixo comportamental"
  ) +
  tema_tcc

ggsave(
  file.path(pasta_resultados, "proporcao_zeros_por_desfecho.png"),
  grafico_zeros, width = 7, height = 4.5, dpi = 600, bg = "white"
)


# Seção 17 — tabela descritiva (N, média, DP, % de zeros), por Tratamento
tabela_descritiva <- dados_longo %>%
  group_by(eixo_var, Tratamento) %>%
  summarise(
    N          = n(),
    Média      = round(mean(proporcao), 3),
    `Desvio padrão` = round(sd(proporcao), 3),
    `% de zeros` = round(mean(proporcao == 0) * 100, 1),
    .groups = "drop"
  ) %>%
  rename(`Eixo comportamental` = eixo_var, `Nível de Visitação` = Tratamento)

write_xlsx(tabela_descritiva, file.path(pasta_resultados, "tabela_descritiva_tratamento.xlsx"))
salvar_tabela_png(tabela_descritiva, file.path(pasta_resultados, "tabela_descritiva_tratamento.png"), largura = 8)

cat("\nTabelas (.xlsx), tabelas em imagem (.png) e gráficos (forest plots, médias marginais, boxplot, efeito marginal, descritivos e métricas de desempenho) salvos em:\n", pasta_resultados, "\n")
