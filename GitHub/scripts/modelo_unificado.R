## =============================================================================
## MODELO UNIFICADO
## =============================================================================
##
## Objetivo: em vez de um modelo separado por pergunta (P1: aberto x fechado;
## P2: perturbação/presença calma; P3: número de visitantes), ajustamos UM
## modelo por desfecho comportamental, com todas as covariáveis relevantes
## entrando juntas. Isso evita inflação de teste múltiplo (não estamos
## rodando 3 modelos independentes para a mesma variável) e dá estimativas
## mutuamente ajustadas (o efeito de "perturbação" já é calculado controlando
## pelo número de visitantes e vice-versa).
##
## VERSÃO RECONSTRUÍDA A PARTIR DA TRANSCRIÇÃO DA CONSULTORIA COM A FERNANDA
## (reunião "Reunião Paulo.m4a"). Cada decisão abaixo cita a fala exata da
## consultora (ou a literatura discutida com o Natã na mesma conversa em que
## este script foi escrito), sem nenhuma suposição extra:
##
##   1. Unidade de análise é a SESSÃO (10 minutos), não o dia. Cada sessão
##      aparece 1 única vez na base (não é medida repetidamente). Por isso,
##      SEM decomposição "entre-dia" / "dentro-dia" dos preditores de
##      visitação (perturbação, presença calma, nº de visitantes): cada um
##      entra como UM único valor por sessão. A Fernanda foi direta sobre
##      isso: "não precisa fazer essa distinção, porque no fundo a sua
##      pergunta é: tem diferença entre o nível de movimento... em relação ao
##      comportamento das pessoas no momento em que a observação está
##      acontecendo? que são as suas sessões."
##
##   2. SEM efeito aleatório de Sessão: "quando você bota, você tira a
##      sessão, tanto é que deve ter dado singular... porque ela não vai
##      explicar nada" — cada sessão só acontece 1 vez, então um efeito
##      aleatório por sessão não tem o que agrupar.
##
##   3. Efeito aleatório final é só por DIA: "(1 | dia)". A fórmula final que
##      a Fernanda escreveu, de próprio punho, no chat da reunião foi:
##      "proporção de tempo ~ tratamento + (1 | dia)". Esse "(1 | dia)" é
##      MANTIDO no modelo sempre, sem nenhum teste formal prévio de
##      necessidade — a consultora deu a fórmula pronta, sem condicionar a
##      nenhum teste. Isso também é consistente com a prática padrão para
##      efeitos aleatórios ligados à estrutura de AMOSTRAGEM do estudo
##      (sessões aninhadas dentro de dias) — diferente de um efeito aleatório
##      exploratório, esses costumam ser mantidos por justificativa de
##      desenho, não por teste de significância (Zuur et al., 2009, "Mixed
##      Effects Models and Extensions in Ecology with R"; Bolker et al.,
##      2009, Trends in Ecology & Evolution, "Generalized linear mixed
##      models: a practical guide for ecology and evolution").
##
##      OBSERVAÇÃO: num momento anterior da mesma reunião, antes de vocês
##      chegarem no desenho final com Tratamento como fator fixo de 3 níveis,
##      a Fernanda chegou a sugerir DUAS variáveis aleatórias (o dia
##      numérico, de 1 a 10, e o tipo de dia — segunda/quarta/domingo). Essa
##      sugestão foi SUPERADA mais adiante na mesma conversa, quando ela
##      fecha com a fórmula única "tratamento + (1 | dia)" — é essa versão
##      final que este script implementa.
##
##   4. "Zoologico_aberto" (2 níveis: Aberto/Fechado) vira "Tratamento" (3
##      níveis: Fechado / visitação moderada / visitação intensa), porque os
##      7 dias abertos do desenho não têm todos a mesma intensidade de
##      visitação pública — quarta-feira e domingo são bem diferentes entre
##      si. A própria Fernanda orienta: "o seu modelo tem que ser, na sua
##      editora, tem que ser fechado, mede o movimento e muito movimento...
##      fica com três barrinhas. O resto é tudo igual, não vai mudar nada."
##      (Rótulo "intenso" em vez de "muito" — decisão de nomenclatura do
##      Natã, não muda nenhum cálculo.)
##
##   5. Distribuição: BETA ORDENADA (`ordbeta()` no glmmTMB), não a Beta
##      padrão (`beta_family()`). A Beta padrão (Ferrari & Cribari-Neto,
##      2004) só é definida no intervalo ABERTO (0, 1) — não aceita valores
##      exatos de 0 ou 1, que é exatamente o que ocorre nestes dados (sessões
##      com 0% de tempo em Sono, Ocultação etc.). A beta ordenada (Kubinec,
##      2022, "Ordered Beta Regression: A Parsimonious, Well-Fitting Model
##      for Continuous Data with Bounded Support", Political Analysis) foi
##      construída justamente para esse caso: modela, numa única
##      especificação, o acúmulo nos extremos (0 e/ou 1) e a variação
##      contínua no meio, sem precisar transformar os zeros artificialmente.
##      É também a família que a própria Fernanda reconhece que vocês já
##      vinham usando ("você usou o modelo de LNM com distribuição beta
##      ordenada... foi bom só na atividade").
##
##      Usamos UM único modelo ordbeta por desfecho (sem testar variantes com
##      dispersão variável por Tratamento nem modelo hurdle/zero-inflado) —
##      decisão explícita do Natã para manter o R² computável de forma
##      direta nos 4 eixos (ver Seção 7 sobre R²).
##
##   6. Significância dos preditores: correção de HOLM (`p.adjust(method =
##      "holm")`), aplicada aos 5 p-valores de Wald de cada modelo (2
##      coeficientes de Tratamento + Perturbação + Presença calma + Nº de
##      visitantes), tratando cada modelo (cada desfecho) como sua própria
##      família de hipóteses. A Fernanda prefere o Teste de Razão de
##      Verossimilhança, mas validou explicitamente que Holm não está errado
##      e não muda as conclusões: "fazendo com P de Home, do jeito que eu
##      falei, vai dar a mesma coisa" e "isso pro seu TCC está muito
##      tranquilo... não está errado" (ela só recomendaria TRV para quem for
##      publicar artigo depois).
##
##   7. Métricas de desempenho: R² MARGINAL e R² CONDICIONAL (Nakagawa &
##      Schielzeth, 2013, 2017), via `performance::r2_nakagawa()`. SEM ICC —
##      "o ICC você não precisa, porque ele é um componente para quem está
##      interessado em personalidade." SEM R² de Nagelkerke — "o R de Nagel
##      Kerc é feito para modelos fixos... GLMM, TNB, ele não funciona."
##      Definições da própria Fernanda: "o R marginal é a porcentagem da
##      variação do comportamento explicada pelos fatores de visitação" e "o
##      condicional é a variação que é explicada só pelas variáveis
##      aleatórias."
##
##   8. Zero não é, por si só, um problema: "zero não é um problema, a não
##      ser que os seus resíduos fiquem bizarros... se o Dharma falou que os
##      resíduos estão belíssimos, não tem problema de ter 70% de zero." O
##      critério é o diagnóstico de resíduos simulados (pacote DHARMa:
##      uniformidade e dispersão), não a % de zeros em si — reportado na
##      Seção 6 deste script, eixo a eixo.
##
##   9. Odds Ratio (OR): mantido nos forest plots/tabelas (a Fernanda não o
##      veta), mas complementado com uma figura de médias marginais
##      estimadas por nível de Tratamento, na escala da proporção (não em
##      log-odds/OR) — ela prefere essa leitura: "eu prefiro a figura, porque
##      eu olho com a figura e entendo que está acontecendo assim." Essa
##      figura complementar usa o pacote `emmeans` (Lenth, 2016, "Least-
##      Squares Means: The R Package lsmeans", Journal of Statistical
##      Software), que é o pacote com suporte direto e documentado pelo
##      próprio glmmTMB para extrair médias marginais estimadas (EMMs) e
##      contrastes par a par já na escala de resposta (via `regrid(type =
##      "response")`), e não exige nenhum pacote adicional.
##
## Agitação (pacing) NÃO tem modelo estatístico: com só 7 ocorrências em 100
## sessões, todas em dias abertos, a Fernanda foi categórica: "não existe
## nenhuma análise que vai resolver o pacing... quando não tinha ninguém, era
## zero. Quando tinha gente, a chance dela fazer pace era 7,70, 10%... porque
## um teste que vai testar variação contra zero, ele vai dar essas loucuras
## aí. Porque zero não varia." Fica só a proporção observada, descritiva.
##
## Pacotes citados na consultoria/discussão: glmmTMB, DHARMa, emmeans,
## broom.mixed, writexl, ggplot2, gridExtra, performance.
## =============================================================================

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
library(writexl)    ## salvar tabelas em .xlsx (nunca em .csv — decisão do Natã)
library(ggplot2)    ## gráficos (forest plots, médias marginais, boxplots)
library(gridExtra)  ## transformar tabelas (data.frame) em imagem PNG
library(performance) ## R² de Nakagawa para os modelos glmmTMB


## =============================================================================
## SEÇÃO 0A — PALETA DE CORES E TEMA GRÁFICO DO TCC (DECISÃO DO NATÃ)
## =============================================================================
## Paleta fixa usada em todas as tabelas/gráficos do TCC. "geral" é a cor
## PRINCIPAL ÚNICA usada em praticamente todos os gráficos e tabelas (forest
## plots, médias marginais, tabelas de coeficientes, tabelas descritivas/
## diagnóstico/métricas, gráfico de barras de zeros, gráfico descritivo de
## Agitação). As demais 5 cores ficam reservadas EXCLUSIVAMENTE para os
## gráficos que comparam grupos/séries e por isso exigem cores visualmente
## distintas: o boxplot por Tratamento e os gráficos de dispersão de efeito
## marginal (ponto observado x linha de tendência).
paleta_tcc <- c(
  geral     = "#8c4327",
  Repouso   = "#8b5a2b",
  Sono      = "#d95a2b",
  Atividade = "#4a5d36",
  Ocultação = "#e6a135",
  Agitação  = "#c2b280"
)

## Tema comum a todos os gráficos do TCC: sem grades de fundo (painel e eixos
## limpos), mantendo só o essencial para leitura (eixo e marcações).
tema_tcc <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x      = element_line(color = "grey40"),
    plot.title       = element_text(face = "bold"),
    plot.margin      = margin(0, 0, 0, 0)
  )

## =============================================================================
## SEÇÃO 0 — CAMINHOS (ABSOLUTOS, FIXOS)
## =============================================================================
## Nenhum caminho relativo neste script: tudo é lido e salvo por caminho
## absoluto, então não importa qual é o diretório de trabalho do RStudio
## nem onde este arquivo está salvo.
pasta_resultados <- "C:/Users/User/OneDrive/Desktop/TCC/resultados/atual"

cat("Pasta de resultados:", pasta_resultados, "\n")
cat("Arquivo dados_descritiva.rds encontrado?", file.exists(file.path(pasta_resultados, "dados_descritiva.rds")), "\n")


## =============================================================================
## SEÇÃO 1 — LEITURA DOS DADOS
## =============================================================================
## Reaproveita o que já foi preparado em analise_descritiva.Rmd:
##   - Zoologico_aberto recodificado para "Aberto"/"Fechado"
##   - Perturbacao_s = GRI_s + BAT_s + ASS_s + PAL_s
##   - Presenca_calma_s = CNV_s + OBS_s
##   - Repouso_prop, Sono_prop, Atividade_prop, Ocultacao_prop (eixos contínuos)
##   - Agitacao_s (duração) e, a partir dela, PAC_ocorreu (ocorreu/não)
##
## Esse .rds precisa estar atualizado (ou seja: analise_descritiva.Rmd
## precisa ter sido rodado de novo no RStudio) antes de rodar este script.

dados <- readRDS(file.path(pasta_resultados, "dados_descritiva.rds"))


## =============================================================================
## SEÇÃO 2 — TRATAMENTO (3 NÍVEIS: FECHADO / MODERADO / INTENSO)
## =============================================================================
## Em vez do binário Zoologico_aberto (Aberto/Fechado), o preditor principal
## passa a distinguir também a intensidade da visitação dentro dos dias
## abertos, porque os 7 dias abertos do desenho não são todos iguais em
## movimento — quarta-feira tem visitação moderada e domingo tem visitação
## bem maior (ver citação completa da Fernanda no topo do script, item 4).
## Como cada dia só tem UM nível de tratamento (não varia sessão a sessão
## dentro do mesmo dia), isso é construído direto a partir de
## Zoologico_aberto + Dia_semana, sem precisar de nenhuma conta nova:
##
##   Segunda-feira  -> Fechado           (zoológico não abre)
##   Quarta-feira    -> Aberto (moderado)
##   Domingo         -> Aberto (intenso)
##
## Usamos grepl() em vez de comparar a string exata porque não controlamos
## a formatação de Dia_semana na planilha original (pode vir com ou sem
## acento, maiúscula, "-feira" etc.) — o stopifnot() logo abaixo garante que
## nenhuma linha fique sem classificação (se a suposição acima sobre os
## valores de Dia_semana estiver errada para a sua planilha, o script para
## aqui com um erro, em vez de classificar algo errado silenciosamente).
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


## =============================================================================
## SEÇÃO 3 — PADRONIZAÇÃO (Z-SCORE) DOS PREDITORES CONTÍNUOS
## =============================================================================
## Sem decomposição entre-dia/dentro-dia (ver item 1 do topo do script).
## Perturbação, presença calma e nº de visitantes entram como UM único valor
## por sessão (o dado bruto da sessão), só padronizado (z-score) para que os
## coeficientes fiquem numa escala comparável entre si ("1 desvio-padrão",
## em vez de "1 segundo" ou "1 visitante").
##
## A média e o desvio-padrão usados na padronização são calculados só com
## os valores de dias abertos (onde essas variáveis de fato variam — em dia
## fechado não há visitante, então o valor bruto já é 0); nos dias fechados,
## o valor padronizado continua sendo 0, para não inventar um "desvio-padrão
## negativo" sem sentido nesses dias.

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


## =============================================================================
## SEÇÃO 4 — DESFECHOS CONTÍNUOS E FÓRMULA UNIFICADA
## =============================================================================
## 4 desfechos contínuos (proporção de tempo da sessão). HIG não entra como
## controle — foi removido por instrução explícita em conversa anterior.

eixos_continuos <- c("Repouso_prop", "Sono_prop", "Atividade_prop", "Ocultacao_prop")

nomes_eixos <- c(
  Repouso_prop   = "Repouso",
  Sono_prop      = "Sono",
  Atividade_prop = "Atividade",
  Ocultacao_prop = "Ocultação"
)

## Mesma fórmula para os 4 desfechos, só trocando a variável resposta
## (preenchida via sprintf mais abaixo). "(1 | Data)" é o efeito aleatório
## único e final (ver item 3 do topo do script) — sem efeito de Sessão (cada
## sessão só acontece 1 vez na base).
formula_unificada <- paste(
  "%s ~ Tratamento +",
  "Perturbacao_z + Presenca_calma_z + N_visitantes_z +",
  "(1 | Data)"
)


## =============================================================================
## SEÇÃO 5 — AJUSTE DOS 4 DESFECHOS CONTÍNUOS (BETA ORDENADA, ÚNICO MODELO)
## =============================================================================
## Um único modelo `ordbeta()` por desfecho, sem testar variantes (sem
## dispersão variando por Tratamento, sem hurdle/zero-inflado) — decisão
## explícita do Natã para manter o R² computável de forma direta e uniforme
## nos 4 eixos (ver item 5 e Seção 7). O diagnóstico de resíduos (Seção 6)
## é reportado para cada modelo tal como ele sai daqui, sem seleção
## automática entre candidatos.

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


## =============================================================================
## SEÇÃO 6 — DIAGNÓSTICO DE RESÍDUOS (DHARMa) — ZERO SÓ É PROBLEMA SE O
## RESÍDUO FALHAR
## =============================================================================
## "Zero não é um problema, a não ser que os seus resíduos fiquem bizarros...
## se o Dharma falou que os resíduos estão belíssimos, não tem problema de
## ter 70% de zero." — por isso, o critério de adequação do modelo é o
## resultado do DHARMa (uniformidade dos resíduos e dispersão), não a
## proporção de zeros nos dados. Essa tabela registra, para os 4 desfechos, o
## resultado exato dos 3 testes do DHARMa (uniformidade/KS, dispersão,
## outliers) do modelo único ajustado na Seção 5 — para citar o número exato
## na redação do TCC, sem precisar voltar a rodar nada, e para decidir com
## transparência se algum eixo precisa de uma ressalva na seção de
## limitações (não para escolher ou descartar um modelo).

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


## =============================================================================
## SEÇÃO 7 — EXTRAÇÃO DE RESULTADOS E CORREÇÃO DE HOLM
## =============================================================================
## Correção de Holm (`stats::p.adjust(method = "holm")`), aplicada aos 5
## p-valores de Wald de cada modelo (2 coeficientes de Tratamento + 3
## preditores contínuos), tratando cada modelo (cada desfecho) como sua
## própria família de hipóteses — não há correção cruzada entre os 4
## desfechos. A Fernanda validou que Holm não está errado e não muda as
## conclusões em relação ao Teste de Razão de Verossimilhança (ver item 6 do
## topo do script) — por isso os 2 coeficientes de Tratamento aqui têm
## p-valores de Holm INDEPENDENTES um do outro (cada um corrigido dentro da
## família de 5 hipóteses do seu modelo), diferente do TRV (que trataria
## Tratamento como um único preditor de 2 graus de liberdade).

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

## Junta o p-valor do teste de dispersão de volta na tabela de coeficientes,
## para que cada linha já carregue, junto do próprio efeito, se o modelo de
## onde ela veio passou ou não na verificação de resíduos.
resultados_tabela <- resultados_tabela %>%
  left_join(
    diagnostico_residuos_tabela %>% select(eixo, dispersao_p, uniformidade_p),
    by = "eixo"
  )

cat("\n--- Resultados (desfechos contínuos), com p de Holm ---\n")
print(resultados_tabela)


## =============================================================================
## SEÇÃO 8 — SALVAR RESULTADOS (SEMPRE EM .xlsx, NUNCA EM .csv)
## =============================================================================
## Os .rds guardam os objetos de modelo completos (para reabrir no R depois,
## sem precisar re-rodar tudo). As tabelas de resultado — as que vão para o
## texto do TCC — são salvas em Excel (.xlsx), por decisão explícita: nada
## de .csv neste script.

saveRDS(resultados_continuos, file.path(pasta_resultados, "modelos_unificados_continuos.rds"))

write_xlsx(resultados_tabela,           file.path(pasta_resultados, "resultados_unificados_continuos.xlsx"))
write_xlsx(diagnostico_residuos_tabela, file.path(pasta_resultados, "resultados_diagnostico_residuos.xlsx"))


## =============================================================================
## SEÇÃO 9 — TABELAS FORMATADAS PARA LEITURA (E EM IMAGEM PNG)
## =============================================================================
## As tabelas da Seção 8 têm a precisão numérica completa (boas para reabrir
## no R e continuar a análise). Para citar no texto do TCC, isso não é
## legível: precisão de 10+ casas decimais, nomes de variável crus
## (ex.: "Perturbacao_z") e p-valores em notação científica.
##
## Esta seção cria uma SEGUNDA versão de cada tabela, só para exibição:
## números arredondados, nomes de termo traduzidos para um rótulo legível,
## e p-valores formatados (ex.: "< 0,001" em vez de "2.2e-16"). Essa versão
## também é salva em .xlsx e, além disso, transformada em imagem .png — o
## formato pedido para colar direto no texto do TCC.

## Dicionário: nome técnico do termo -> rótulo legível para a tabela.
## glmmTMB nomeia cada nível não-referência de um fator como
## "<nome_da_variável><nome_do_nível>" colado, por isso os dois primeiros
## rótulos abaixo não têm espaço entre "Tratamento" e o nível.
rotulos_termo <- c(
  "TratamentoAberto (moderado)" = "Aberto, visitação moderada (vs. Fechado)",
  "TratamentoAberto (intenso)"  = "Aberto, visitação intensa (vs. Fechado)",
  "Perturbacao_z"               = "Perturbação",
  "Presenca_calma_z"            = "Presença calma",
  "N_visitantes_z"              = "Nº de visitantes"
)

## Formata um p-valor: menor que 0.001 vira "< 0.001"; senão, 3 casas decimais.
formatar_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
}

## --- Tabela 1: resultados dos 4 desfechos contínuos -------------------------
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

## --- Tabela 2: diagnóstico de resíduos (DHARMa), por desfecho ---------------
tabela_diagnostico_exibicao <- diagnostico_residuos_tabela %>%
  transmute(
    `Eixo comportamental` = eixo,
    `Família`            = familia,
    `% de zeros`         = round(prop_zeros * 100, 1),
    `Uniformidade (p)`   = round(uniformidade_p, 3),
    `Dispersão (p)`      = round(dispersao_p, 3),
    `Outliers (p)`       = round(outliers_p, 3)
  )

## Salva as 2 versões de exibição também em .xlsx (além das brutas da Seção 8).
write_xlsx(tabela_continuos_exibicao,   file.path(pasta_resultados, "tabela_resultados_continuos_exibicao.xlsx"))
write_xlsx(tabela_diagnostico_exibicao, file.path(pasta_resultados, "tabela_diagnostico_exibicao.xlsx"))

## Função auxiliar: desenha um data.frame como tabela colorida (paleta do
## TCC) e salva em .png. Por decisão do Natã, TODAS as tabelas (descritivas,
## metodológicas, de métricas e também as de coeficientes por desfecho)
## usam a mesma cor principal única `paleta_tcc[["geral"]]` no cabeçalho —
## não há cor por desfecho nas tabelas. `largura` controla a largura da
## imagem, em polegadas; a altura é calculada automaticamente a partir do
## nº de linhas.
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
  ## Usa a largura/altura REAIS do grob (soma das colunas/linhas), em vez de
  ## um canvas de tamanho fixo maior que a tabela — era isso que sobrava
  ## como borda branca ao redor do conteúdo. `largura` vira um teto
  ## opcional: se o grob for mais estreito que `largura`, usa o tamanho do
  ## próprio grob.
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


## =============================================================================
## SEÇÃO 10 — FOREST PLOTS EM RAZÃO DE CHANCES (OR), PALETA PRÓPRIA DO TCC
## =============================================================================
## Forest plots em escala de Razão de Chances (OR = exp(estimativa)).
## Motivo: a família beta ordenada usa link logit para a média, então
## exp(estimativa) tem uma leitura de razão de chances (o quanto a razão
## média/(1-média) multiplica por unidade do preditor). O IC95% do OR é
## obtido por exponenciação direta do IC95% em log-odds — não é recalculado
## por outro método — então a informação estatística é idêntica; muda só a
## escala de leitura. A Fernanda não veta o OR (ver item 9 do topo do
## script); essa figura é complementada, não substituída, pela figura de
## médias marginais da Seção 11.
##
## Paleta de cores fixa do TCC (Seção 0A); sem grades de fundo (tema
## `tema_tcc`); valores de OR, IC95% e p dentro de uma caixinha posicionada
## acima da bolinha de cada preditor.

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

## A significância (destaque visual) usa o p de Wald corrigido por Holm
## (p_holm) — ver Seção 7.
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

## Tabela de coeficientes em OR, uma imagem por desfecho contínuo (usada no
## documento de Métodos/Resultados, uma tabela logo antes de cada forest
## plot).
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


## =============================================================================
## SEÇÃO 11 — FIGURA COMPLEMENTAR: MÉDIAS MARGINAIS POR TRATAMENTO (emmeans)
## =============================================================================
## Complemento ao forest plot em OR (Seção 10), na escala da proporção (não
## em log-odds/OR) — a Fernanda disse preferir essa leitura: "eu prefiro a
## figura, porque eu olho com a figura e entendo que está acontecendo
## assim." Usamos `emmeans` (Lenth, 2016) para extrair, por nível de
## Tratamento, a proporção de tempo média estimada pelo modelo (já
## retrotransformada para a escala de resposta via `regrid(type =
## "response")`) com IC95%, e os contrastes par a par entre níveis (também
## na escala de resposta, em pontos percentuais), com a mesma correção de
## Holm usada na Seção 7. Essa figura só se aplica ao preditor categórico
## (Tratamento) — para os 3 preditores contínuos, o efeito já é mostrado na
## escala de OR pelo forest plot da Seção 10.

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

  ## Médias marginais estimadas por nível de Tratamento, na escala de
  ## resposta (proporção de tempo).
  emm <- emmeans(r$modelo, ~ Tratamento, type = "response")
  emm_df <- as.data.frame(emm)

  ## glmmTMB com família ordbeta devolve a coluna de média estimada com o
  ## nome "response" quando type = "response" é usado.
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

  ## Contrastes par a par entre níveis de Tratamento, na escala de resposta
  ## (pontos percentuais), com correção de Holm — mesma lógica de correção
  ## usada na Seção 7 para os coeficientes do modelo.
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


## =============================================================================
## SEÇÃO 12 — R² MARGINAL E R² CONDICIONAL (SEM ICC, SEM NAGELKERKE)
## =============================================================================
## R² marginal e R² condicional (Nakagawa & Schielzeth, 2013, 2017), via
## `performance::r2_nakagawa()`. SEM ICC — "o ICC você não precisa, porque
## ele é um componente para quem está interessado em personalidade." SEM R²
## de Nagelkerke — "o R de Nagel Kerc é feito para modelos fixos... GLMM,
## TNB, ele não funciona." Definições da Fernanda: "o R marginal é a
## porcentagem da variação do comportamento explicada pelos fatores de
## visitação" e "o condicional é a variação que é explicada só pelas
## variáveis aleatórias."
##
## Como os 4 modelos agora usam a beta ordenada SIMPLES (sem dispersão
## variando por Tratamento — Seção 5), a decomposição de variância do
## `performance::r2_nakagawa()` (que assume variância residual constante)
## se aplica sem nenhuma restrição.
##
## OBSERVAÇÃO (ajuste singular de (1 | Data)): como Tratamento já é uma
## variável de NÍVEL DE DIA (cada dia tem um único valor de Tratamento),
## depois de incluí-lo como efeito fixo sobra muito pouca variação "entre
## dias" para o `(1 | Data)` explicar — nos 4 eixos, o glmmTMB estima essa
## variância como ~0 (ajuste singular). Por padrão, `performance::
## r2_nakagawa()` se recusa a calcular o R² condicional nesse caso e
## devolve NA (ver `?performance::r2_nakagawa`, argumento `tolerance`).
## Isso não significa que o R² condicional "verdadeiro" seja outro valor
## escondido: com variância do efeito aleatório ~0, R² condicional
## converge para o mesmo valor do R² marginal — a fórmula de Nakagawa &
## Schielzeth soma a variância do efeito aleatório só no numerador do
## condicional, e essa variância é ~0. Usamos `tolerance = 1e-10` (bem
## abaixo do padrão do pacote) para que ele aceite esse valor pequeno em
## vez de recusar o cálculo — isso não infla nem altera a variância
## estimada pelo glmmTMB em nenhum momento, só evita que o pacote devolva
## NA por considerá-la "numericamente indistinguível de zero".

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


## =============================================================================
## SEÇÃO 13 — AGITAÇÃO (PACING): SÓ DESCRITIVO, SEM MODELO ESTATÍSTICO
## =============================================================================
## Com só 7 ocorrências em 100 sessões, todas em dias abertos, a Fernanda foi
## categórica: "não existe nenhuma análise que vai resolver o pacing...
## quando não tinha ninguém, era zero. Quando tinha gente, a chance dela
## fazer pace era 7,70, 10%... porque um teste que vai testar variação
## contra zero, ele vai dar essas loucuras aí. Porque zero não varia." Fica
## só o dado observado: quantas sessões tiveram Agitação, do total, por
## nível de Tratamento.

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

## Cópia só para exibição (tabela exportada): renomeia a coluna Tratamento
## para o rótulo pedido, sem alterar `taxa_agitacao` (usada no gráfico logo
## abaixo, que referencia a coluna original `Tratamento`).
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


## =============================================================================
## SEÇÃO 14 — GRÁFICOS DESCRITIVOS: TRATAMENTO (DADOS BRUTOS)
## =============================================================================
## Antes dos resultados do modelo, vale mostrar o dado bruto: como os 4
## desfechos contínuos se distribuem entre os 3 níveis de Tratamento, sem
## nenhum ajuste estatístico por trás. Boxplot (mediana e quartis) sobreposto
## a pontos individuais (cada ponto = 1 sessão) mostra tanto o padrão
## central quanto a dispersão real e o tamanho da amostra em cada grupo —
## mais informativo do que só a média. Com 3 grupos, usamos 3 cores
## distintas da paleta (gráfico de comparação entre grupos).

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


## =============================================================================
## SEÇÃO 15 — GRÁFICOS DE EFEITO MARGINAL (DADO BRUTO + SUAVIZAÇÃO)
## =============================================================================
## O forest plot (Seção 10) mostra o tamanho do efeito na escala do modelo
## (padronizada, abstrata). Complementamos com um gráfico de efeito
## marginal: o dado real (não padronizado, direto da sessão) no eixo x, a
## proporção do desfecho no eixo y, com linha de tendência suavizada
## (loess). Mostra a FORMA da relação, na unidade original (segundos de
## perturbação, número de visitantes) — mais fácil de explicar no texto do
## que "1 desvio-padrão", e é esse tipo de leitura visual direta que a
## Fernanda disse preferir.
##
## Repouso x Perturbação e Sono x Nº de visitantes são os pares que fizeram
## mais sentido na análise anterior; depois de rodar o modelo, vale conferir
## se os mesmos dois efeitos continuam sendo os significativos no forest
## plot da Seção 10 e, se não, trocar aqui para o par que saiu significativo.

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


## =============================================================================
## SEÇÃO 16 — APOIO VISUAL À % DE ZEROS EXATOS (DESCRITIVO, NÃO É RESSALVA
## AUTOMÁTICA)
## =============================================================================
## Mostra, para os 4 desfechos contínuos lado a lado, a % de sessões com
## valor exato zero — dado bruto e descritivo. Por si só NÃO indica problema
## algum (ver Seção 6: o critério de ressalva é o diagnóstico de resíduos,
## não a % de zeros). Este gráfico serve só de apoio visual para contextuar
## os números, ao lado da tabela de diagnóstico de resíduos.

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


## =============================================================================
## SEÇÃO 17 — TABELA DESCRITIVA (N, MÉDIA, DP, % DE ZEROS), POR TRATAMENTO
## =============================================================================
## Tabela clássica de estatística descritiva, por desfecho e por nível de
## Tratamento — o tipo de tabela que normalmente abre a seção de Resultados
## de um TCC, antes de entrar nos modelos. Complementa o boxplot (Seção 14)
## com os números exatos por trás dele.

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
