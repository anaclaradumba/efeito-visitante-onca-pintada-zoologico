# Orçamento Comportamental de *Panthera onca* sob Efeito da Visitação Pública

Este repositório contém a base de dados e a rotina de análise estatística em R utilizadas no trabalho de conclusão de curso sobre a fêmea de onça-pintada (*Panthera onca*) mantida no Jardim Zoológico de Belo Horizonte.

## 📂 Estrutura do Repositório

- `dados/Dados.xlsx`: planilha com as variáveis registradas por amostragem focal contínua no *software* BORIS (sessão, data, Tratamento/nível de visitação, temperatura, clima, nº de visitantes e os tempos por categoria comportamental do etograma).
- `scripts/modelo_unificado.R`: script único com todo o pipeline — preparação dos dados, padronização dos preditores de visitação, ajuste da Regressão Beta Ordenada por eixo comportamental (`glmmTMB`), diagnóstico residual (`DHARMa`), correção de Holm, cálculo de R² (`performance::r2_nakagawa`), tabelas e figuras (forest plots de Odds Ratio, médias marginais via `emmeans`, boxplots e estatística descritiva).

## 🛠️ Pacotes Utilizados no R

- `glmmTMB`
- `DHARMa`
- `broom.mixed`
- `emmeans`
- `performance`
- `ggplot2` / `gridExtra`
- `dplyr` / `tidyr`
- `writexl`

## 📧 Contato e Citação

- **Autora:** Ana Clara Dumbá Silva (UFMG)
- **Orientador:** Prof. Dr. Adriano Pereira Paglia
- **Coorientadora:** Me. Fernanda da Paixão Alves
