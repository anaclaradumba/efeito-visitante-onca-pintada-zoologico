# Orçamento Comportamental de *Panthera onca* sob Efeito da Visitação Pública

Dados e script de análise estatística em R do Trabalho de Conclusão de Curso sobre o orçamento comportamental da fêmea de onça-pintada (*Panthera onca*) mantida no Jardim Zoológico de Belo Horizonte, sob efeito da visitação pública.

## Estrutura

- `dados/dados-observaçao-onca-pintada.xlsx`: planilha com as variáveis registradas por amostragem focal contínua no BORIS (sessão, data, tratamento/nível de visitação, nº de visitantes e os tempos por categoria comportamental do etograma).
- `scripts/01_preparacao_dados.R`: lê a planilha, trata ausentes, deriva os eixos comportamentais (Repouso, Sono, Atividade, Agitação, Vigilância, Ocultação, Em pé) e os preditores de visitação, e salva os objetos usados pelo script seguinte.
- `scripts/modelo_unificado.R`: pipeline estatístico completo — regressão Beta Ordenada por eixo comportamental (`glmmTMB`), diagnóstico residual (`DHARMa`), correção de Holm, R² (`performance::r2_nakagawa`), tabelas e figuras (forest plots de Odds Ratio, médias marginais via `emmeans`, boxplots e estatística descritiva).

## Pacotes utilizados no R

`readxl`, `glmmTMB`, `DHARMa`, `broom.mixed`, `emmeans`, `performance`, `ggplot2`, `gridExtra`, `dplyr`, `tidyr`, `stringr`, `writexl`, `scales`.


## 📧 Contato e Citação

- **Autora:** Ana Clara Dumbá Silva ( anaclaradumba7@gmail.com )
- **Orientador:** Prof. Dr. Adriano Pereira Paglia ( apaglia@ufmg.br )
- **Coorientadora:** Prof. Dr. Paulo Enrique Cardoso Peixoto e (pauloenrique@gmail.com ) Me. Fernanda da Paixão Alves (fernandadapaixaoalves@gmail.com )
