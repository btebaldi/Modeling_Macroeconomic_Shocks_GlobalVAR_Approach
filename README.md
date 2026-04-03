
# Modeling how macroeconomic Shocks Affect Regional Employment: Analyzing the Brazilian Formal Labor Market Using the Global VAR Approach
GitHub repository: Modeling_Macroeconomic_Shocks_GlobalVAR_Approach

## 📘 Repositório de Reprodução - Artigo "Modeling how macroeconomic Shocks Affect Regional Employment: Analyzing the Brazilian Formal Labor Market Using the Global VAR Approach"

Este repositório contém os materiais necessários para reprodução e extensão dos resultados apresentados no artigo publicado na revista *Estudos Econômicos* ([https://revistas.usp.br/ee](https://revistas.usp.br/ee)):

**Título:** *Modelando como os choques macroeconômicos afetam o emprego regional: analisando o mercado de trabalho formal brasileiro usando a abordagem GVAR*  
**Title:** *Modeling how macroeconomic Shocks Affect Regional Employment: Analyzing the Brazilian Formal Labor Market Using the Global VAR Approach*  
**Data Publicação:** 25-03-2026  
**DOI:** https://doi.org/10.1590/1980-53575615bep  
**Link:** [revistas.usp.br/ee/pt_BR/article/view/227096](https://revistas.usp.br/ee/pt_BR/article/view/227096)

**Autores:**
- Bruno Tebaldi Q. Barbosa  
Professor na Escola de Economia de São Paulo - FGV;  
orcid: [0000-0001-9638-5651](https://orcid.org/0000-0001-9638-5651)

- Emerson Fernandes Marçal  
Professor na Escola de Economia de São Paulo - FGV;  
orcid: [0000-0002-0841-5644](https://orcid.org/0000-0002-0841-5644)

- Pedro L. Valls Pereira  
Professor na Escola de Administração de Empresas de São Paulo - FGV;  
orcid: [0000-0002-3709-9564](https://orcid.org/0000-0002-3709-9564)

### Citação

Para citar este repositório, arquivos ou materiais associados, utilize o formato apresentado abaixo, adaptando conforme necessário para o contexto do seu trabalho. 

Barbosa, B. T. (2025). *Repositorio do artigo Modeling how macroeconomic Shocks Affect Regional Employment: Analyzing the Brazilian Formal Labor Market Using the Global VAR Approach*. https://github.com/btebaldi/Modeling_Macroeconomic_Shocks_GlobalVAR_Approach/.

> @misc{barbosa2025RepositoryModeling_Macroeconomic_Shocks_GlobalVAR_Approach,  
author       = {Barbosa, Bruno Tebaldi},  
title        = {Repositório do Artigo {Modeling how macroeconomic Shocks Affect Regional Employment: Analyzing the Brazilian Formal Labor Market Using the Global VAR Approach}},  
year         = {2026},  
howpublished = {\url{https://github.com/btebaldi/Modeling_Macroeconomic_Shocks_GlobalVAR_Approach}}  
}


### 🧩 Conteúdo do Repositório

* 📂 `Base de dados/`: Bases de dados usadas no artigo, como séries temporais de admissões/demissões por mesorregião, produção industrial e a matriz de conexões regionais (IBGE).
* 📂 `Export/`: Resultados principais do artigo, incluindo tabelas, gráficos e saídas do modelo.
* 📂 `GVAR_Toolbox2.0/`: Códigos em MATLAB utilizados para estimar o modelo GVAR, incluindo scripts auxiliares e de pré-processamento de dados. Utiliza a toolbox **GVAR Toolbox 2.0**.
* 📂 `scripts/`: Scripts adicionais. Utilizado na construção de base de dados e análises de resultados.
* 📂 `Videos/`: Vídeos explicativos com visualizações das previsões regionais de emprego sob diferentes cenários macroeconômicos.

### 🔁 Reprodutibilidade

Todas as etapas do artigo — da preparação dos dados à estimação do modelo e à simulação de cenários — são totalmente reproduzíveis por meio dos scripts disponíveis. Consulte o arquivo `README.txt` em `GVAR_Toolbox2.0/` para instruções detalhadas de execução.


### 📌 Modelos Disponíveis

Este repositório disponibiliza diferentes versões do modelo estimado, com variações que permitem explorar alternativas metodológicas e aprofundar a análise. Abaixo, uma breve descrição de cada versão:

* **`Meso17`**: Modelo principal utilizado no artigo, estimado com dados mensais de 2004 a 2016 para 137 mesorregiões brasileiras. Utiliza a configuração padrão descrita no texto final publicado.
* **`Meso17_AllTests`**: Versão idêntica ao modelo `Meso17`, porém configurada para gerar e imprimir todos os testes estatísticos disponíveis no **GVAR Toolbox**.
* **`Meso17_FullSample`**: Variante do modelo `Meso17` que utiliza toda a amostra disponível sem divisão entre períodos de estimação e previsão.
* **`Meso19`**: Extensão do modelo que trata explicitamente a mesorregião metropolitana de São Paulo como uma unidade dominante, incorporando seu impacto direto sobre as demais regiões de forma diferenciada.

### 📌 Destaques do Artigo

* Estimação de um modelo GVAR para 137 mesorregiões brasileiras, com base em dados mensais de 2004 a 2016.
* Uso de uma matriz de pesos baseada em conexões econômicas entre municípios brasileiros (IBGE, 2008).
* Simulação de diferentes trajetórias de recuperação econômica após a recessão de 2014–2016, com foco em impactos regionais no emprego formal.
* Identificação de regiões mais e menos resilientes a choques macroeconômicos.
