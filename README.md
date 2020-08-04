# Pipeline de Farmacogenética da Genomika 🧬🖥️💊
## Pipeline para definição de haplótipos utilizado no exame de farmacogenética da Genomika

Ess repositório contém o código  (além da sua respectiva documentação) do pipeline de farmacogenética utilizado na Genomika. Esse pipeline tem como objetivo a definição de haplótipos para um conjunto específico de genes, conseguindo assim criar recomendações de metabolização para um dado grupo de medicamentos.

----

### Arquitetura Proposta
Esse pipeline se divide essencialmente em quatro partes:

1. A descoberta de pequenas variantes (SNPs e Indels) de acordo com as [melhores práticas propostas pelo GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-);
2. Análise de cobertura usando o `DepthOfCoverage` do GATK;
3. Definição de haplótipos e criação de recomendação de metabolização usando o [**Stargazer**](https://stargazer.gs.washington.edu/stargazerweb/);
4. Organização de arquivos resultantes por amostra processada;

Como o Stargazer é responsável pelo produto final do processamento, a parametrização e versão dos softwares utilizados aqui foi feita de acordo com o que é proposto pela ferramenta [`pipeline`](https://stargazer.gs.washington.edu/stargazerweb/res/documentation.html#pipeline) do Stargazer. Dessas propostas é importante notar que foi recomendado a utilização do **GATK 3** além do **GRCh37** como genoma de referência.

Todo o processamento é orquestrado pelo [**Nextflow**](https://www.nextflow.io/), e cada etapa é executada dentro de um container [Docker](https://www.docker.com/) específico.

----

### Requisitos para execução

Para executar esse pipeline é necessário um ambiente de execução com os seguintes requisitos:

* Nextflow instalado (versão recomendada `20.04.1.5335`) ([guia de instalação](https://www.nextflow.io/index.html#GetStarted));
* Docker instalado ([guia de instalação](https://www.docker.com/get-started));
* Imagens Docker disponíveis na máquina;
* Diretório contendo o genoma de referência GRCh37, além de seus respectivos arquivos de índice e dicionário;
* Diretório com o pacote de recursos do genoma `b37` da Broad Institute.

#### Preparando as imagens Docker


Para construir a imagem necessária para a etapa de alinhamento (imagem contém ambos [`bwa`](https://github.com/lh3/bwa) (v0.7.17) e [`samtools`](https://github.com/samtools/samtools) (1.10)), vá para o diretório `dockerfiles/alignment` e execute o comando:

```bash
$ docker build . -t alignment:v0.1.0
```

Para construir a imagem necessária para executar o [`Picard`](https://broadinstitute.github.io/picard/), execute o comando:

```bash
$ docker pull broadinstitute/picard:2.22.8
```

Para construir a imagem necessária para executar o [GATK 3](https://github.com/broadgsa/gatk/releases/tag/3.8-1), execute o comando:

```bash
$ docker pull broadinstitute/gatk3:3.8-0
```

Para construir a imagem necessária para executar o Stargazer, vá para o diretório `dockerfiles/stargazer` execute o comando:

```bash
$ docker build . -t stargazer:v1.0.8
```

Para construir a imagem necessária para importar o framework [`pandas`](https://pandas.pydata.org/) dentro do python, vá para o diretório `dockerfiles/pandas` execute o comando:

```bash
$ docker build . -t pandas:1.0.5
```

Para construir a imagem necessária para usar o utilitário de linha de comando [`pdfunite`](https://manpages.debian.org/buster/poppler-utils/pdfunite.1.en.html), vá para o diretório `dockerfiles/poppler` execute o comando:

```bash
$ docker build . -t poppler:0.82.0-r1
```

Após construir e baixar todas as imagens acima, execute o seguinte comando para limpar seu ambiente:

```bash
$ docker rmi $(docker images --filter dangling=true -q)
```

#### Preparando o genoma de referência

1. Faça [download do arquivo `.fasta`](ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/) do genoma de referência GRCh37;
2. Certifique-se que o nome do arquivo baixado é `Homo_sapiens.GRCh37.dna.primary_assembly.fa`;
3. Coloque o arquivo dentro de um diretório específico, e dentro desse diretório crie os índices do genoma utilizando o comando
```bash
$ samtools faidx <fasta>
```
4. Ainda dentro desse diretório, crie o arquivo de dicionário do genoma de referência utilizando o comando
```bash
$ picard CreateSequenceDictionary R=<fasta> O=<nome base>.dict
```

#### Pacotes de Recurso da Broad Institute

Os genomas de referência além dos pacotes de recursos utilizados na Broad são disponibilizados através do seu servidor de ftp. As intruções de como baixar os arquivos dos pacotes podem ser encontradas [nesse artigo](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

Faça download do pacote de recursos do genoma `b37` e coloque os arquivos dentro de um diretório específico. Certifique-se que os seguintes arquivos estão dentro desse diretório:

1. `dbsnp_138.b37.vcf`
2. `Mills_and_1000G_gold_standard.indels.b37.vcf`
3. `1000G_omni2.5.b37.vcf`
4. `1000G_phase1.indels.b37.vcf`

Note que os arquivos precisam estar descompactados.

----

### Executando o pipeline

Uma vez preparado o ambiente, podemos rodar o pipeline. Como dito anteriormente, toda a execução é orquestrada pelo Nextflow, e as instruções de processamento estão presentes no arquivo `main.nf`. Esse arquivo recebe os seguintes argumentos de linha de comando:

* `--samples`: Caminho de diretório contendo arquivos `.fastq.gz` da batch de amostras que será processada;
* `--reference`: Caminho do diretório onde se encontra o `fasta` do genoma de referência, além dos arquivos de dicionário e índice;
* `--resources`: Caminho do diretório onde se encontram os arquivos do pacote de recursos criado pela Broad do genoma `b37`;
* `--results`: Caminho de um diretório alvo para gerar os arquivos resultantes do processamento

Exemplo de comando para execução do pipeline:

```bash
nextflow run main.nf -with-report --samples /home/watson/wilder/pharmacogenetics-analyses/nextflow/fastqs/NS20191022/ --reference /home/watson/wilder/reference-genomes/GRCh37/ --resources /home/watson/wilder/reference-genomes/b37-resource-bundle/ --results /home/watson/wilder/pharmacogenetics-analyses/nextflow/20190346/
```

#### Observações
* Para cada amostra processada, o pipeline produzirá um diretório (dentro do diretório especificado no parâmetro `--results`), e nesse diretório publicará os seguintes artefatos:

    * Um arquivo de alinhamento com seu respectivo índice (`<amostra>.sorted.duplicate_marked.recalibrated.bam` e `<amostra>.sorted.duplicate_marked.recalibrated.bai`);
    * Um arquivo de chamada de variantes e seu respectivo índice (`<amostra>.filtered.vcf` e `<amostra>.filtered.vcf.idx`);
    * Um arquivo tabulado com os resultados produzidos pelo Stargazer (`<amostra>.haplotypes.tsv`);
    * Um relatório em pdf (também produzido pelo Stargazer) com gráficos de informações de número de cópias dos genes analisados (`<amostra>.CNV-report.pdf`).

* O pipeline tenta fazer a definição de haplótipos para os 54 genes suportados pelo Stargazer (`CACNA1S`, `CFTR`, `CYP1A1`, `CYP1A2`, `CYP1B1`, `CYP2A6`, `CYP2A7`, `CYP2A13`, `CYP2B6`, `CYP2B7`, `CYP2C8`, `CYP2C9`, `CYP2C19`, `CYP2D6`, `CYP2D7`, `CYP2E1`, `CYP2F1`, `CYP2J2`, `CYP2R1`, `CYP2S1`, `CYP2W1`, `CYP3A4`, `CYP3A5`, `CYP3A7`, `CYP3A43`, `CYP4B1`, `CYP26A1`, `CYP4F2`, `CYP19A1`, `DPYD`, `G6PD`, `GSTM1`, `GSTP1`, `GSTT1`, `IFNL3`, `NAT1`, `NAT2`, `NUDT15`, `POR`, `RYR1`, `SLC15A2`, `SLC22A2`, `SLCO1B1`, `SLCO1B3`, `SLCO2B1`, `SULT1A1`, `TBXAS1`, `TPMT`, `UGT1A1`, `UGT1A4`, `UGT2B7`, `UGT2B15`, `UGT2B17` e `VKORC1`). Contudo, caso não haja cobertura em algum desses genes em uma dada amostra a execução do Stargazer encontrará um erro. O pipeline continuará a ser executado, mas no final não haverá informação para esse gene nessa amostra;

* O Nextflow criará um arquivo `.nextflow.log`, e dois diretórios `work` e `.nextflow` dentro do diretório de execução dos pipelines. Esses artefatos podem ser excluídos apenas quando o pipeline terminar o processamento com sucesso;

* Também será criado um arquivo `report.html` no diretório de execução do Nextflow. Esse relatório é importante para análise de métricas de processamento do pipeline, mas não necessariamente deverá ser armazenado;

* Caso haja um erro na execução, podemos consertá-lo e continuar a execução do pipeline de onde a mesma parou executando o mesmo comando inicial e adicionando a opção `-resume`.

----

### Considerações Finais

Esse pipeline foi desenvolvido pela equipe de bioinformática da Genomika. Para dúvidas, sugestões ou comentários entre em contato com algum integrante do time: Marcel (marcel@genomika.com.br), Rodrigo Bertollo (rodrigo@genomika.com.br), George (george@genomika.com.br) ou Wilder (wilder@genomika.com.br).
