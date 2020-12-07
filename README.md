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

A partir da versão 2.0 esse pipeline se disponibiliza para ser executado em dois ambientes: *cloud* (usando o [Amazon Web Services](https://aws.amazon.com/)) e *on-premise*.

Para executar esse pipeline em um ambiente de execução local são necessários os seguintes requisitos:

* Nextflow instalado (versão recomendada `20.04.1.5335`) ([guia de instalação](https://www.nextflow.io/index.html#GetStarted));
* Docker instalado ([guia de instalação](https://www.docker.com/get-started));
* Diretório contendo o genoma de referência GRCh37, além de seus respectivos arquivos de índice e dicionário;
    * Esse diretório deve estar especificado no arquivo de configuração do pipeline (`nextflow.config`), no escopo `watson` e parametro `reference`
* Diretório com o pacote de recursos do genoma `b37` da Broad Institute.
    * Esse diretório deve estar especificado no arquivo de configuração do pipeline (`nextflow.config`), no escopo `watson` e parametro `resources`

Para executar esse pipeline na AWS deve-se preparar o ambiente de acordo com [a documentação oficial do Nextflow](https://www.nextflow.io/docs/latest/awscloud.html). Além desse ambiente é necessário disponibilizar as imagens Docker presentes no diretório `dockerfiles` em algum repositório remoto. Elas já se encontram disponíveis no [Elastic Container Registry](https://console.aws.amazon.com/ecr/repositories?region=us-east-1) da Genomika, em caso de mudanças seguir os passos abaixo e, se necessário, trocar nos processos do pipeline para a imagem atualizada.

#### Preparando as imagens Docker

Para construir a imagem necessária para executar o processo de alinhamento (usando o [`bwa`](https://github.com/lh3/bwa) e [`samtools`](https://github.com/samtools/samtools)), vá para o diretório `dockerfiles/alignment` execute o comando:

```bash
$ docker build . -t alignment:v0.1.1
```

Para construir a imagem necessária para executar o [`Picard`](https://broadinstitute.github.io/picard/), execute o comando:

```bash
$ docker pull broadinstitute/picard:2.22.8
```

Para construir a imagem necessária para executar o [GATK 3](https://github.com/broadgsa/gatk/releases/tag/3.8-1), execute o comando:

```bash
$ docker pull broadinstitute/gatk3:3.8-0
```

Para construir a imagem necessária para executar o [`mosdepth`](https://github.com/brentp/mosdepth), execute o comando:

```bash
$ docker pull quay.io/biocontainers/mosdepth:0.2.4--he527e40_0
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

Esses passos deixarão a imagem local, mas é necessário que essas imagens também estejam disponíveis no ECR da Genomika para que não seja mitigada a possibilidade de execução do pipeline no AWS. Com as imagens locais, veja [documentação oficial do ECR](https://docs.aws.amazon.com/AmazonECR/latest/userguide/docker-push-ecr-image.html) como fazer o *`push`* dessas imagens para a nuvem.

O passo acima só é necessário para as imagens que foram construidas localmente e não estão disponíveis em repositório remoto.

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
* `--threads`: Número de threads usadas no processo de alinhamento (recomendação - local: 4, AWS: 8);
* `--sdtype`: Tipo de sequenciamento usado nas amostras ('ts' para paineis e exoma 'wgs' para genomas);
* `--reference`: Caminho do diretório onde se encontra o `fasta` do genoma de referência, além dos arquivos de dicionário e índice;
* `--resources`: Caminho do diretório onde se encontram os arquivos do pacote de recursos criado pela Broad do genoma `b37`;
* `--results`: Caminho de um diretório alvo para gerar os arquivos resultantes do processamento

Exemplo de comando para execução do pipeline para execução local:

```bash
nextflow run main.nf -profile watson --samples /home/watson/wilder/pharmacogenetics-analyses/nextflow/fastqs/NS20191022/ --sdtype ts --threads 4 -with-report
```

Exemplo de comando para execução do pipeline para execução no AWS:

```bash
nextflow run main.nf -profile amazon --samples s3://genomika-samples-data/NS20191022/ -bucket-dir s3://genomika-pharmacogenetics-pipeline/work/ -with-report --threads 4 --sdtype ts
```

#### Observações

* Para executar o pipeline no AWS é necessário anteriormente fazer o upload dos fastqs das amostras para um diretório no bucket do S3 criado especificamente para esse propósito (`genomika-samples-data`);

```
aws s3 cp --recursive /home/watson/wilder/NS20191022 s3://genomika-samples-data/NS20191022
```

> As amostras ficam no `s3` do console, na pasta `genomika-samples-data`, a senha e login estão no email "Usuário do AWS".

* Além disso, para executar o pipeline no AWS também temos que especificar o parâmetro adicional `-bucket-dir` que será usado pelo Nextflow e pode ser usado de acordo com o exemplo;

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
