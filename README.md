## Pre-processing scRNA-seq reads using kallisto|bustools

This repository contains the scripts I used for alignment and quantification of 10X Genomics Chromium scRNA-seq data using the kallisto|bustools pipeline.  

10X Genomics provide their own software - [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) - for pre-processing scRNA-seq data generated using one of their protocols. Cell Ranger is well-documented and really easy to run but it's also quite slow. Kallisto|bustools is an alternative that pseudo-aligns reads to equivalence classes (more details can be found in this [preprint](https://www.biorxiv.org/content/10.1101/673285v2) by Melsted, Booeshaghi et al., BioRxiv 2019), which really speeds up the process. According to the authors, their pipeline is "*up to 51 times faster than Cell Ranger*" and apparently it's even better for the [environment](https://twitter.com/lpachter/status/1217148183052111872?s=20)! The pipeline is modular which makes it easier to hack but the authors have also recently released [kb-python](https://github.com/pachterlab/kb_python) - a wrapper around the pipeline that really simplifies the workflow. For these reasons, I decided to use kallisto|bustools to pre-process my scRNA-seq reads but I encountered a few hurdles along the way which I've described here. Also, I found that the results that kallisto|bustools returns are a bit scarse, compared to Cell Ranger which returns too much information, so this is my Goldilocks solution!

<p align="center">

<img src="https://www.thenational.ae/image/policy:1.677154:1511147194/image.jpg?f=16x9&q=0.6&w=1200&$p$f$q$w=70c86c9" width="300">

</p>

-----

#### Step 1: Building a reference... 

##### ....(gtf files can GTFO)

The first step in the process is to build a reference transcriptome index for kallisto to pseudoalign to. This was probably the most painful step for several reasons. Firstly, you would expect that the cdna.all fasta file from Ensembl (which contains the sequences for all transcripts resulting from Ensembl gene predictions) would be a good place to start. Indeed, earlier versions of kallisto|bustools documentation and tutorials suggested using this file to build the index and then use the matching gtf file (which annotates ensembl transcript ids with lots of information) from Ensembl to build a t2g.txt file, which will then be used to map transcripts to gene symbols. The problem with this approach is that the cdna fasta file and gtf file from Ensembl don't exactly match. There are transcripts in the cdna file that aren't in the gtf and vice versa (more information in [this very helpful blogpost](https://fromsystosys.netlify.com/2020/01/31/comparing-ensembl-gtf-and-cdna/)). To get around this, the kb-python package has a handy `kb ref` command to build the index file. The `kb ref` command takes a gtf file and a dna fasta file as input. These are used to reconstruct a cdna.fa file that will contain all the transcripts in the input gtf file which can then be used to construct the transcriptome index. The t2g.txt file is also constructed from the input gtf file. In addition the `kb ref` command has a `-d` option to download a pre-built index and t2g.txt file for human or mouse transcriptomes. The files that are downloaded when running `kb ref -d human` are the exact same as what you would get if you used `kb ref` to build an index with the most recent Ensembl [gtf](ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz) and [dna.fa](ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz). The problem with using this as your index is that the Ensembl gtf contains annotations for ~60,700 genes, a large proportion of which are non-coding RNAs, pseudogenes and TEC (to be experimentally confirmed) genes. The gene biotypes that are included in the Ensembl gtf are shown in the screenshot below:

<p align="center">

<img src="https://github.com/Sarah145/scRNA_pre_process/blob/master/imgs/ens_biotypes.png?raw=true" width="500">

</p>

Given that my cells were sequenced with the 10X genomics protocol, the RNAs should be polyA selected and so shouldn't contain many of these gene biotypes. Including all these genes in the reference will confuse the alignment and you'll inevitably get reads falsely mapping to them so you end up with counts for over 60,000 genes - obviously not ideal. This prompted me to look into what Cell Ranger include in their reference, shown in this screenshot of the [Cell Ranger documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references):

<p align="center">

<img src="https://github.com/Sarah145/scRNA_pre_process/blob/master/imgs/cr_biotypes.png?raw=true" width="500">

</p>

The Cell Ranger reference contains much fewer biotypes but includes some biotypes that are not included in the Ensembl gtf, such as lincRNA (which is polyadenylated). The reference files that Cell Ranger uses for humans (hg38) can be downloaded from [their website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and the gtf file can be found within the downloaded folder:  `refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf` (there's also a copy of it in this repository: cellranger_genes.gtf.gz). This gtf contains annotations for ~33,500 genes which is much more reasonable so I decided to use this gtf to build my reference with `kb ref` (see `build_ref.sh` script in this repository).  Ideally, the `kb ref -d` option would allow you to download a pre-built index that was built using the Cell Ranger gtf but alas...

<p align="center">

<img src="https://github.com/Sarah145/scRNA_pre_process/blob/master/imgs/alas.gif?raw=true">

</p>

------

#### Step 2: Generate a count matrix



