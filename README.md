## Pre-processing scRNA-seq reads using kallisto|bustools

This repository contains the scripts I used for alignment and quantification of 10X Genomics Chromium scRNA-seq data using the kallisto|bustools pipeline.  

10X Genomics provide their own software - [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) - for pre-processing scRNA-seq data generated using one of their protocols. Cell Ranger is well-documented and really easy to run but it's also quite slow. Kallisto|bustools is an alternative that pseudo-aligns reads to equivalence classes (more details can be found in this [preprint](https://www.biorxiv.org/content/10.1101/673285v2) by Melsted, Booeshaghi et al., BioRxiv 2019), which really speeds up the process. According to the authors, their pipeline is "*up to 51 times faster than Cell Ranger*" and apparently it's even [better for the environment](https://twitter.com/lpachter/status/1217148183052111872?s=20)! The pipeline is modular which makes it easier to hack but the authors have also recently released [kb-python](https://github.com/pachterlab/kb_python) - a wrapper around the pipeline that really simplifies the workflow. For these reasons, I decided to use kallisto|bustools to pre-process my scRNA-seq reads but I encountered a few hurdles along the way which I've described here. Also, I found that the results that kallisto|bustools returns are a bit scarse, compared to Cell Ranger which returns too much information, so this is my Goldilocks solution!

<p align="center">

<img src="https://www.thenational.ae/image/policy:1.677154:1511147194/image.jpg?f=16x9&q=0.6&w=1200&$p$f$q$w=70c86c9" width="300">

</p>

<sub>**Note**: To reproduce this analysis and skip my rant about how I got here, here's the [TL;DR](https://github.com/Sarah145/scRNA_pre_process#to-reproduce-this-analysis).</sub>

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

The Cell Ranger reference contains much fewer biotypes but includes some biotypes that are not included in the Ensembl gtf, such as lincRNA (which is polyadenylated). The reference files that Cell Ranger uses for humans (hg38) can be downloaded from [their website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and the gtf file can be found within the downloaded folder:  `refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf` (there's also a copy of it in this repository: cellranger_genes.gtf.gz). This gtf contains annotations for ~33,500 genes which is much more reasonable so I decided to use this gtf to build my reference with `kb ref` (see [build_ref.sh](https://github.com/Sarah145/scRNA_pre_process/blob/master/scripts/build_ref.sh) script in this repository).  Ideally, the `kb ref -d` option would allow you to download a pre-built index that was built using the Cell Ranger gtf but alas...

<p align="center">

<img src="https://github.com/Sarah145/scRNA_pre_process/blob/master/imgs/alas.gif?raw=true">

</p>

------

#### Step 2: Generate a raw count matrix

I found the best way to generate a raw count matrix is to use `kb count`, which (from what I can tell) is a wrapper around `kallisto bus`, `bustools correct`, `bustools sort` and `bustools count`. Briefly, this generates an output.bus file which contains a binary representation of equivalence class counts for barcodes and UMIs (more details in [this paper](https://academic.oup.com/bioinformatics/article/35/21/4472/5487510)). The barcodes are then corrected against a whitelist of valid barcodes to account for single-base sequencing errors in the reads. Then the bus file is sorted and `bustools count` generates three files: one with the count matrix, one with the gene names and one with the barcodes. After running `kb count`, the results directory will look something like this:

```bash
Sample1_kb_out/
├── counts_unfiltered
│   ├── cells_x_genes
│   ├── cells_x_genes.barcodes.txt
│   ├── cells_x_genes.genes.txt
│   └── cells_x_genes.mtx
├── inspect.json
├── matrix.ec
├── output.bus
├── output.unfiltered.bus
├── run_info.json
└── transcripts.txt
```

The cells_x_genes.mtx file is (as the name would suggest) a sparse matrix of counts with cells as rows and genes as columns. The default output from Cell Ranger is the transpose of this matrix (i.e. genes x cells) and since many of the downstream tools I'll be using expect 10X data in this format, I decided to reformat the output from `kb count` to make life easier for myself. Also the cells_x_genes.genes.txt file that `kb count` outputs contains only Ensembl gene IDs and not gene symbols which is not ideal. My [kb_count.sh](https://github.com/Sarah145/scRNA_pre_process/blob/master/scripts/kb_count.sh) script runs `kb_count` and then runs a [reformat.R](https://github.com/Sarah145/scRNA_pre_process/blob/master/scripts/reformat.R) script which mainly uses the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to reformat the `kb count` output to a (genes x cells) matrix.mtx file and a genes.tsv file with gene symbols (which it gets from the t2g.txt file) instead of Ensembl gene IDs.  

<sub>**Note**:  When running `kb count` you have to specify which sequencing technology (e.g. 10xv2, 10xv3) was used to sequence the cells - I've included code in my script that will automatically detect which technology was used based on the length of the R1 (26bp for v2, 28bp for v3), but this will only work for 10X reads. </sub>

------

#### Step 3: Filter the raw count matrix

The count matrix that was generated in the previous step is a raw count matrix, meaning that many of the 'cells' in this file are probably not cells but correspond to empty droplets where ambient RNA in the input cell suspension has been captured in a droplet and tagged with a barcode. Barcodes from these 'cells' will have a very low number of counts associated with them so Cell Ranger filters them out by determining an inflection point in the number of counts per barcode, above which barcodes are assumed to be tagging actual cells and not empty droplets.

-----

### To reproduce this analysis...

1. Clone this repository and navigate into it:

   ```bash
   git clone https://github.com/Sarah145/scRNA_pre_process
   cd scRNA_pre_process
   ```

2. Create conda environment (assuming you have [Anaconda](https://www.anaconda.com/distribution/#download-section) installed) from the [scRNA_pre_process.yml](https://github.com/Sarah145/scRNA_pre_process/blob/master/scRNA_pre_process.yml) file and activate it:

   ```bash
   conda env create -f scRNA_pre_process.yml
   conda activate scRNA_pre_process
   ```

   <sub> **Note:** This will take several minutes. </sub>

3. Download all necessary files:

   ```bash
   # reference files
   cd ref
   wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
   
   # data files 
   cd data
   # copy in your own fastq files
   cp /path/to/your/fastqs/*.fastq.gz . 
   		## OR ## 
   # download sample dataset from 10X (may take a while)
   curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
   tar -xf pbmc_1k_v3_fastqs.tar
   
   cd ..
   ```

4. Build the transcriptome index for pseudoalignment:

   ```bash
   cd scripts
   ./build_ref.sh
   ```

5. Open the [kb_count.sh](https://github.com/Sarah145/scRNA_pre_process/blob/master/scripts/kb_count.sh) script (using nano or vim or whatever floats your boat) and edit the first two lines with sample information and the path to your fastq files. Make sure to specify your fastq files in R1/R2 pairs - example:

   ```
   # Assign a sample_id
   sample_id="Sample1"
   
   # Point to fastq files in R1/R2 pairs
   fastqs=("../data/Sample_R1.fastq.gz" "../data/Sample1_L1_R2.fastq.gz" "../data/Sample1_L2_R1.fastq.gz" "../data/Sample1_L2_R2.fastq.gz") 
   
   ```

   