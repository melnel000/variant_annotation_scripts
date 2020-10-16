### :desktop_computer: Ensembl Variant Effect Predictor (VEP):

This script runs [vep](https://www.ensembl.org/info/docs/tools/vep/index.html) to perform annotation of variants using a local cache (101_GRCh38) containing Ensembl transcripts.

--everything flag is used to add the following annotations: --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --var_synonyms, --variant_class, --mane

The following additional fields are added to the output file: --nearest symbol, --overlaps, --cell_type (optional, see masterscript variables below)

The input is a vcf file and the output is a tab delimited text file which facilitates filtering and display of data.

The script is configured to reduce each variant to a single row in the output file to simplify filtering. This is achieved by displaying only a single consequence for each variant using the most biologically relevant transcript. Flag "--pick_order" defines the order of priority by which various kinds of transcript support evidence are applied to select the transcript:

1. MANE transcript status
2. canonical status of transcript
3. APPRIS isoform annotation
4. transcript support level
5. biotype of transcript ("protein_coding" preferred)
6. CCDS status of transcript
7. consequence rank according to this table
8. translated, transcript or feature length (longer preferred)

MANE transcript status is prioritized and refers to the Matched Annotation from the NCBI and EMBL-EBI (MANE) project which aims to release a genome-wide transcript set that contains one well-supported transcript per protein-coding locus. All transcripts in the MANE set will perfectly align to GRCh38 and will represent 100% identity (5’UTR, coding sequence, 3’UTR) between the RefSeq (NM) and corresponding Ensembl (ENST) transcript.

### :earth_africa: Frequency data:

The following population frequency data is added:

dataset | sample information | population groups | data source
------------ | ------------- | ------------ | -------------
[1000 Genomes (phase 3) <sup> 1](https://www.internationalgenome.org/) | 2548 WGS samples across 26 populations from the 1000 Genomes Project | AFR (indigenous African, African American, AMR (admixed American), EAS (East Asian), EUR (European), SAS (South Asian) | vep cache
[gnomAD exomes (version 2.1) <sup> 2](https://gnomad.broadinstitute.org/) | 125,748  WES samples from various large scale sequencing projects | AFR (indigenous African, African American), AMR (admixed American), ASJ (Ashkenazi Jewish), EAS (East Asian), FIN (Finnish), NFE (non-Finnish European), SAS (South Asian), OTH (other) | vep cache
[NHLBI-ESP (V2-SSA137)](https://evs.gs.washington.edu/EVS/)| 6503 WES samples | African and European American populations | vep cache
[gnomAD genomes (version 3)](https://gnomad.broadinstitute.org/blog/2019-10-gnomad-v3-0/) | 71,702 WGS samples | AMR (admixed American), NFE (non-Finnish European), EAS (East Asian), FIN (Finnish), ASJ (Ashkenazi Jewish), SAS (South Asian), AFR (indigenous African, African American), OTH (other), AMI (Amish) | custom annotation

1. Lowy-Gallego et al., 2019. *Variant calling on the GRCh38 assembly with the data from phase three of the 1000 Genomes Project*.
2. Karczewski et al., 2019. *Variation across 141,456 human exomes and genomes reveals the spectrum of loss-of-function intolerance across human protein-coding genes*.

### :desktop_computer: Plugins:

The script is configured to run the following plugins:

#### <ins> Function prediction tools for splice region variants <ins>

Tool | Annotation 
------------ | ------------- 
[SpliceRegion](https://www.biorxiv.org/content/10.1101/256636v2.full) | Adds the following consequences for non-canonical positions in splice sites: splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant
[GeneSplicer](https://academic.oup.com/nar/article/29/5/1185/2381159) | Predicted splicing regions that overlap the variant are reported with the following states: no_change, diff, gain, loss including score and confidence level (low, medium,high)
[MaxEntScan (MES)](https://pubmed.ncbi.nlm.nih.gov/30475984/) | Splice site predictions based on reference, alternate and difference (REF-ALT) maximum entropy scores
[SpliceAI](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf) | Deep neural network which predicts splice altering variants based on probability scores
[dbscSNV](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267638/) | Ensemble score based on 8 in silico tools: PWM, MES, NNSplice, HSF, GeneSplicer, GENSCAN, NetGene2, SplicePredictor

#### <ins> Function prediction tools for missense variants <ins>

<ins> Masterscript variables: </ins>

$1: path to the working directory (without trailing slash) <br>
$2: naming prefix for output files <br>
$3: path to the input vcf file <br>
$4: --cell_type followed by comma separated list of cell type(s) [listed here](/Pipeline/Files/vep_cell_types) to report only regulatory regions that are found in the given cell types <br>
$5: path to gnomADv3 genomes file
