# autoseq

## Automated testing on travis-ci

For automated testing, a test reference genome and a test datas set with relevant data are supplied. 

### Reference genome

The test reference genome and assets is available for download at `https://export.uppmax.uu.se/b2010040/test-genome.tar.gz`. This archive contains a sliced version of a full set of genome files for autoseq, including various key genes. 

The whole chromosomes 3, 10, 17, X and Y are selected, after which everything except the following regions have been masked (to speed up alignment):

~~~
3	178863388	179014224	PIK3CA_150k
10	83068546	96283182	PTEN_13M
17	7558477	7589399	TP53_30k
X	66782057	66796840	14k_AR_exon
Y	6810425	6825985	15k_on_Y
~~~

From these regions, key exons and various other regions have been selected to mimic a small exome. 

### The Test Dataset

A sythetic tumor/normal/plasma dataset has been created for testing purpuses. From the illumina platinum 200x WGS sample from NA12877, read pairs from the seleted targets have been extracted. These reads have then been randomly assigned to create a virtual normal sample with ≈50x coverage, and remaining reads (≈150x coverage) have been put aside. To create a virtual tumor and a virtual plasma sample, variants have been spiked into the 150x data in the following positions: 

* TP53 insertion: MU2185182, chr17:g.7578475->G
* TP53 deletion: MU25947, chr17:g.7577558G>-
* TP53 DNV: MU52971976, chr17:g.7574003GG>AA
* PIK3CA hotspot E545K, MU5219, chr3:g.178936091G>A
* PTEN hotspot R130Q, MU29098, chr10:g.89692905G>A
* PTEN hotspot R233*, MU589331, chr10:g.89717672C>T
* AR intron variant, MU50988553, chrX:g.66788924G>A

In the virtual tumor, the target variant allele fraction (VAF) is 30% and in the virtual plasma sample the target VAF is 20%. 

The variants have been selected from ICGC simple somatic mutations v20 with the aim to cover common small variants, including SNVs, deletions, insertions and DNVs. Note that the tests does not address the issue of global sensitivity and PPV of the pipeline, but are only intented to ensure that variants of all kinds are detected by the pipeline. 
