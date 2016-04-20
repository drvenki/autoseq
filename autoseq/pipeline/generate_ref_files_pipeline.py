import json
import logging

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from pypedream.runners.shellrunner import Shellrunner

from autoseq.tools.genes import FilterGTFChromosomes, GTF2GenePred, FilterGTFGenes
from autoseq.tools.indexing import BwaIndex, SamtoolsFaidx, GenerateChrSizes
from autoseq.tools.intervals import SlopIntervalList, IntervalListToBed, MsiSensorScan, IntersectMsiSites
from autoseq.tools.picard import PicardCreateSequenceDictionary
from autoseq.tools.qc import *
from autoseq.tools.unix import Gunzip, Curl, Copy
from autoseq.tools.variantcalling import VcfFilter, CurlSplitAndLeftAlign, InstallVep
from autoseq.util.path import stripsuffix, normpath

__author__ = 'dankle'


class GenerateRefFilesPipeline(PypedreamPipeline):
    outdir = None
    maxcores = None

    def __init__(self, genome_resources, outdir, maxcores=1, runner=Shellrunner()):
        PypedreamPipeline.__init__(self, normpath(outdir), runner=runner)

        self.genome_resources = genome_resources
        self.input_reference_sequence = "{}/human_g1k_v37_decoy.fasta.gz".format(genome_resources)
        self.cosmic_vcf = "{}/CosmicCodingMuts_v71.vcf.gz".format(genome_resources)
        self.qdnaseq_background = "{}/qdnaseq_background.Rdata".format(genome_resources)
        self.outdir = outdir
        self.maxcores = maxcores
        self.reference_data = dict()

        self.exac_remote = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz"
	self.dbsnp_remote = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160408.vcf.gz"
        self.clinvar_remote = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive/2016/clinvar_20160203.vcf.gz"
        self.icgc_somatic_remote = "https://dcc.icgc.org/api/v1/download?fn=/release_20/Summary/simple_somatic_mutation.aggregated.vcf.gz"
        self.ypar_remote = "https://raw.githubusercontent.com/dakl/genome/master/ypars.bed"
        self.ensembl_version = "75"
        self.ensembl_gtf_remote = "ftp://ftp.ensembl.org/pub/release-" + self.ensembl_version + \
                                  "/gtf/homo_sapiens/Homo_sapiens.GRCh37." + self.ensembl_version + ".gtf.gz"
        self.mitranscriptome_remote = "http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz"

        self.prepare_reference_genome()
        self.prepare_genes()
        self.prepare_intervals()
        self.prepare_variants()

        install_vep = InstallVep()
        install_vep.ensembl_version = self.ensembl_version
        install_vep.output_cache_dir = "{}/vep/cache/".format(self.outdir)
        install_vep.output_bin_dir = "{}/vep/bin/".format(self.outdir)
        self.add(install_vep)

        self.reference_data['vep_bin_dir'] = install_vep.output_bin_dir
        self.reference_data['vep_cache_dir'] = install_vep.output_cache_dir

        with open("{}/autoseq-genome.json".format(self.outdir), "w") as output_file:
            json.dump(self.reference_data, output_file, indent=4, sort_keys=True)

    def prepare_variants(self):
        curl_dbsnp = CurlSplitAndLeftAlign()
        curl_dbsnp.input_reference_sequence = self.reference_data['reference_genome']
        curl_dbsnp.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_dbsnp.remote = self.dbsnp_remote
        curl_dbsnp.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.dbsnp_remote))
        curl_dbsnp.is_intermediate = True
        self.add(curl_dbsnp)

        filter_dbsnp = VcfFilter()
        filter_dbsnp.input = curl_dbsnp.output
        filter_dbsnp.filter = "\"! ( SAO = 3 | SAO = 2 )\""
        filter_dbsnp.output = "{}/variants/dbsnp142-germline-only.vcf.gz".format(self.outdir)
        self.add(filter_dbsnp)

        curl_cosmic = CurlSplitAndLeftAlign()
        curl_cosmic.input_reference_sequence = self.reference_data['reference_genome']
        curl_cosmic.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_cosmic.remote = "file://" + self.cosmic_vcf
        curl_cosmic.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.cosmic_vcf))
        self.add(curl_cosmic)

        curl_clinvar = CurlSplitAndLeftAlign()
        curl_clinvar.input_reference_sequence = self.reference_data['reference_genome']
        curl_clinvar.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_clinvar.remote = self.clinvar_remote
        curl_clinvar.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.clinvar_remote))
        self.add(curl_clinvar)

        curl_exac = CurlSplitAndLeftAlign()
        curl_exac.input_reference_sequence = self.reference_data['reference_genome']
        curl_exac.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_exac.remote = self.exac_remote
        curl_exac.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.exac_remote))
        self.add(curl_exac)

        curl_icgc = CurlSplitAndLeftAlign()
        curl_icgc.input_reference_sequence = self.reference_data['reference_genome']
        curl_icgc.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_icgc.remote = self.icgc_somatic_remote
        curl_icgc.output = "{}/variants/{}".format(self.outdir,
                                                   "icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz")
        self.add(curl_icgc)

        self.reference_data['dbSNP'] = filter_dbsnp.output
        self.reference_data['cosmic'] = curl_cosmic.output
        self.reference_data['exac'] = curl_exac.output
        self.reference_data['clinvar'] = curl_clinvar.output
        self.reference_data['icgc'] = curl_icgc.output

    def prepare_intervals(self):
        self.reference_data['targets'] = {}
        target_intervals_dir = "{}/target_intervals/".format(self.genome_resources)
        input_files = [f for f in os.listdir(target_intervals_dir) if f.endswith(".interval_list")]

        scan_for_microsatellites = MsiSensorScan()
        scan_for_microsatellites.input_fasta = self.reference_data['reference_genome']
        scan_for_microsatellites.homopolymers_only = True
        scan_for_microsatellites.output = "{}/intervals/msisensor-microsatellites.tsv".format(self.outdir)
        self.add(scan_for_microsatellites)

        for f in input_files:
            file_full_path = "{}/target_intervals/{}".format(self.genome_resources, f)
            logging.debug("Parsing intervals file {}".format(file_full_path))
            kit_name = stripsuffix(f, ".interval_list")
            self.reference_data['targets'][kit_name] = {}

            copy_file = Copy()
            copy_file.input = file_full_path
            copy_file.output = "{}/intervals/targets/{}".format(self.outdir, os.path.basename(file_full_path))
            self.add(copy_file)

            slop_interval_list = SlopIntervalList()
            slop_interval_list.input = copy_file.output
            slop_interval_list.output = stripsuffix(copy_file.output, ".interval_list") + ".slopped20.interval_list"
            self.add(slop_interval_list)

            interval_list_to_bed = IntervalListToBed()
            interval_list_to_bed.input = slop_interval_list.output
            interval_list_to_bed.output = stripsuffix(slop_interval_list.output, ".interval_list") + ".bed"
            self.add(interval_list_to_bed)

            intersect_msi = IntersectMsiSites()
            intersect_msi.input_msi_sites = scan_for_microsatellites.output
            intersect_msi.target_bed = interval_list_to_bed.output
            intersect_msi.output_msi_sites = stripsuffix(interval_list_to_bed.output, ".bed") + ".msisites.tsv"
            self.add(intersect_msi)

            cnvkit_ref_file = stripsuffix(file_full_path, ".interval_list") + ".cnn"
            if os.path.exists(cnvkit_ref_file):
                copy_cnvkit_ref = Copy()
                copy_cnvkit_ref.input = cnvkit_ref_file
                copy_cnvkit_ref.output = "{}/intervals/targets/{}".format(self.outdir,
                                                                          os.path.basename(cnvkit_ref_file))
                self.add(copy_cnvkit_ref)
                self.reference_data['targets'][kit_name]['cnvkit-ref'] = copy_cnvkit_ref.output
            else:
                self.reference_data['targets'][kit_name]['cnvkit-ref'] = None

            self.reference_data['targets'][kit_name]['targets-interval_list'] = copy_file.output
            self.reference_data['targets'][kit_name]['targets-interval_list-slopped20'] = slop_interval_list.output
            self.reference_data['targets'][kit_name]['targets-bed-slopped20'] = interval_list_to_bed.output
            self.reference_data['targets'][kit_name]['msisites'] = intersect_msi.output_msi_sites

    def prepare_genes(self):
        curl_ensembl_gtf = Curl()
        curl_ensembl_gtf.remote = self.ensembl_gtf_remote
        curl_ensembl_gtf.output = "{}/genes/{}".format(self.outdir, os.path.basename(self.ensembl_gtf_remote))
        curl_ensembl_gtf.jobname = "curl-ensembl-gtf"
        curl_ensembl_gtf.is_intermediate = True
        self.add(curl_ensembl_gtf)

        gunzip_ensembl_gtf = Gunzip()
        gunzip_ensembl_gtf.input = curl_ensembl_gtf.output
        gunzip_ensembl_gtf.output = stripsuffix(curl_ensembl_gtf.output, ".gz")
        gunzip_ensembl_gtf.is_intermediate = True
        self.add(gunzip_ensembl_gtf)

        filt_ensembl_gtf_chrs = FilterGTFChromosomes()
        filt_ensembl_gtf_chrs.input = gunzip_ensembl_gtf.output
        filt_ensembl_gtf_chrs.output = stripsuffix(gunzip_ensembl_gtf.output, ".gtf") + ".filtered.gtf"
        self.add(filt_ensembl_gtf_chrs)

        gtf2genepred_ensembl = GTF2GenePred()
        gtf2genepred_ensembl.input = filt_ensembl_gtf_chrs.output
        gtf2genepred_ensembl.output = stripsuffix(filt_ensembl_gtf_chrs.output, ".gtf") + ".genepred"
        self.add(gtf2genepred_ensembl)

        filt_genes_ensembl_gtf_genes = FilterGTFGenes()
        filt_genes_ensembl_gtf_genes.input = filt_ensembl_gtf_chrs.output
        filt_genes_ensembl_gtf_genes.output = stripsuffix(filt_ensembl_gtf_chrs.output, ".gtf") + ".genes-only.gtf"
        self.add(filt_genes_ensembl_gtf_genes)

        self.reference_data['ensemblVersion'] = self.ensembl_version
        self.reference_data['genesGtf'] = filt_ensembl_gtf_chrs.output
        self.reference_data['genesGenePred'] = gtf2genepred_ensembl.output
        self.reference_data['genesGtfGenesOnly'] = filt_genes_ensembl_gtf_genes.output

    def prepare_reference_genome(self):
        genome_unzipped = stripsuffix(os.path.basename(self.input_reference_sequence), ".gz")

        gunzip_ref = Gunzip()
        gunzip_ref.input = self.input_reference_sequence
        gunzip_ref.output = "{}/genome/{}".format(self.outdir, genome_unzipped)
        self.add(gunzip_ref)

        copy_ref_to_bwa = Copy()
        copy_ref_to_bwa.input = gunzip_ref.output
        copy_ref_to_bwa.output = "{}/bwa/{}".format(self.outdir, os.path.basename(gunzip_ref.output))
        self.add(copy_ref_to_bwa)

        bwa_index = BwaIndex()
        bwa_index.input_fasta = copy_ref_to_bwa.output
        bwa_index.output = copy_ref_to_bwa.output + ".bwt"
        bwa_index.algorithm = "bwtsw"
        self.add(bwa_index)

        create_dict = PicardCreateSequenceDictionary()
        create_dict.input = gunzip_ref.output
        create_dict.output_dict = gunzip_ref.output.replace(".fasta", "") + ".dict"
        self.add(create_dict)

        samtools_faidx = SamtoolsFaidx()
        samtools_faidx.input_fasta = gunzip_ref.output
        samtools_faidx.output = gunzip_ref.output + ".fai"
        self.add(samtools_faidx)

        create_chrsizes = GenerateChrSizes()
        create_chrsizes.input_fai = samtools_faidx.output
        create_chrsizes.output = gunzip_ref.output.replace(".fasta", "") + ".chrsizes.txt"
        self.add(create_chrsizes)

        copy_qdnaseq_bg = Copy()
        copy_qdnaseq_bg.input = self.qdnaseq_background
        copy_qdnaseq_bg.output = "{}/genome/{}".format(self.outdir, os.path.basename(self.qdnaseq_background))
        self.add(copy_qdnaseq_bg)

        self.reference_data['reference_genome'] = gunzip_ref.output
        self.reference_data['reference_dict'] = create_dict.output_dict
        self.reference_data['chrsizes'] = create_chrsizes.output
        self.reference_data['bwaIndex'] = bwa_index.input_fasta
        self.reference_data['qdnaseq_background'] = copy_qdnaseq_bg.output
