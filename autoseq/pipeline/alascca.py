import logging

from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import CNVkit, AlasccaCNAPlot
from autoseq.tools.contamination import ContEst, ContEstToContamCaveat
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.picard import PicardCollectHsMetrics
from autoseq.tools.picard import PicardCollectInsertSizeMetrics
from autoseq.tools.picard import PicardCollectOxoGMetrics
from autoseq.tools.qc import *
from autoseq.tools.reports import CompileMetadata, CompileAlasccaGenomicJson, WriteAlasccaReport
from autoseq.tools.variantcalling import Freebayes, VcfAddSample, call_somatic_variants
from autoseq.util.library import get_libdict, find_fastqs
from autoseq.util.path import normpath, stripsuffix

__author__ = 'dankle'


class AlasccaPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, outdir, libdir, analysis_id=None, maxcores=1, scratch="/scratch/tmp/tmp/",
                 referral_db_conf="tests/referrals/referral-db-config.json",
                 addresses="tests/referrals/addresses.csv",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, outdir, libdir, analysis_id,
                                 maxcores, scratch, **kwargs)

        logging.debug("Unnormalized outdir is {}".format(outdir))
        logging.debug("self.outdir is {}".format(self.outdir))
        self.referral_db_conf = referral_db_conf
        self.addresses = addresses
        self.targets_name = get_libdict(self.sampledata['panel']['T'])['capture_kit_name']
        self.panel_tumor_fqs = find_fastqs(self.sampledata['panel']['T'], self.libdir)
        self.panel_normal_fqs = find_fastqs(self.sampledata['panel']['N'], self.libdir)

        panel_bams = self.analyze_panel()

        ################################################
        # QC

        # per-bam qc
        self.qc_files += self.run_panel_bam_qc([panel_bams['tbam'], panel_bams['nbam']])

        # per-fastq qc
        fqs = self.get_all_fastqs()
        self.qc_files += self.run_fastq_qc(fqs)

        multiqc = MultiQC()
        multiqc.input_files = self.qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.analysis_id)
        multiqc.jobname = "multiqc/{}-{}".format(self.sampledata['panel']['T'],
                                                 self.sampledata['panel']['N'])
        self.add(multiqc)

    def get_all_fastqs(self):
        fqs = []
        if self.sampledata['panel']['T']:
            fqs.extend(self.panel_tumor_fqs[0])
            fqs.extend(self.panel_tumor_fqs[1])
        if self.sampledata['panel']['N']:
            fqs.extend(self.panel_normal_fqs[0])
            fqs.extend(self.panel_normal_fqs[1])

        return [fq for fq in fqs if fq is not None]

    def analyze_panel(self):
        if self.sampledata['panel']['T'] is None or self.sampledata['panel']['T'] == "NA":
            logging.info("No panel data found.")
            return {}

        tbam = align_library(self,
                             fq1_files=self.panel_tumor_fqs[0],
                             fq2_files=self.panel_tumor_fqs[1],
                             lib=self.sampledata['panel']['T'],
                             ref=self.refdata['bwaIndex'],
                             outdir=self.outdir + "/bams/panel",
                             maxcores=self.maxcores)

        nbam = align_library(self,
                             fq1_files=self.panel_normal_fqs[0],
                             fq2_files=self.panel_normal_fqs[1],
                             lib=self.sampledata['panel']['N'],
                             ref=self.refdata['bwaIndex'],
                             outdir=self.outdir + "/bams/panel",
                             maxcores=self.maxcores)

        vep = False
        if self.refdata['vep_dir']:
            vep = True

        somatic_vcfs = call_somatic_variants(self, tbam, nbam,
                                             tlib=self.sampledata['panel']['T'],
                                             nlib=self.sampledata['panel']['N'],
                                             target_name=self.targets_name,
                                             refdata=self.refdata,
                                             outdir=self.outdir,
                                             callers=['vardict'],
                                             vep=vep,
                                             min_alt_frac=0.02)

        germline_vcf = self.call_germline_variants(nbam, library=self.sampledata['panel']['N'])

        libdict = get_libdict(self.sampledata['panel']['N'])
        rg_sm = "{}-{}-{}".format(libdict['sdid'], libdict['type'], libdict['sample_id'])

        hzconcordance = HeterzygoteConcordance()
        hzconcordance.input_vcf = germline_vcf
        hzconcordance.input_bam = tbam
        hzconcordance.reference_sequence = self.refdata['reference_genome']
        hzconcordance.target_regions = self.refdata['targets'][self.targets_name]['targets-interval_list-slopped20']
        hzconcordance.normalid = rg_sm
        hzconcordance.filter_reads_with_N_cigar = True
        hzconcordance.jobname = "hzconcordance/{}".format(self.sampledata['panel']['T'])
        hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(self.outdir, self.sampledata['panel']['T'],
                                                                        self.sampledata['panel']['N'])
        self.add(hzconcordance)
        self.qc_files.append(hzconcordance.output)

        vcfaddsample = VcfAddSample()
        vcfaddsample.input_bam = tbam
        vcfaddsample.input_vcf = germline_vcf
        vcfaddsample.samplename = self.sampledata['panel']['T']
        vcfaddsample.filter_hom = True
        vcfaddsample.scratch = self.scratch
        vcfaddsample.output = "{}/variants/{}-and-{}.germline.vcf.gz".format(self.outdir,
                                                                             self.sampledata['panel']['T'],
                                                                             self.sampledata['panel']['N'])
        vcfaddsample.jobname = "vcf-add-sample/{}".format(self.sampledata['panel']['T'])
        self.add(vcfaddsample)

        msisensor = MsiSensor()
        msisensor.msi_sites = self.refdata['targets'][self.targets_name]['msisites']
        msisensor.input_normal_bam = nbam
        msisensor.input_tumor_bam = tbam
        msisensor.output = "{}/msisensor.tsv".format(self.outdir)
        msisensor.threads = self.maxcores
        msisensor.scratch = self.scratch
        msisensor.jobname = "msisensor/{}".format(self.sampledata['panel']['T'])
        self.add(msisensor)

        cnvkit = CNVkit(input_bam=tbam,
                        output_cnr="{}/cnv/{}.cnr".format(self.outdir,
                                                          self.sampledata['panel']['T']),
                        output_cns="{}/cnv/{}.cns".format(self.outdir,
                                                          self.sampledata['panel']['T']),
                        scratch=self.scratch
                        )
        # If we have a CNVkit reference
        if self.refdata['targets'][self.targets_name]['cnvkit-ref']:
            cnvkit.reference = self.refdata['targets'][self.targets_name]['cnvkit-ref']
        else:
            cnvkit.targets_bed = self.refdata['targets'][self.targets_name]['targets-bed-slopped20']

        cnvkit.jobname = "cnvkit/{}".format(self.sampledata['panel']['T'])
        self.add(cnvkit)

        # TODO: Make function that assigns the correct inputs to ContEst() to avoid repetition of code for both T & N? Or keep like this?
        contest_tumor = ContEst()
        contest_tumor.reference_genome = self.refdata['reference_genome']
        contest_tumor.input_eval_bam = tbam
        contest_tumor.input_genotype_bam = nbam
        contest_tumor.population_af_vcf = self.get_pop_af_vcf()
        # TODO: Is it necessary to create the output subdir contamination somewhere? Check how it's done for e.g. cnvkit.
        contest_tumor.output = "{}/contamination/{}.contest.txt".format(self.outdir, self.sampledata['panel']['T'])  # TODO: Should the analysis id also be in name of out file?
        contest_tumor.jobname = "contest_tumor/{}".format(self.sampledata['panel']['T'])  # TODO: Is it ok that the job name does not contain analysis id, i.e. may not be unique?
        # only run the job if a population allele frequency vcf is implemented for the capture kits used for T & N:
        if contest_tumor.population_af_vcf:
            self.add(contest_tumor)

        contest_normal = ContEst()
        contest_normal.reference_genome = self.refdata['reference_genome']
        contest_normal.input_eval_bam = nbam
        contest_normal.input_genotype_bam = tbam
        contest_normal.population_af_vcf = self.get_pop_af_vcf()
        contest_normal.output = "{}/contamination/{}.contest.txt".format(self.outdir, self.sampledata['panel']['N'])  # Should the analysis id also be in name of out file?
        contest_normal.jobname = "contest_normal/{}".format(self.sampledata['panel']['N']) #Is it ok that the job name does not contain analysis id, i.e. may not be unique?
        # only run the job if a population allele frequency vcf is implemented for the capture kits used for T & N:
        if contest_normal.population_af_vcf:
            self.add(contest_normal)

        # Generate ContEst contamination QC call JSON files from the ContEst
        # outputs:
        process_contest_tumor = ContEstToContamCaveat()
        process_contest_tumor.input_contest_results = contest_tumor.output
        process_contest_tumor.output = "{}/qc/{}-contam-qc-call.json".format(self.outdir, self.sampledata['panel']['T'])
        if contest_tumor.population_af_vcf:
            # Only add the contest output processing if contest is to be run
            # for the tumor sample:
            self.add(process_contest_tumor)

        alascca_cna = AlasccaCNAPlot()
        alascca_cna.input_somatic_vcf = somatic_vcfs['vardict']
        alascca_cna.input_germline_vcf = vcfaddsample.output
        alascca_cna.input_cnr = cnvkit.output_cnr
        alascca_cna.input_cns = cnvkit.output_cns
        alascca_cna.chrsizes = self.refdata['chrsizes']
        alascca_cna.output_cna = "{}/variants/{}-alascca-cna.json".format(self.outdir,
                                                                          self.sampledata['panel']['T'])
        alascca_cna.output_purity = "{}/variants/{}-alascca-purity.json".format(self.outdir,
                                                                                self.sampledata['panel']['T'])
        alascca_cna.output_png = "{}/qc/{}-alascca-cna.png".format(self.outdir, self.sampledata['panel']['T'])
        alascca_cna.jobname = "alascca-cna/{}".format(self.sampledata['panel']['T'])
        self.add(alascca_cna)

        tlib = get_libdict(self.sampledata['panel']['T'])
        nlib = get_libdict(self.sampledata['panel']['N'])
        blood_barcode = nlib['sample_id']
        tumor_barcode = tlib['sample_id']
        metadata_json = "{}/report/{}-{}.metadata.json".format(self.outdir, blood_barcode, tumor_barcode)
        compile_metadata_json = CompileMetadata(self.referral_db_conf, blood_barcode, tumor_barcode,
                                                output_json=metadata_json,
                                                addresses=self.addresses)
        compile_metadata_json.jobname = "compile-metadata/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_metadata_json)

        genomic_json = "{}/report/{}-{}.genomic.json".format(self.outdir, blood_barcode, tumor_barcode)
        # FIXME: This is getting quite nasty (linking outputs and inputs, some conditional and others
        # generated in a different scope). Not sure how Daniel intended this sort of thing
        # to be done...
        contam_qc = None
        if contest_tumor.population_af_vcf:
            contam_qc = process_contest_tumor.output
        tumor_prefix = stripsuffix(os.path.basename(tbam), ".bam")
        normal_prefix = stripsuffix(os.path.basename(nbam), ".bam")
        tcov_qc = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, tumor_prefix)
        ncov_qc = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, normal_prefix)

        compile_genomic_json = CompileAlasccaGenomicJson(input_somatic_vcf=somatic_vcfs['vardict'],
                                                         input_cn_calls=alascca_cna.output_cna,
                                                         input_msisensor=msisensor.output,
                                                         input_purity_qc=alascca_cna.output_purity,
                                                         input_contam_qc=contam_qc,
                                                         input_tcov_qc=tcov_qc,
                                                         input_ncov_qc=ncov_qc,
                                                         output_json=genomic_json)
        compile_genomic_json.jobname = "compile-genomic/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_genomic_json)

        pdf = "{}/report/AlasccaReport-{}-{}.pdf".format(self.outdir, blood_barcode, tumor_barcode)
        writeAlasccaPdf = WriteAlasccaReport(input_genomic_json=compile_genomic_json.output_json,
                                             input_metadata_json=compile_metadata_json.output_json,
                                             output_pdf=pdf)
        writeAlasccaPdf.jobname = "writeAlasccaPdf/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(writeAlasccaPdf)

        return {'tbam': tbam, 'nbam': nbam}

    def call_germline_variants(self, bam, library):
        """
        Call germline variants from a bam
        :param bam:
        :return:
        """
        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.somatic_only = False
        freebayes.params = None
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][self.targets_name]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.output = "{}/variants/{}.freebayes-germline.vcf.gz".format(self.outdir, library)
        freebayes.jobname = "freebayes-germline/{}".format(library)
        self.add(freebayes)
        return freebayes.output

    def run_fastq_qc(self, fastq_files):
        """
        Run QC on fastq files
        :param fastq_files:
        :return:
        """
        qc_files = []
        for fq in fastq_files:
            basefn = stripsuffix(os.path.basename(fq), ".fastq.gz")
            fastqc = FastQC()
            fastqc.input = fq
            fastqc.outdir = "{}/qc/fastqc/".format(self.outdir)
            fastqc.output = "{}/qc/fastqc/{}_fastqc.zip".format(self.outdir, basefn)
            fastqc.jobname = "fastqc/{}".format(basefn)
            qc_files.append(fastqc.output)
            self.add(fastqc)
        return qc_files

    def run_panel_bam_qc(self, bams):
        """
        Run QC on panel bams
        :param bams: list of bams
        :return: list of generated files
        """

        qc_files = []
        for bam in bams:
            logging.debug("Adding QC jobs for {}".format(bam))
            basefn = stripsuffix(os.path.basename(bam), ".bam")
            isize = PicardCollectInsertSizeMetrics()
            isize.input = bam
            isize.output_metrics = "{}/qc/picard/panel/{}.picard-insertsize.txt".format(self.outdir, basefn)
            isize.jobname = "picard-isize/{}".format(basefn)
            self.add(isize)

            oxog = PicardCollectOxoGMetrics()
            oxog.input = bam
            oxog.reference_sequence = self.refdata['reference_genome']
            oxog.output_metrics = "{}/qc/picard/panel/{}.picard-oxog.txt".format(self.outdir, basefn)
            oxog.jobname = "picard-oxog/{}".format(basefn)
            self.add(oxog)

            hsmetrics = PicardCollectHsMetrics()
            hsmetrics.input = bam
            hsmetrics.reference_sequence = self.refdata['reference_genome']
            hsmetrics.target_regions = self.refdata['targets'][self.targets_name][
                'targets-interval_list-slopped20']
            hsmetrics.bait_regions = self.refdata['targets'][self.targets_name][
                'targets-interval_list-slopped20']
            hsmetrics.bait_name = self.targets_name
            hsmetrics.output_metrics = "{}/qc/picard/panel/{}.picard-hsmetrics.txt".format(self.outdir, basefn)
            hsmetrics.jobname = "picard-hsmetrics/{}".format(basefn)
            self.add(hsmetrics)

            sambamba = SambambaDepth()
            sambamba.targets_bed = self.refdata['targets'][self.targets_name]['targets-bed-slopped20']
            sambamba.input = bam
            sambamba.output = "{}/qc/sambamba/{}.sambamba-depth-targets.txt".format(self.outdir, basefn)
            sambamba.jobname = "sambamba-depth/{}".format(basefn)
            self.add(sambamba)

            alascca_coverage_hist = CoverageHistogram()
            if 'alascca_targets' in self.refdata['targets']:
                alascca_coverage_hist.input_bed = self.refdata['targets']['alascca_targets']['targets-bed-slopped20']
            else:
                alascca_coverage_hist.input_bed = self.refdata['targets'][self.targets_name][
                    'targets-bed-slopped20']
            alascca_coverage_hist.input_bam = bam
            alascca_coverage_hist.output = "{}/qc/{}.coverage-histogram.txt".format(self.outdir, basefn)
            alascca_coverage_hist.jobname = "alascca-coverage-hist/{}".format(basefn)
            self.add(alascca_coverage_hist)

            coverage_qc_call = CoverageCaveat()
            coverage_qc_call.input_histogram = alascca_coverage_hist.output
            coverage_qc_call.output = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, basefn)
            coverage_qc_call.jobname = "coverage-qc-call/{}".format(basefn)
            self.add(coverage_qc_call)

            qc_files += [isize.output_metrics, oxog.output_metrics,
                         hsmetrics.output_metrics, sambamba.output, alascca_coverage_hist.output,
                         coverage_qc_call.output]

        return qc_files

    def get_pop_af_vcf(self):
        """Get the path to the population allele frequency vcf to use in ContEst.
        Currently only available for panels CB (big design), CS (clinseq v3) & CZ (clinseq v4)"""

        tpanel = get_libdict(self.sampledata['panel']['T'])['capture_kit_name']
        npanel = get_libdict(self.sampledata['panel']['N'])['capture_kit_name']

        # TODO: Move this selection of which vcf to use somewhere else in the pipeline? Maybe possible to specify on command line?
        if tpanel == "big_design" and npanel == "big_design":
            return self.refdata['contest_vcfs']['big'] # "path/to/big_swegene_contest.vcf"
        elif tpanel in ["clinseq_v3_targets", "clinseq_v4"] and npanel in ["clinseq_v3_targets", "clinseq_v4"]:
            return self.refdata['contest_vcfs']['clinseqV3V4'] # "path/to/clinseqV3V4_exac_contest.vcf"
        elif tpanel in ["big_design", "clinseq_v3_targets", "clinseq_v4"] and npanel in ["big_design", "clinseq_v3_targets", "clinseq_v4"]:
            return self.refdata['contest_vcfs']['clinseqV3V4big'] # "path/to/clinseqV3V4big_intersection_exac_contest.vcf"
        elif tpanel == "test-regions" and npanel == "test-regions":
            return self.refdata['contest_vcfs']['test-regions'] # "path/to/test-regions_contest.vcf"
        else:
            raise ValueError("Invalid tpanel/npanel combination: {}/{}".format(tpanel, npanel))
