from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.cnvcalling import LiqbioCNAPlot
from autoseq.util.clinseq_barcode import *

__author__ = 'dankle'


class LiqBioPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, job_params, outdir, libdir, maxcores=1, scratch="/scratch/tmp/tmp",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, job_params, outdir, libdir,
                                 maxcores, scratch, **kwargs)

        # Set the min alt frac value:
        self.default_job_params["vardict-min-alt-frac"] = 0.01
        self.default_job_params["vardict-min-num-reads"] = None
        self.default_job_params["vep-additional-options"] = " --pick --filter_common "

        # Remove clinseq barcodes for which data is not available:
        self.check_sampledata()

        # Configure alignment and merging of fastq data for all clinseq barcodes:
        self.configure_align_and_merge()

        # Configure all panel analyses:
        self.configure_panel_analyses()

        # Configure liqbio-specific panel analyses:
        self.configure_panel_analyses_liqbio()

        # Configure additional msings analysis:
        self.configure_panel_msings_analyses()

        # Configure QC of all panel data:
        self.configure_all_panel_qcs()

        # Configure fastq QCs:
        self.configure_fastq_qcs()

        # Configure the low-pass whole genome analysis:
        self.configure_lowpass_analyses()

        # Configure low-pass whole genome data QC:
        self.configure_all_lowpass_qcs()

        # Configure MultiQC:
        self.configure_multi_qc()

    def configure_panel_analyses_liqbio(self):
        for normal_capture in self.get_mapped_captures_normal():
            for cancer_capture in self.get_mapped_captures_cancer():
                self.configure_panel_analysis_cancer_vs_normal_liqbio(
                    normal_capture, cancer_capture)

    def configure_liqbio_cna(self, normal_capture, cancer_capture):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)]
        tumor_results = self.capture_to_results[cancer_capture]

        # Configure the liqbio frankenplots:
        # NOTE: Get PureCN outputs:
        pureCN_outputs = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].pureCN_outputs

        cancer_str = compose_lib_capture_str(cancer_capture)

        liqbio_cna = LiqbioCNAPlot()
        liqbio_cna.tumor_cnr = self.capture_to_results[cancer_capture].cnr
        liqbio_cna.tumor_cns = self.capture_to_results[cancer_capture].cns
        liqbio_cna.normal_cnr = self.capture_to_results[normal_capture].cnr
        liqbio_cna.normal_cns = self.capture_to_results[normal_capture].cns
        liqbio_cna.het_snps_vcf = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].vcf_addsample_output
        liqbio_cna.purecn_csv = pureCN_outputs["csv"]
        liqbio_cna.purecn_genes_csv = pureCN_outputs["genes_csv"]
        liqbio_cna.purecn_loh_csv = pureCN_outputs["loh_csv"]
        liqbio_cna.purecn_variants_csv = pureCN_outputs["variants_csv"]
        liqbio_cna.svcaller_T_DEL = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_T_DUP = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_T_INV = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_T_TRA = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_N_DEL = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_N_DUP = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_N_INV = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.svcaller_N_TRA = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.germline_mut_vcf = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.somatic_mut_vcf = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.plot_png = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.output_cna = # XXX CONTINUE HERE: FILL THESE IN
        liqbio_cna.output_purity = # XXX CONTINUE HERE: FILL THESE IN

        return liqbio_cna.output_cna, liqbio_cna.output_purity

    def configure_panel_analysis_cancer_vs_normal_liqbio(self, normal_capture, cancer_capture):
        cancer_capture_str = compose_lib_capture_str(cancer_capture)

        # XXX CONFIGURE SVCALLER, AND STORE THE OUTPUT FILES SOMEHOW -> 

        if self.refdata['targets'][cancer_capture_str]['purecn_targets']:
            self.configure_purecn(normal_capture, cancer_capture)

            self.configure_liqbio_cna(normal_capture, cancer_capture)