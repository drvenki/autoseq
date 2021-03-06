from pypedream.job import Job


class Svcaller(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bam = None
        self.event_type = None
        self.output_bam = None
        self.output_gtf = None
        self.reference_sequence = None
        self.scratch = None
        self.jobname = "svcaller-run-all"

    def command(self):
        activate_env_cmd = "source activate svcallerenv"

        run_all_cmd = ("svcaller run-all --tmp-dir {scratch} " +
                      "--event-type {event_type} " +
                      "--fasta-filename {reference_seq} " +
                      "--filter-event-overlap --events-gtf {output_gtf} "
                      "--events-bam {output_bam} {input_bam}").format(
                          scratch=self.scratch,
                          event_type=self.event_type,
                          reference_seq=self.reference_sequence,
                          output_gtf=self.output_gtf,
                          output_bam=self.output_bam,
                          input_bam=self.input_bam,
                      )

        deactivate_env_cmd = "source deactivate"

        return "{} && {} && {}".format(
            activate_env_cmd,
            run_all_cmd,
            deactivate_env_cmd,
        )


class Sveffect(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_del_gtf = None
        self.input_dup_gtf = None
        self.input_inv_gtf = None
        self.input_tra_gtf = None
        self.ts_regions = None
        self.ar_regions = None
        self.fusion_regions = None
        self.output_combined_bed = None
        self.output_effects_json = None
        self.jobname = "sveffect"

    def command(self):
        activate_env_cmd = "source activate svcallerenv"

        make_bed_cmd = ("sveffect make-bed " +
                       "--del-gtf {del_gtf} " +
                       "--dup-gtf {dup_gtf} " +
                       "--inv-gtf {inv_gtf} " +
                       "--tra-gtf {tra_gtf} " +
                       "{output_combined_bed}").format(
                           del_gtf=self.input_del_gtf,
                           dup_gtf=self.input_dup_gtf,
                           inv_gtf=self.input_inv_gtf,
                           tra_gtf=self.input_tra_gtf,
                           output_combined_bed=self.output_combined_bed
                       )

        predict_cmd = ("sveffect predict " +
                      "--ts-regions {ts_regions} " +
                      "--ar-regions {ar_regions} " +
                      "--fusion-regions {fusion_regions} " +
                      "--effects-filename {output_effects_json} {combined_effects_bed}").format(
                          ts_regions=self.ts_regions,
                          ar_regions=self.ar_regions,
                          fusion_regions=self.fusion_regions,
                          output_effects_json=self.output_effects_json,
                          combined_effects_bed=self.output_combined_bed,
                      )

        deactivate_env_cmd = "source deactivate"

        return "{} && {} && {} && {}".format(
            activate_env_cmd,
            make_bed_cmd,
            predict_cmd,
            deactivate_env_cmd,
        )