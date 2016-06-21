from pypedream.job import Job, required


class CompileMetadata(Job):
    def __init__(self, referral_db_conf, blood_barcode, tumor_barcode, output_json, addresses):
        Job.__init__(self)
        self.referral_db_conf = referral_db_conf
        self.blood_barcode = blood_barcode
        self.tumor_barcode = tumor_barcode
        self.output_json = output_json
        self.addresses = addresses

    def command(self):
        # compileMetadata 3098121 3098849 --db_config $HOME/repos/reportgen/tests/referral-db-config.json \
        #  --output /dev/stdout  --address_table_file reportgen/assets/addresses.csv
        return "compileMetadata" + \
               required('', self.blood_barcode) + \
               required('', self.tumor_barcode) + \
               required('--db_config ', self.referral_db_conf) + \
               required('--address_table_file ', self.addresses) + \
               required('--output ', self.output_json)


class CompileAlasccaGenomicJson(Job):
    def __init__(self, input_somatic_vcf, input_cn_calls, input_msisensor, output_json):
        Job.__init__(self)
        self.input_somatic_vcf = input_somatic_vcf
        self.input_cn_calls = input_cn_calls
        self.input_msisensor = input_msisensor
        self.output_json = output_json

    def command(self):
        return 'compileAlasccaGenomicReport ' + \
            required('', self.input_somatic_vcf) + \
            required('', self.input_cn_calls) + \
            required('', self.input_msisensor) + \
            required('--output ', self.output_json)

