import unittest
from autoseq.tools.intervals import *


class TestIntervals(unittest.TestCase):
    def test_slop_interval_list(self):
        slop_interval_list = SlopIntervalList()
        slop_interval_list.input = "test_input"
        slop_interval_list.output = "test_output"
        cmd = slop_interval_list.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_interval_list_to_bed(self):
        interval_list_to_bed = IntervalListToBed()
        interval_list_to_bed.input = "test_input"
        interval_list_to_bed.output = "test_output"
        cmd = interval_list_to_bed.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_msi_sensor_scan(self):
        msi_sensor_scan = MsiSensorScan()
        msi_sensor_scan.input_fasta = "input.fa"
        msi_sensor_scan.output = "test_output"
        cmd = msi_sensor_scan.command()
        self.assertIn('input.fa', cmd)
        self.assertIn('test_output', cmd)

    def test_intersect_msi_sites(self):
        msi_sensor_scan = IntersectMsiSites()
        msi_sensor_scan.input_msi_sites = "test_input"
        msi_sensor_scan.target_bed = "test.bed"
        msi_sensor_scan.output_msi_sites = "test_output"
        cmd = msi_sensor_scan.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test.bed', cmd)
        self.assertIn('test_output', cmd)

    def test_msi_sensor(self):
        msi_sensor = MsiSensor()
        msi_sensor.msi_sites = "dummy_msi_sites"
        msi_sensor.input_normal_bam = "normal.bam"
        msi_sensor.input_tumor_bam = "tumor.bam"
        msi_sensor.output = "test_output"
        cmd = msi_sensor.command()
        self.assertIn('dummy_msi_sites', cmd)
        self.assertIn('normal.bam', cmd)
        self.assertIn('tumor.bam', cmd)
        self.assertIn('test_output', cmd)
