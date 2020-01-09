#!/usr/bin/env python3

import os
import shutil
import subprocess
import time
import unittest


class MSTSIntegration(unittest.TestCase):

    def setUp(self):
        self.project_path = "../../bin"
        self.data_path = "../../test-data"

        # create dir/env
        rundir = 'integration'
        os.mkdir(rundir)
        os.chdir(rundir)

    def tearDown(self):

        os.chdir('..')
#        shutil.rmtree('integration')

    def step1__MSTS_converter(self):
        '''run MSTS_converter.py'''

        f = 'MSTS_converter.py'
        args = ['{}/mapping.bam'.format(self.data_path),"--bed", "--size", "--wig", "-m", "fragment-middle", "-w", "20","-p", "TEST" ]
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("TEST.bed"))
        self.assertTrue(os.path.isfile("TEST.size"))
        self.assertTrue(os.path.isfile("TEST.wig"))

    def step2__wigToBigWig(self):
        '''run wigToBigWig'''

        args = ['TEST.wig', '{}/mapping.genome'.format(self.data_path), 'TEST.bw']
        cmd = ['wigToBigWig']
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("TEST.bw"))

    def step3__MSTS_phasogram(self):
        '''run MSTS_phasogram.py'''

        f = 'MSTS_phasogram.py'
        args = ['TEST.bw', '-o', 'TEST_PHASO.png', '--regression']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("TEST_PHASO.png"))

    def step4__MSTS_feature_phasogram(self):
        '''run MSTS_feature_phasogram'''

        f = 'MSTS_feature_phasogram.py'
        args = ['TEST.bw', '{}/annotations.gff3'.format(self.data_path), '-o', 'TEST_FEAT_PHASO.png', '--context','--norm', '--GaussianSmoothing','-v','2']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("TEST_FEAT_PHASO.png"))

    def step5__MSTS_count_TPM(self):
        '''run MSTS_count_TPM'''

        stdout = open('counts.TPM', 'w')
        f = 'MSTS_count_TPM.py'
        args = ['{}/mapping.RNA.bam'.format(self.data_path), '{}/annotations.gff3'.format(self.data_path), '-s', 'reverse' ]
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        self.assertTrue(os.path.isfile("counts.TPM"))

    def step6__MSTS_merge_phasogram(self):
        '''run merge phasogram'''

        stdout = open('50.tpm', 'w')
        p1 = subprocess.Popen(["tail","-n+2","counts.TPM"], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", 'BEGIN{FS="\t"}{if($5 > 50){print $1;}}'], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        stdout.write(output.decode())
        stdout.close()
        self.assertTrue(os.path.isfile("50.tpm"))

        stdout = open('50-5.tpm', 'w')
        p1 = subprocess.Popen(["tail","-n+2","counts.TPM"], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", 'BEGIN{FS="\t"}{if($5 <= 50 && $5 > 5){print $1}}'], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        stdout.write(output.decode())
        stdout.close()
        self.assertTrue(os.path.isfile("50-5.tpm"))

        stdout = open('5-1.tpm', 'w')
        p1 = subprocess.Popen(["tail","-n+2","counts.TPM"], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", 'BEGIN{FS="\t"}{if($5 <= 5 && $5 > 1){print $1}}'], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        stdout.write(output.decode())
        stdout.close()
        self.assertTrue(os.path.isfile("5-1.tpm"))

        stdout = open('1.tpm', 'w')
        p1 = subprocess.Popen(["tail","-n+2","counts.TPM"], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", 'BEGIN{FS="\t"}{if($5 < 1){print $1}}'], stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        output = p2.communicate()[0]
        stdout.write(output.decode())
        stdout.close()
        self.assertTrue(os.path.isfile("1.tpm"))

        f = 'MSTS_feature_phasogram.py'
        stdout = open('50.tpm.phaso', 'w')
        args = ['TEST.bw', '{}/annotations.gff3'.format(self.data_path), '-o', 'PHASO_50.png', '--context', '--GaussianSmoothing', '-ft', 'mRNA', '--flush','-l', '50.tpm']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        stdout.close()
        self.assertTrue(os.path.isfile("50.tpm.phaso"))

        stdout = open('50-5.tpm.phaso', 'w')
        args = ['TEST.bw', '{}/annotations.gff3'.format(self.data_path), '-o', 'PHASO_50-5.png', '--context', '--GaussianSmoothing', '-ft', 'mRNA', '--flush','-l', '50-5.tpm']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        stdout.close()
        self.assertTrue(os.path.isfile("50-5.tpm.phaso"))

        stdout = open('5-1.tpm.phaso', 'w')
        args = ['TEST.bw', '{}/annotations.gff3'.format(self.data_path), '-o', 'PHASO_5-1.png', '--context', '--GaussianSmoothing', '-ft', 'mRNA', '--flush','-l', '5-1.tpm']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        stdout.close()
        self.assertTrue(os.path.isfile("5-1.tpm.phaso"))

        stdout = open('1.tpm.phaso', 'w')
        args = ['TEST.bw', '{}/annotations.gff3'.format(self.data_path), '-o', 'PHASO_1.png', '--context', '--GaussianSmoothing', '-ft', 'mRNA', '--flush','-l', '1.tpm']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        stdout.close()
        self.assertTrue(os.path.isfile("1.tpm.phaso"))

        phaso_list = open('list.txt', 'w')
        phaso_list.write("1.tpm.phaso\t1TPM>X\n")
        phaso_list.write("5-1.tpm.phaso\t5TPM>X>1TPM\n")
        phaso_list.write("50-5.tpm.phaso\t50TPM>X>5TPM\n")
        phaso_list.write("50.tpm.phaso\tX>50TPM\n")
        phaso_list.close()
        self.assertTrue(os.path.isfile("list.txt"))

        f = 'MSTS_merge_phasograms.py'
        args = ['list.txt']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("phaso.merge.png"))

    def step7__MSTS_detect_nucleosomes(self):
        '''run MSTS_detect_nucleosomes'''

        f = 'MSTS_detect_nucleosomes.py'
        args = ['TEST.bw', '-p', 'DETECT', '--bed', '--wig', '-v', '2', '--refine' ]
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("DETECT.nucleosomes.txt"))

    def step8__bedToBigBed(self):
        '''run wigToBigWig'''

        stdout = open('TEST.sort.bed', 'w')
        args = ['-k1,1', '-k2,2n', 'TEST.bed']
        cmd = ['sort']
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        stdout.close()
        self.assertTrue(os.path.isfile("TEST.sort.bed"))

        args = ['TEST.sort.bed', '{}/mapping.genome'.format(self.data_path), 'TEST.sort.bb']
        cmd = ['bedToBigBed']
        cmd.extend(args)
        subprocess.call(cmd)
        self.assertTrue(os.path.isfile("TEST.sort.bb"))


    def step9__MSTS_dinuc_frequency(self):

        stdout = open('dinuc.freq', 'w')
        f = 'MSTS_dinuc_frequency.py'
        args = ['{}/genome.fasta'.format(self.data_path),'TEST.sort.bb', '--flush', '--pFreqNormMix']
        cmd = ['python', '{}/{}'.format(self.project_path,f)]
        cmd.extend(args)
        subprocess.call(cmd, stdout=stdout)
        stdout.close()
        self.assertTrue(os.path.isfile("dinuc.freq"))
        self.assertTrue(os.path.isfile("freq_ATGC_Normalized.png"))



    def steps(self):
        '''
        Generates the step methods from their parent object
        '''
        for name in sorted(dir(self)):
            if name.startswith('step'):
                yield name, getattr(self, name)

    def test_steps(self):
        '''
        Run the individual steps associated with this test
        '''
        # Create a flag that determines whether to raise an error at
        # the end of the test
        failed = False

        # An empty string that the will accumulate error messages for 
        # each failing step
        fail_message = ''
        for name, step in self.steps():
            try:
                step()
            except Exception as e:
              # A step has failed, the test should continue through
              # the remaining steps, but eventually fail
                failed = True

              # get the name of the method -- so the fail message is
              # nicer to read :)
                name = name.split('__')[1]
              # append this step's exception to the fail message
                fail_message += "\n\nFAIL: {}\n {} failed ({}: {})".format(name, step, type(e), e)

      # check if any of the steps failed
            if failed is True:
            # fail the test with the accumulated exception message
                self.fail(fail_message)

if __name__ == "__main__":

    suite = unittest.TestLoader().loadTestsFromTestCase(MSTSIntegration)
    unittest.TextTestRunner(verbosity=2).run(suite)
