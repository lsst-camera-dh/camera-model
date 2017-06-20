"Unit tests for simulation.fake_raft module"
import os
import unittest
import simulation.fake_raft as fake_raft

class FakeRaftTestCase(unittest.TestCase):
    "Test case class for fake_raft module."
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_make_outfile_path(self):
        "Unit test for make_outfile_path"
        slot_name = 'S00'
        file_string = 'ITL-3800C-023-Dev_fe55_fe55_000_4549D_20170128205319.fits'
        jobname = 'fe55_raft_acq_sim'
        outfile_path = fake_raft.make_outfile_path(outpath='.',
                                                   slot_name=slot_name,
                                                   file_string=file_string,
                                                   jobname=jobname)
        expected_fn = ('ITL-3800C-023-Dev_fe55_fe55_000_4549D_'
                       + jobname + '_20170128205319.fits')
        self.assertEqual(os.path.join('.', slot_name, expected_fn),
                         outfile_path)

if __name__ == '__main__':
    unittest.main()
