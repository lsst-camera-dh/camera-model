import os
import unittest
import camera_components

class CameraComponentsTestCase(unittest.TestCase):
    "Test case class for camera_components module."
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_create_from_etrav(self):
        """
        Expose bug deriving database name when using LIMS_URL for Prod
        eTraveler tables (LSSTTD-970).
        """
        os.environ['LCATR_UNIT_TYPE'] = 'LCA-11021_RTM'
        os.environ['LCATR_LIMS_URL'] = 'http://lsst-camera.slac.stanford.edu:80/eTraveler/Prod'

        raft_id = 'LCA-11021_RTM-003_ETU2'

        user = 'jchiang'
        db_name = None
        prodServer = True
        htype = 'LCA-11021_RTM'
        noBatched = 'false'

        raft = camera_components.Raft.create_from_etrav(raft_id, user=user,
                                                        db_name=db_name,
                                                        prodServer=prodServer,
                                                        htype=htype,
                                                        noBatched=noBatched)

        sensor = raft.sensor('S00')
        self.assertEqual(sensor.sensor_id, 'ITL-3800C-145')
        self.assertEqual(sensor.manufacturer_sn, '20429')

        sensor = raft.sensor('S21')
        self.assertEqual(sensor.sensor_id, 'ITL-3800C-146')
        self.assertEqual(sensor.manufacturer_sn, '20263')

        # Also test using Dev tables.
        os.environ['LCATR_LIMS_URL'] = 'http://lsst-camera.slac.stanford.edu:80/eTraveler/Dev'
        raft_id = 'LCA-11021_RTM-004_ETU2-Dev'
        raft = camera_components.Raft.create_from_etrav(raft_id, user=user,
                                                        db_name=db_name,
                                                        prodServer=prodServer,
                                                        htype=htype,
                                                        noBatched=noBatched)

        sensor = raft.sensor('S11')
        self.assertEqual(sensor.sensor_id, 'ITL-3800C-017-Dev')
        self.assertEqual(sensor.manufacturer_sn, '20459')

        sensor = raft.sensor('S22')
        self.assertEqual(sensor.sensor_id, 'ITL-3800C-103-Dev')
        self.assertEqual(sensor.manufacturer_sn, '20261')

if __name__ == '__main__':
    unittest.main()
