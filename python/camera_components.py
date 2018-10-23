"""
Abstractions of rafts and sensors.
"""
from __future__ import print_function, absolute_import, division
import os
import copy
import yaml
import eTraveler.clientAPI.connection

__all__ = ['Raft', 'Sensor', 'REB', 'ROOT_FOLDER', 'camera_info']

ROOT_FOLDER = os.environ.get('LCATR_DATACATALOG_FOLDER',
                             'LSST/mirror/SLAC-prod/prod')
USER = os.environ['USER']

def parse_etraveler_response(rsp, validate):
    """ Convert the response from an eTraveler clientAPI query to a
    key,value pair

    Parameters
    ----------
    rsp : return type from
        eTraveler.clientAPI.connection.Connection.getHardwareHierarchy
        which is an array of dicts information about the 'children' of a
        particular hardware element.
    validate : dict
        A validation dictionary, which contains the expected values
        for some parts of the rsp.  This is here for sanity checking,
        for example requiring that the parent element matches the
        input element to the request.

    Returns
    ----------
    slot_name,child_esn:
    slot_name  : str
        A string given to the particular 'slot' for each child
    child_esn : str
        The sensor id of the child, e.g., E2V-CCD250-104
    """
    for key, val in validate.items():
        try:
            rsp_val = rsp[key]
            if isinstance(val, list):
                if rsp_val not in val:
                    errmsg = "eTraveler response does not match expectation for key %s: " % (key)
                    errmsg += "%s not in %s" % (rsp_val, val)
                    raise ValueError(errmsg)
            else:
                if rsp_val != val:
                    errmsg = "eTraveler response does not match expectation for key %s: " % (key)
                    errmsg += "%s != %s" % (rsp_val, val)
                    raise ValueError(errmsg)
        except KeyError:
            raise KeyError("eTraveler response does not include expected key %s" % (key))

    child_esn = rsp['child_experimentSN']
    slot_name = rsp['slotName']
    return slot_name, child_esn


class Sensor:
    '''
    A simple class to carry around some information about sensors in a raft.

    Parameters
    ----------
    sensor_id : str
        Name of the sensor, e.g., 'E2V-CCD250-104'
    raft_id : str
        Name of the associated raft
    manufacturer_sn : str
        Manufacturer's serial number.
    '''
    def __init__(self, sensor_id, raft_id, manufacturer_sn):
        """
        Class constructor.
        """
        self.__sensor_id = str(sensor_id)
        self.__raft_id = str(raft_id)
        self._manufacturer_sn = manufacturer_sn

    @property
    def sensor_id(self):
        """ Return the name of the sensor, e.g., 'E2V-CCD250-104' """
        return self.__sensor_id

    @property
    def raft_id(self):
        """ Return the name of the raft, e.g., 'RAFT-000' """
        return self.__raft_id

    @property
    def manufacturer_sn(self):
        "The manufacturer's serial number."
        return self._manufacturer_sn


class REB:
    '''
    Class to contain information on a REB (raft electronics board)
    that's been extracted from the eTraveler database tables.
    '''
    def __init__(self, reb_id, manufacturer_sn, firmware_version):
        '''
        Parameters
        ----------
        reb_id : str
            LSST ID number of the REB.
        manufacturer_sn : str
            Manufacturer's serial number.  This should be the hex
            representation of the internal hardware id.
        firmware_version : str
            The firmware version that is supposed to be installed.
        '''
        self._reb_id = reb_id
        self._manufacturer_sn = manufacturer_sn
        self._firmware_version = firmware_version

    @property
    def reb_id(self):
        return self._reb_id

    @property
    def manufacturer_sn(self):
        return self._manufacturer_sn

    @manufacturer_sn.setter
    def manufacturer_sn(self, value):
        self._manufacturer_sn = value

    @property
    def firmware_version(self):
        return self._firmware_version

    @staticmethod
    def get_rebs(raft_id, conn=None, htype='LCA-11021_RTM', db_name='Prod'):
        """
        Factory method for creating a dictionary of REB objects given
        a raft_id.

        Parameters
        ----------
        raft_id : str
            The LSST ID of the raft.
        conn : eTraveler.clientAPI.connection.Connection, optional
            If None, then a connection object will be created.
        htype : str, optional
            LSST hardware type. Default: 'LCA-11021_RTM'
        db_name : str, optional
            eTraveler database to use.  Default: 'Prod'

        Returns
        -------
        dict of REB objects keyed by slot name.
        """
        if conn is None:
            conn = eTraveler.clientAPI.connection.Connection(USER, db_name)
        resp = conn.getHardwareHierarchy(experimentSN=raft_id, htype=htype)
        rebs = dict()
        for item in resp:
            if item['slotName'].startswith('REB'):
                slot = item['slotName']
                reb_id = item['child_experimentSN']
                htype = item['child_hardwareTypeName']
                manufacturer_sn = conn.getManufacturerId(experimentSN=reb_id,
                                                         htype=htype)
                firmware_version = None   # There is no interface to this yet.
                rebs[slot] = REB(reb_id, manufacturer_sn, firmware_version)
        return rebs

class Raft:
    '''
    A simple class to carry around some information about a raft.

    Parameters
    ----------
    raft_id : str
        Name of the raft
    sensor_type : str
        Type of sensors in the raft, either 'e2v-CCD' or 'ITL-CCD'
    sensor_dict : dict
        Dictionary for slot to Sensor
    rebs : dict
        dict of REB objects keyed by slot name.
    '''
    def __init__(self, raft_id, sensor_type, sensor_dict, rebs=None):
        """
        Class constructor.
        """
        self.__raft_id = raft_id
        self.__sensor_type = sensor_type
        self.__sensor_dict = sensor_dict
        self._rebs = dict()
        if rebs is not None:
            self._rebs = rebs

    @staticmethod
    def create_from_yaml(yamlfile):
        """ Create a Raft object from a yaml file """
        input_dict = yaml.safe_load(open(yamlfile))
        raft_id = input_dict['raft_id']
        sensor_type = input_dict['sensor_type']
        sensors = input_dict['sensors']
        sensor_dict = {}
        for slot_name, sensor_name in sensors.items():
            sensor_dict[slot_name] = Sensor(sensor_name, raft_id, None)
        return Raft(raft_id, sensor_type, sensor_dict)

    @staticmethod
    def create_from_etrav(raft_id, **kwargs):
        """ Create a Raft object from query to the eTraveler

        Parameters
        ----------
        raft_id : str
            Name of the raft, this must match the 'parent_experimentSN' field
            in the eTraveler db.

        Keyword Arguments
        ----------
        user   : str
            Expected by the eTraveler interface
        db_name : str [None]
            Version of the eTraveler to query.
            If None, then derive db_name from the LCATR_LIMS_URL
            environment variable.
        prodServer : bool [True]
        htype : str ['LCA-10753-RSA_sim']
            Hardware type, this must match the 'parent_hardware_type' field
            in the eTraveler db.
        noBatched : str ['false']

        Returns
        ----------
        Newly created Raft object
        """
        user = kwargs.get('user', USER)
        db_name = kwargs.get('db_name', None)
        prod_server = kwargs.get('prod_server', True)
        htype = kwargs.get('htype', 'LCA-11021_RTM')
        no_batched = kwargs.get('no_batched', 'false')
        if db_name is None:
            db_name = os.path.split(os.environ['LCATR_LIMS_URL'])[-1]
            if db_name not in 'Prod Dev Test Raw'.split():
                # This case occurs when using the fake_eT server.
                db_name = os.environ.get('LCATR_ET_DB_NAME', 'Dev')

        my_conn = eTraveler.clientAPI.connection.Connection(user, db_name,
                                                            prod_server)
        return Raft.create_from_connection(my_conn, raft_id, htype,
                                           no_batched=no_batched)

    @staticmethod
    def create_from_connection(connection, raft_id, htype,
                               no_batched='false'):
        """ Create a Raft object from query to the eTraveler

        Parameters
        ----------
        connection : 'eTraveler/clientAPI/connection.Connection'
            Object that wraps connection to eTraveler database
        raft_id : str
            Name of the raft, this must match the 'parent_experimentSN' field
            in the eTraveler db.
        htype : str
            Hardware type, this must match the 'parent_hardwareTypeName' field
            in the eTraveler db.
        no_batched : str ['false']

        Returns
        ----------
        Newly created Raft
        """
        rsp = connection.getHardwareHierarchy(experimentSN=raft_id,
                                              htype=htype,
                                              noBatched=no_batched)
        sensor_dict = {}

        ccd_types = ['e2v-CCD', 'ITL-CCD']

        validate_dict = dict(child_hardwareTypeName=ccd_types)

        sensor_type = None

        for rsp_item in rsp:
            if rsp_item['child_hardwareTypeName'] in ccd_types:
                if sensor_type is None:
                    sensor_type = rsp_item['child_hardwareTypeName']
                slot, c_esn = parse_etraveler_response(rsp_item, validate_dict)
                manu_sn = connection.getManufacturerId(experimentSN=c_esn,
                                                       htype=sensor_type)
                sensor_dict[str(slot)] = Sensor(c_esn, raft_id, manu_sn)

        rebs = REB.get_rebs(raft_id, conn=connection, htype=htype)
        return Raft(raft_id, sensor_type, sensor_dict, rebs)

    @property
    def raft_id(self):
        """ The name of this raft """
        return self.__raft_id

    @property
    def sensor_type(self):
        """ The type of sensors in this raft.  'e2v-CCD' or 'ITL-CCD' """
        return self.__sensor_type

    @property
    def slot_names(self):
        """ The names of the 'slots' associated with the sensors """
        slots = list(self.__sensor_dict.keys())
        slots.sort()
        return slots

    @property
    def sensor_names(self):
        """ The names of the sensors in this raft, sorted to match the
        slot names """
        return [self.__sensor_dict[slot].sensor_id for slot in self.slot_names]

    def items(self):
        """ Iterator over slot_name, sensor_name pairs """
        return zip(self.slot_names, self.sensor_names)

    def sensor(self, slot):
        """ Sensor associated with a particular slot """
        return self.__sensor_dict[slot]

    @property
    def rebs(self):
        return self._rebs

class CameraInfo:
    """
    This class will eventually be an interface to the camera info
    available from the obs_lsst package.  For now, use hard-coded
    values for things like the generic detector and raft names.
    """
    _det_names = ('R01_S00', 'R01_S01', 'R01_S02', 'R01_S10', 'R01_S11', 'R01_S12', 'R01_S20', 'R01_S21', 'R01_S22', 'R02_S00', 'R02_S01', 'R02_S02', 'R02_S10', 'R02_S11', 'R02_S12', 'R02_S20', 'R02_S21', 'R02_S22', 'R03_S00', 'R03_S01', 'R03_S02', 'R03_S10', 'R03_S11', 'R03_S12', 'R03_S20', 'R03_S21', 'R03_S22', 'R10_S00', 'R10_S01', 'R10_S02', 'R10_S10', 'R10_S11', 'R10_S12', 'R10_S20', 'R10_S21', 'R10_S22', 'R11_S00', 'R11_S01', 'R11_S02', 'R11_S10', 'R11_S11', 'R11_S12', 'R11_S20', 'R11_S21', 'R11_S22', 'R12_S00', 'R12_S01', 'R12_S02', 'R12_S10', 'R12_S11', 'R12_S12', 'R12_S20', 'R12_S21', 'R12_S22', 'R13_S00', 'R13_S01', 'R13_S02', 'R13_S10', 'R13_S11', 'R13_S12', 'R13_S20', 'R13_S21', 'R13_S22', 'R14_S00', 'R14_S01', 'R14_S02', 'R14_S10', 'R14_S11', 'R14_S12', 'R14_S20', 'R14_S21', 'R14_S22', 'R20_S00', 'R20_S01', 'R20_S02', 'R20_S10', 'R20_S11', 'R20_S12', 'R20_S20', 'R20_S21', 'R20_S22', 'R21_S00', 'R21_S01', 'R21_S02', 'R21_S10', 'R21_S11', 'R21_S12', 'R21_S20', 'R21_S21', 'R21_S22', 'R22_S00', 'R22_S01', 'R22_S02', 'R22_S10', 'R22_S11', 'R22_S12', 'R22_S20', 'R22_S21', 'R22_S22', 'R23_S00', 'R23_S01', 'R23_S02', 'R23_S10', 'R23_S11', 'R23_S12', 'R23_S20', 'R23_S21', 'R23_S22', 'R24_S00', 'R24_S01', 'R24_S02', 'R24_S10', 'R24_S11', 'R24_S12', 'R24_S20', 'R24_S21', 'R24_S22', 'R30_S00', 'R30_S01', 'R30_S02', 'R30_S10', 'R30_S11', 'R30_S12', 'R30_S20', 'R30_S21', 'R30_S22', 'R31_S00', 'R31_S01', 'R31_S02', 'R31_S10', 'R31_S11', 'R31_S12', 'R31_S20', 'R31_S21', 'R31_S22', 'R32_S00', 'R32_S01', 'R32_S02', 'R32_S10', 'R32_S11', 'R32_S12', 'R32_S20', 'R32_S21', 'R32_S22', 'R33_S00', 'R33_S01', 'R33_S02', 'R33_S10', 'R33_S11', 'R33_S12', 'R33_S20', 'R33_S21', 'R33_S22', 'R34_S00', 'R34_S01', 'R34_S02', 'R34_S10', 'R34_S11', 'R34_S12', 'R34_S20', 'R34_S21', 'R34_S22', 'R41_S00', 'R41_S01', 'R41_S02', 'R41_S10', 'R41_S11', 'R41_S12', 'R41_S20', 'R41_S21', 'R41_S22', 'R42_S00', 'R42_S01', 'R42_S02', 'R42_S10', 'R42_S11', 'R42_S12', 'R42_S20', 'R42_S21', 'R42_S22', 'R43_S00', 'R43_S01', 'R43_S02', 'R43_S10', 'R43_S11', 'R43_S12', 'R43_S20', 'R43_S21', 'R43_S22')
    def __init__(self):
#        # Get the detector names from the obs_lsst package.
#        self._det_names = [det.getName() for det
#                           in obs_lsst.LsstCamMapper().camera
#                           if det.getType()==cameraGeom.SCIENCE]
        self.raft_names = set()
        self.slot_names = set()
        for det_name in self._det_names:
            raft, slot = det_name.split('_')
            self.raft_names.add(raft)
            self.slot_names.add(slot)
        self.raft_names = sorted(list(self.raft_names))
        self.slot_names = sorted(list(self.slot_names))
        # Remove "_" from detector names to match CCS naming convention
        # for focalplane-level FITS files.
        self.det_names = [x.replace('_', '') for x in self._det_names]

    def get_det_names(self):
        """Return a copy of the list of detector names."""
        return copy.copy(self.det_names)

    def get_raft_names(self):
        """Return a copy of the list of raft names."""
        return copy.copy(self.raft_names)

    def get_slot_names(self):
        """Return a copy of the list of slot names."""
        return copy.copy(self.slot_names)

# Create a single concrete instance of the CameraInfo class that
# clients can import directly.
camera_info = CameraInfo()
