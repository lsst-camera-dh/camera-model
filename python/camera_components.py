"""
Abstractions of rafts and sensors.
"""
from __future__ import print_function, absolute_import, division
import os
import copy
import pickle
import yaml
import eTraveler.clientAPI.connection
from lsst.obs.lsst import LsstCamMapper

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
    def __init__(self, user='ccs', db_name='Prod',
                 hardware_info_file='lca-10134_cryostat-0001_hw_info.pkl'):
        # Get the detector names from the obs_lsst package.
        self.camera_object = LsstCamMapper().camera
        self._det_names = [det.getName() for det in self.camera_object]
        self.raft_names = set()
        self.slot_names = set()
        for det_name in self._det_names:
            raft, slot = det_name.split('_')
            self.raft_names.add(raft)
            self.slot_names.add(slot)
        self.raft_names = sorted(list(self.raft_names))
        self.slot_names = sorted(list(self.slot_names))

        self._science_rafts = None
        self._corner_rafts = None
        self._conn_args = user, db_name
        self._hardware_info_file = hardware_info_file

    def _get_installed_locations(self):
        if os.path.isfile(self._hardware_info_file):
            with open(self._hardware_info_file, 'rb') as fd:
                resp = pickle.load(fd)
        else:
            conn = eTraveler.clientAPI.connection.Connection(*self._conn_args)
            cryostat_id = 'LCA-10134_Cryostat-0001'
            htype = 'LCA-10134_Cryostat'
            resp = conn.getHardwareHierarchy(experimentSN=cryostat_id,
                                             htype=htype)
            with open(self._hardware_info_file, 'wb') as fd:
                pickle.dump(resp, fd)

        self._science_rafts = []
        self._corner_rafts = []
        for item in resp:
            if item['child_hardwareTypeName'] == 'LCA-11021_RTM':
                self._science_rafts.append('R' + item['slotName'][-2:])
            if item['child_hardwareTypeName'] == 'LCA-10692_CRTM':
                self._corner_rafts.append('R' + item['slotName'][-2:])

    def get_det_names(self):
        """Return a copy of the list of detector names."""
        return copy.copy(self._det_names)

    def get_raft_names(self):
        """Return a copy of the list of raft names."""
        return copy.copy(self.raft_names)

    @property
    def installed_science_rafts(self):
        """A list of installed science rafts."""
        if self._science_rafts is None:
            self._get_installed_locations()
        return self._science_rafts

    @property
    def installed_corner_rafts(self):
        """A list of installed corner rafts."""
        if self._corner_rafts is None:
            self._get_installed_locations()
        return self._corner_rafts

    def get_installed_raft_names(self):
        """Return the list of installed raft names."""
        return self.installed_science_rafts + self.installed_corner_rafts

    def get_installed_det_names(self):
        """Return the detector names based on the installed rafts."""
        det_names = self.get_det_names()
        raft_names = self.get_installed_raft_names()
        return [_ for _ in det_names if _[:3] in raft_names]

    def get_installed_science_raft_det_names(self):
        """Return the detector names based on the installed rafts."""
        det_names = self.get_det_names()
        return [_ for _ in det_names if _[:3] in self.installed_science_rafts]

    def get_installed_corner_raft_det_names(self):
        """Return the detector names based on the installed rafts."""
        det_names = self.get_det_names()
        return [_ for _ in det_names if _[:3] in self.installed_corner_rafts]

    def get_slot_names(self):
        """Return a copy of the list of slot names."""
        return copy.copy(self.slot_names)

# Create a single concrete instance of the CameraInfo class that
# clients can import directly.
camera_info = CameraInfo()
