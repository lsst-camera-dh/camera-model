"""
Makes simulated data in the format that ccs will produce for multi-raft images
by copying input data from previously taken runs.
"""

from __future__ import print_function, absolute_import, division

import sys
import os
import glob
import yaml
import time

import astropy.io.fits as fits


ROOT_FOLDER = '/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive-test/LCA-11021_RTM'

SENSOR_SLOTS = ['S00', 'S01', 'S02',
                'S10', 'S11', 'S12',
                'S20', 'S21', 'S22']

def make_output_topdir_path(**kwargs):
    """ Build the output path for a particular set of data
    """
    return os.path.join(kwargs.get('root_folder_out', '.'))

def make_output_path(**kwargs):
    """ Build the output path for a particular set of data
    """
    kwargs['ordinal_str'] = "%06i" % kwargs['ordinal']
    return os.path.join(kwargs.get('root_folder_out', '.'),
                        kwargs['image_type'],
                        "MC_{date}_{ordinal_str}_R{raft_slot}S{sensor_slot}.fits".format(**kwargs))

def find_latest_acq_dir(**kwargs):
    """ Find the directory corresponding to the lastest run with particular set of criteria.

    This is done using glob, sorting the list of directories and return the last
    one on the list.

    Two criteria must be specified:
    raft_name: e.g., 'LCA-11021_RTM-010-Dev'
    acq_type_in: e.g., 'dark_raft_acq'

    Three criteria may be specified:
    root_folder_in: defaults to '/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive-test/LCA-11021_RTM'
    run_number_in: defaults to '*' (i.e., take all runs, and select the highest numbered one)
    jh_version: defaults to 'v0'
    activity_id_in: defaults to '*' (i.e., take all activity ids, and select the highest numbered one)

    """
    xpath = os.path.join(kwargs.get('root_folder_in', ROOT_FOLDER),
                         kwargs['raft_name'],
                         kwargs.get('run_number_in', '*'),
                         kwargs['acq_type_in'],
                         kwargs.get('jh_version', 'v0'),
                         kwargs.get('activity_id_in', '*'))
    dirlist = sorted(glob.glob(xpath))
    if len(dirlist) == 0:
        return None
    return dirlist[-1]


def get_image_type_from_path(filepath):
    """ Get the type of image file by parsing the input filepath
    """
    tokens = os.path.basename(filepath).split('_')[1:-2]
    s = ''
    first = True
    for token in tokens:
        if first:
            first = False
        else:
            s += '_'
        s += token
    return s


def get_slots_from_path(filepath):
    """ Get the raft and sensor slots by parsing the output filepath
    """
    token = os.path.basename(filepath).split('_')[-1]
    raft_slot = token[1:3]
    sensor_slot = token[4:6]
    return raft_slot, sensor_slot


def get_files_for_slot(acq_dir, slot_id):
    """ Find the input files for a particular slot in a raft
    """
    xpath = os.path.join(acq_dir, slot_id, '*.fits')
    d_out = {}
    files = glob.glob(xpath)
    for f in files:
        image_type = get_image_type_from_path(f)
        d_out[image_type] = f
    return d_out

def get_files_for_raft(acq_dir):
    """ Find the input files for all thes slot in a raft.
    """
    d_out = {}
    for slot in SENSOR_SLOTS:
        d_out[slot] = get_files_for_slot(acq_dir, slot)
    return d_out


def find_files_for_raft(**kwargs):
    """ Find the latest file for a particular acquistion type for a given raft.
    """
    acq_dir = find_latest_acq_dir(**kwargs)
    if acq_dir is None:
        return {}
    return get_files_for_raft(acq_dir)


def find_files_for_rafts(raft_map, **kwargs):
    """ Find the latest file for a particular acquistion type for all of the mapped rafts.
    """
    raft_args = kwargs.copy()
    d_out = {}
    for k, v in raft_map.items():
        raft_args.update(dict(raft_name=v['dev'],
                              run_number_in=v.get('run_number', '*')))
                              
        d_out[k] = find_files_for_raft(**raft_args)
    return d_out


def reorder_acq_dict(acq_dict):
    """ Re-order the dictionary of files associated with an acquisition to a format
    that is better suited for multi-raft outputs.
    """
    d_out = {}
    for raft, v in acq_dict.items():
        for sensor, vv in v.items():
            for image_type, fits_file in vv.items():
                if image_type not in d_out:
                    d_out[image_type] = {}
                if raft not in d_out[image_type]:
                    d_out[image_type][raft] = {}
                d_out[image_type][raft][sensor] = fits_file
    return d_out


def make_output_files_for_image(f_dict, image_type, **kwargs):
    """ Make a map from output file names to input
    files names for all the files in a single multi-raft image.
    """
    image_args = kwargs.copy()
    image_args['image_type'] = image_type
    fout_map = {}
    for raft, r_dict in f_dict.items():
        image_args['raft_slot'] = raft[1:]
        for sensor, finpath in r_dict.items():
            image_args['sensor_slot'] = sensor[1:]
            outfile = make_output_path(**image_args)
            fout_map[outfile] = finpath
    return fout_map


def make_output_files(image_dict, **kwargs):
    """ Make a map from output file names to input files
    names for all the files in a series of multi-raft images.
    """
    image_map = {}
    for image_type, file_dict in sorted(image_dict.items()):
        image_map[image_type] = make_output_files_for_image(file_dict, image_type, **kwargs)
        kwargs['ordinal'] += 1
    return image_map


def update_primary_header(raft_slot, sensor_slot, hdu):
    """
    Update the primary header for a FITS file will a sensor level image

    Parameters
    ----------
    raft_slot : str
        Name of the raft slot within the mulit-raft setup
    sensor_slot : str
        Name of the sensor slot within the raft
    hdu : fits.Image
        FITS image whose header is being updated
    """
    pass


def update_image_header(raft_slot, sensor_slot, hdu):
    """
    Update the image header for one of the readout segments.

    Parameters
    ----------
    raft_slot : str
        Name of the raft slot within the mulit-raft setup
    sensor_slot : str
        Name of the sensor slot within the raft
    hdu : fits.Image
        FITS image whose header is being updated
        """
    pass


def copy_sensor_image(input_sensor_file, output_sensor_file,
                      raft_slot, sensor_slot, **kwargs):
    """
    Copies a FITS file, with hooks to update the various image headers

    Parameters
    ----------
    input_sensor_file : str
        Name of the file to be copied
    output_sensor_file : str
        Destination
    raft_slot : str
        Name of the raft slot within the mulit-raft setup
    sensor_slot : str
        Name of the sensor slot within the raft

    Keyword arguments
    -----------------
    overwrite : bool, optional
        Flag indicating whether to overwrite an existing output file
    dry_run : bool, optional
        If true, just print output file names, but do not copy files
    """

    overwrite = kwargs.get('overwrite', True)
    dry_run = kwargs.get('dry_run', False)

    if dry_run:
        os.system("touch %s"% output_sensor_file)
        return
    hdulist = fits.open(input_sensor_file)

    update_primary_header(raft_slot, sensor_slot, hdulist[0])

    for ext_num in range(1, 16):
        update_image_header(raft_slot, sensor_slot, hdulist[ext_num])

    hdulist.writeto(output_sensor_file, clobber=overwrite)
    hdulist.close()


def copy_files(image_map, **kwargs):
    """ Copy the files from a mapped of single raft images to a series of multi-raft images.
    """
    for image_type, image_dict in sorted(image_map.items()):
        sys.stdout.write("%s: copying %i files.\n" % (image_type, len(image_dict)))
        first = True
        for f_to, f_from in sorted(image_dict.items()):
            raft_slot, sensor_slot = get_slots_from_path(f_to)
            if first:
                first = False
                d_to = os.path.dirname(f_to)
                if not os.path.exists(d_to):
                    os.makedirs(d_to)
            copy_sensor_image(f_from, f_to, raft_slot, sensor_slot, **kwargs)


def read_input_camera_yaml(yamlpath):
    """ Read the raft slot to sensor ID mapping from a yamlfile
    """
    with open(yamlpath, 'r') as yaml_input:
        yaml_dict = yaml.load(yaml_input)
    return yaml_dict


class FakeCamera(object):
    '''
    Class that can copy a raft's worth of images.

    Parameters
    ----------
    raft_map : dict
        Dictionary mapping raft slots to raft ids.

    '''

    def __init__(self, raft_map, **kwargs):
        """
        Class constructor.
        """
        self.raft_map = raft_map
        self.datestring = kwargs.get('date', "%4i%02i%02i" % (time.localtime()[0:3]))
        self.ordinal = kwargs.get('ordinal', 0)
        

    @classmethod
    def create_from_yaml(cls, yamlpath):
        """ Create a CameraImages object from an input yaml file
        """
        raft_map = read_input_camera_yaml(yamlpath)
        return cls(raft_map)


    def build_acq_dict(self, acq_type, **kwargs):
        """ Create a dictionary of all the files for a set of acquistions
        corresopnding to the rafts in the raft map
        """
        arg_dict = kwargs.copy()
        arg_dict['acq_type'] = acq_type
        in_acq_dict = find_files_for_rafts(self.raft_map, **arg_dict)
        out_acq_dict = reorder_acq_dict(in_acq_dict)
        return out_acq_dict


    def build_output_filemap(self, **kwargs):
        """ Create a dictionary mapping output files to input files
        for all of the images for all of the rafts in the raft map
        """
        acq_type = kwargs['acq_type_in']
        kwargs.setdefault('ordinal', self.ordinal)
        kwargs.setdefault('date', self.datestring)
        acq_dict = self.build_acq_dict(acq_type, **kwargs)
        filemap = make_output_files(acq_dict, **kwargs)
        return filemap


    def run(self, **kwargs):
        """ Chain everyting together
        """
        filemap = self.build_output_filemap(**kwargs)
        copy_files(filemap, **kwargs)


def test():
    """ Hook for testing """
    import argparse

    usage = "usage: %(prog)s [options]"
    description = "Copy raft-level test data into the multi-raft format"

    parser = argparse.ArgumentParser(usage, description=description)
    parser.add_argument('-c', "--camera",
                        type=argparse.FileType('r'),
                        help="Yaml file with the raft slot mapping")
    parser.add_argument('-a', "--acq",
                        type=argparse.FileType('r'),
                        help="Yaml file with the input and output run information")
    args = parser.parse_args(sys.argv[1:])

    fake_cam = FakeCamera.create_from_yaml(args.camera.name)
    fake_cam.run(args.acq.name, overwrite=True)


if __name__ == '__main__':
    test()
