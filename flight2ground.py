#!/usr/bin/env python

# Required packages
import argparse
from configparser import ConfigParser

import requests
from astropy.time import Time
from astropy.io import fits
import numpy as np


def DMS_to_detector(data, detector):
    """Transformations from Robert Jedzrejewski
    https://github.com/STScI-JWST/jwst/blob/master/jwst/refpix/reference_pixels.py#L690
    """
    if detector == 'NRS1':
        # NRS1 is just flipped over the line X=Y
        data = np.swapaxes(data, 2, 3)

    if detector == 'NRS2':
        # NRS2 is flipped over the line Y=X, then rotated 180 degrees
        data = np.swapaxes(data, 2, 3)[:, :, ::-1, ::-1]

    if detector == 'NRCA1':
        # NRCA1 is just flipped in X
        data = data[:, :, :, ::-1]

    if detector == 'NRCA2':
        # NRCA2 is just flipped in Y
        data = data[:, :, ::-1]

    if detector == 'NRCA3':
        # NRCA3 is just flipped in X
        data = data[:, :, :, ::-1]

    if detector == 'NRCA4':
        # NRCA4 is just flipped in Y
        data = data[:, :, ::-1]

    if detector == 'NRCALONG':
        # NRCA3 is just flipped in X
        data = data[:, :, :, ::-1]

    if detector == 'NRCB1':
        # NRCB1 is just flipped in Y
        data = data[:, :, ::-1]

    if detector == 'NRCB2':
        # NRCB2 is just flipped in X
        data = data[:, :, :, ::-1]

    if detector == 'NRCB3':
        # NRCB1 is just flipped in Y
        data = data[:, :, ::-1]

    if detector == 'NRCB4':
        # NRCB4 is just flipped in X
        data = data[:, :, :, ::-1]

    if detector == 'NRCBLONG':
        # NRCB1 is just flipped in Y
        data = data[:, :, ::-1]

    if detector == 'NIS':
        # NIRISS has a 180 degree rotation followed by a flip across the line
        # X=Y
        data = np.swapaxes(data[:, :, ::-1, ::-1], 2, 3)

    if detector == 'GUIDER1':
        # GUIDER1 is flipped in X and Y
        data = data[:, :, ::-1, ::-1]

    if detector == 'GUIDER2':
        # GUIDER2 is just flipped in X
        data = data[:, :, :, ::-1]

    # MIRI data doesn't need transforming

    return data

def detector_to_DMS(data, detector):
    if detector == 'NRS1':
        # Just flip back
        data = np.swapaxes(data, 2, 3)

    if detector == 'NRS2':
        # The inverse is to rotate 180 degrees, then flip over the line Y=X
        data = np.swapaxes(data[:, :, ::-1, ::-1], 2, 3)

    if detector == 'NRCA1':
        # Just flip back
        data = data[:, :, ::-1, ::-1]

    if detector == 'NRCA2':
        # Just flip back
        data = data[:, :, ::-1]

    if detector == 'NRCA3':
        # Just flip back
        data = data[:, :, :, ::-1]

    if detector == 'NRCA4':
        # Just flip back
        data = data[:, :, ::-1]

    if detector == 'NRCALONG':
        # Just flip back
        data = data[:, :, :, ::-1]

    if detector == 'NRCB1':
        # Just flip back
        data = data[:, :, ::-1]

    if detector == 'NRCB2':
        # Just flip back
        data = data[:, :, :, ::-1]

    if detector == 'NRCB3':
        # Just flip back
        data = data[:, :, ::-1]

    if detector == 'NRCB4':
        # Just flip back
        data = data[:, :, :, ::-1]

    if detector == 'NRCBLONG':
        # Just flip back
        data = data[:, :, ::-1]

    if detector == 'NIS':
        # Just flip and rotate back
        data = np.swapaxes(data, 2, 3)[:, :, ::-1, ::-1]
    
    if detector == 'GUIDER1':
        # Just flip back
        data = data[:, :, ::-1, ::-1]

    if detector == 'GUIDER2':
        # Just flip back
        data = data[:, :, :, ::-1]

    # MIRI data doesn't need transforming

    return data

def main(args):

    config = ConfigParser()
    config.read(args.config_file)

    old_hdulist = fits.open(args.input_file)

    new_hdulist = fits.HDUList()
    new_hdulist.append(fits.PrimaryHDU())
    new_hdulist[0].header = old_hdulist[0].header

    # get the exposure start and end times
    start_time = Time(old_hdulist[0].header['EXPSTART'], format='mjd').isot
    end_time = Time(old_hdulist[0].header['EXPEND'], format='mjd').isot

    params = {'sTime' : start_time, 'eTime' : end_time}

    s = requests.Session()

    # jwdmsdevwsvm1 is for testing. The actual DB host will be different.
    url_base = 'http://jwdmsdevwsvm1.stsci.edu/JWDMSEngSpAcc_CV2CV3/TlmMnemonicDataSrv.svc/Data/'

    for keyword, mnemonic in config['mnemonics'].items():

        # get request to server.
        url = url_base + mnemonic

        r = s.get(url, params=params, verify=False)

        # Parse json
        parsed_json = r.json()

        # json ObsTime has format '/Date(1358619814230+0000)/' which is 1358619814.230 in UNIX time
        # isotime = Time(float(parsed_json['Data'][0]['ObsTime'][6:-7])/1000., format='unix').isot

        # just take the first value of the series
        new_hdulist[0].header[keyword] = (parsed_json['Data'][0]['EUValue'], mnemonic.upper())

    # add the Engineering Mnemonics section heading
    new_hdulist[0].header.set('', 'Engineering Mnemonics', before=config['mnemonics'].keys()[0])
    new_hdulist[0].header.set('', '', before=config['mnemonics'].keys()[0])
    new_hdulist[0].header.set('', '', before=config['mnemonics'].keys()[0])

    # transform from DMS to detector orientation
    pixel_data = DMS_to_detector(old_hdulist['SCI'].data, old_hdulist['PRIMARY'].header['DETECTOR'])

    # collapse from 4D to 3D
    nints, ngroups, nx, ny = pixel_data.shape
    new_hdulist[0].data = pixel_data.reshape((nints*ngroups, nx, ny))
    new_hdulist[0].header.set('', '', before='BSCALE')

    # remove the NEXTEND keyword since there is only one extension now
    new_hdulist[0].header.remove('NEXTEND')

    # Write out
    new_hdulist.writeto(args.output_file, clobber=True)

if __name__ == '__main__':
    # Command line argument handler.
    parser = argparse.ArgumentParser(
        description='Convert JWST data from DMS format to FITSWriter format',
        epilog='example: flight2ground tlm.cfg input.fits output.fits')
    parser.add_argument('config_file', help='config file with Telemetry FITS keyword/mnemonic pairs')
    parser.add_argument('input_file', help='level 1B data file to reformat')
    parser.add_argument('output_file', help='name of output file')
    args = parser.parse_args()
    main(args)