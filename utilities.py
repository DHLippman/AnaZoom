"""
Author:         David Henry Lippman
File:           utilities.py
Date created:   07/22/20
Date modified:  07/22/20

"""

import numpy as np
import os
import datetime
import pickle


def rand_rng(min_val, max_val, sign=1):

    """


    Returns a random value within a range with the desired sign


    min_val:        minimum value of range (absolute value)

    max_val:        maximum value of range (absolute value)

    sign:           the sign of the random value

                    1 for positive
                   -1 for negative

    """

    return sign * np.random.uniform(np.abs(min_val), np.abs(max_val))


def find_roc(efl, n):

    """


    Finds the radii of curvature to form a equiconvex/equiconcave thin lens of
    a given focal length.


    efl:            thin lens effective focal length [mm]

    n:              refractive index of thin lens

    """

    R1 = 2 * efl * (n - 1)
    R2 = -R1

    return R1, R2


def get_time_str():

    """

    Gets a timestamp in the form of a string

    """

    # Convert to string

    timeStr = str(datetime.datetime.now())

    # Replace invalid characters

    timeStr = timeStr.replace('-', '_')
    timeStr = timeStr.replace(' ', '_')
    timeStr = timeStr.replace(':', '_')

    # Remove seconds decimal places

    timeStr = timeStr[0: timeStr.find('.')]

    # Put year after month and day

    timeStr = timeStr[5: 11] + timeStr[0: 5] + timeStr[11:]

    return timeStr


def save_obj(obj, filename=None):

    """


    Saves a Python object to a file. Used to be loaded at a future time to cut
    down on repeated computation time


    obj:        Python object to save

    filename:   filename to save the object with; optional

    """

    if filename is None:
        filename = get_time_str()

    path = 'src/obj/' + filename

    print('\nSaving object:  ' + path +'\n')

    with open(path, 'wb') as fid:
        pickle.dump(obj, fid)

    return obj


def load_obj(filename=None):

    """


    Loads a saved Python object from a file. Used to cut down on repeated
    computation time


    filename:   filename containing object to load

    """

    # If file name isn't provided, print list of options

    if filename is None:

        # Print list of targets to analyze

        path = 'src/obj/'
        filename_list = os.listdir(path)
        for fid in filename_list:
            if fid[0] == '.':
                filename_list.remove(fid)  # delete hidden files from list
        template = '\t{0:0d}\t{1}'
        print()
        for i, fid in enumerate(filename_list):
            print(template.format(i + 1, fid))

        # User selects desired target

        result = input('\nEnter the object you would like to load (#): ')
        ind = eval(result) - 1
        filename = filename_list[ind]

    # Create path

    path = 'src/obj/' + filename

    # Load object(s)

    with open(path, 'rb') as fid:
        obj = pickle.load(fid)

    print('\nLoading object:  ' + path + '\n')

    # Add filename as attribute to object(s)

    if type(obj) is list:
        for o in obj:
            o.filename = filename
    else:
        obj.filename = filename

    return obj


def folder_exist(path, erase=False):

    """


    Checks to see if folder exists and creates the folder if it does not exist


    path:           path of folder to check/create

    erase:          flag for erasing current folder contents

    """

    # Add folder if it does not exists already

    if not os.path.exists(path):
        os.makedirs(path)

    # Delete folder and contents, if desired

    if erase:
        for file in os.listdir(path):
            filePath = os.path.join(path, file)
            if os.path.isfile(filePath):
                os.unlink(filePath)
        os.rmdir(path)


def format_time_str(t):

    """


    Returns a time string formatted for the most appropriate unit


    t:      time in seconds

    """

    # Determine best time unit

    if t > 3600 * 24:
        return '{0:0.2f} days'.format(t / 3600 / 24)
    elif t > 3600:
        return '{0:0.2f} hours'.format(t / 3600)
    elif t > 60:
        return '{0:0.2f} minutes'.format(t / 60)
    else:
        return '{0:0.2f} seconds'.format(t)
