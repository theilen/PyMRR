# -*- coding: utf-8 -*-
"""
Extract sequence parameters from dicom header.

@author: Sebastian Theilenberg
"""

__version__ = '0.82'
# $Source$


import dicom

from .siemens_csa import get_phoenix_protocol


# tags to be read out of MrPhoenixProtocol and their corresponding names
_PHOENIX_TAGS = {
    # bvalue
    "sDiffusion.alBValue[0]": ["bvalue", lambda i: float(i)],
    # length of the gradients (one ramp + plateau)
    "sWiPMemBlock.adFree[7]": ["delta", lambda i: float(i)],
    # time between two gradients (end to start)
    "sWiPMemBlock.adFree[8]": ["Delta", lambda i: float(i)],
    # start of the first gradient after the maximum of the excitation pulse
    "sWiPMemBlock.adFree[9]": ["gradstart", lambda i: float(i)],
    # gradient amplitude
    "sWiPMemBlock.adFree[12]": ["maxGrad", lambda i: float(i)],
    # post Trigger Fill Time (ms)
    "sWiPMemBlock.alFree[6]": ["PTFT", lambda i: float(i)*1e-3],
    # PFTF decrement (ms)
    "sWiPMemBlock.alFree[20]": ["PTFT_decr", lambda i: float(i)*1e-3],
    # PTFT averages
    "sWiPMemBlock.alFree[21]": ["PTFT_aver", lambda i: int(i)],
    # number of images (repetitions = images-1)
    "lRepetitions": ["TotalImages", lambda i: int(i)+1],
}


def parse_parameters(dcm):
    """
    Parse specific fields out of the dicom-files CSA header.

    Returns a dictionary containing (name : value) pairs of the following
    parameters:
        protocol : str
            the protocol name of the measurement
        delta : float
            length of motion sensitizing gradient  in ms
        Delta : float
            time between motion sensitizing gradients
        gradstart : float
            time from the maximum of the excitation-pulse to the start of the
            first gradient in ms
        bvalue : float
            specified b-value in s/mm^2
        maxGrad : float
            gradient's strength in mT/m
        PTFT : float
            time the trigger signal got shifted towards the excitation in ms
        PTFT_aver : int (if present)
            number of images per PTFT
        PTFT_decr : float (if present)
            decrement in PTFT
        echotime : float
            echo time in ms
    """
    dcm = _check_for_dicom_data(dcm)

    parameters = {}
    parameters["protocol"] = dcm.ProtocolName
    parameters["echotime"] = dcm.EchoTime

    if not check_sequence(dcm):
        return parameters

    mrp = get_phoenix_protocol(dcm)

    for tag, specifier in _PHOENIX_TAGS.items():
        name, func = specifier
        try:
            parameters[name] = func(mrp[tag])
        except KeyError:
            pass

    # currently sWiPMemBlock.alFree[6] is only set if unequal to zero
    if "PTFT" not in parameters:
        parameters["PTFT"] = 0.

    # (Re-)calculate parameters
    if "Delta" in parameters:
        parameters["Delta"] += parameters["delta"]  # to match Bernstein
    # recalculate PTFT in case of variable seq
    if "PTFT_aver" in parameters and "PTFT_decr" in parameters:
        parameters["PTFT"] = _calc_ptft(dcm.InstanceNumber,
                                        parameters["PTFT"],
                                        parameters["PTFT_aver"],
                                        parameters["PTFT_decr"])

    return parameters


def _check_for_dicom_data(dcm, header_only=True):
    if not isinstance(dcm, dicom.dataset.Dataset):
        dcm = dicom.read_file(dcm, stop_before_pixels=header_only)
    return dcm


def variable_ptft(dcm):
    "Check, whether dcm is a dicom with variable PTFT."
    dcm = _check_for_dicom_data(dcm)
    if "PTFT_aver" in get_phoenix_protocol(dcm):
        return True
    else:
        return False


def check_sequence(dcm):
    dcm = _check_for_dicom_data(dcm)
    return dcm.ProtocolName in ["nin_ep2d_diff_vb10r"]


def _calc_ptft(index, fill, aver, decr):
    "Calculate individual ptft."
    if index < 2:
        return fill
    t = (index-2)//aver
    return fill - t*decr
