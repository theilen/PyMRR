# -*- coding: utf-8 -*-
"""
Extract sequence parameters from dicom header.

@author: Sebastian Theilenberg
"""

__version__ = '0.81'
# $Source$


import dicom

from .siemens_csa import parse_csa_header, parse_protocol_data
from ..bvalue import trapezoid_G


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
    # pause after the trigger
    "sWiPMemBlock.alFree[6]": ["POST", lambda i: float(i)*1e-3]
}


def read_parameters(dicom_file):
    dc = dicom.read_file(dicom_file, stop_before_pixels=True)
    return parse_parameters(dc)


def parse_parameters(dicom_data):
    """
    Parse specific fields out of the dicom-files CSA header.

    Returns a dictionary containing (name : value) pairs of the following
    parameters:
        delta : float
            length of motion sensitizing gradient  in ms
        Delta : float
            time between motion sensitizing gradients
        gradient start : float
            time from the maximum of the excitation-pulse to the start of the
            first gradient in ms
        bvalue : float
            specified b-value in s/mm^2
        maxGrad : float
            gradient's strength in mT/m
        POST : float
            time the trigger signal got shifted towards the excitation in ms
    """
    # Get private CSA header from dicom-file and parse it.
    # This should be tag (0x0029, 0x1020), but may be one of the following,
    # too (Actually the private header seems to be present multiple times in
    # the dicom header): (0x0029, 0x1010), (0x0029, 0x1210), (0x0029, 0x1110),
    # (0x0029, 0x1220), (0x0029, 0x1120)
    for tag in [(0x0029, 0x1020), (0x0029, 0x1120), (0x0029, 0x1220),
                (0x0029, 0x1010), (0x0029, 0x1110), (0x0029, 0x1210)
                ]:
        data = dicom_data[tag].value
        if data:
            csa = parse_csa_header(data)
            if "MrPhoenixProtocol" in csa.keys():
                break
    assert csa

    # parse MrPhoenixProtocol, that contains the magic
    mrp = parse_protocol_data(csa["MrPhoenixProtocol"])
    assert mrp

    parameters = {}
    for tag, specifier in _PHOENIX_TAGS.items():
        name, func = specifier
        try:
            parameters[name] = func(mrp[tag])
        except KeyError:
            parameters[name] = None

    # currently sWiPMemBlock.alFree[6] is only set if unequal to zero
    if parameters["POST"] is None:
        parameters["POST"] = 0.

    # (Re-)calculate parameters
    parameters["Delta"] += parameters["delta"]  # to match Bernstein
    parameters["maxGrad"] = trapezoid_G(parameters['bvalue']*1e6,
                                        parameters['delta']*1e-3,
                                        parameters['Delta']*1e-3)*1e3

    return parameters
