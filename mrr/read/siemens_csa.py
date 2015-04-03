# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:42:12 2015

Functions to parse Siemens private dicom CSA header.

Code in this file has been taken from siemens.py at the VeSPA project
(trunk.common.util.dicom.siemens.py rev 2566)
taken from:
https://scion.duhs.duke.edu/vespa/project/browser/trunk/common/util/dicom/siemens.py
on Feb 18 2015.
I only modified the code slighlty to fit my needs, all credit on how to parse
these CSA headers therefore go to the original author.

@author: Sebastian Theilenberg
"""


__version__ = '0.81'
# $Source$


import struct


def get_phoenix_protocol(dicom_data):
    "Get the parsed 'MrPhoenixProtocol' data of dicom_data"
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
            csa = _parse_csa_header(data)
            if "MrPhoenixProtocol" in csa.keys():
                break
    assert csa

    # parse MrPhoenixProtocol, that contains the magic
    mrp = _parse_protocol_data(csa["MrPhoenixProtocol"])
    assert mrp

    return mrp


def _parse_csa_header(tag, little_endian=True):
    """The CSA header is a Siemens private tag that should be passed as
    a string. Any of the following tags should work: (0x0029, 0x1010),
    (0x0029, 0x1210), (0x0029, 0x1110), (0x0029, 0x1020), (0x0029, 0x1220),
    (0x0029, 0x1120).

    The function returns a dictionary keyed by element name.
    """
    # Let's have a bit of fun, shall we? A Siemens CSA header is a mix of
    # binary glop, ASCII, binary masquerading as ASCII, and noise masquerading
    # as signal. It's also undocumented, so there's no specification to which
    # to refer.

    # The format is a good one to show to anyone who complains about XML being
    # verbose or hard to read. Spend an afternoon with this and XML will
    # look terse and read like a Shakespearean sonnet.

    # The algorithm below is a translation of the GDCM project's
    # CSAHeader::LoadFromDataElement() inside gdcmCSAHeader.cxx. I don't know
    # how that code's author figured out what's in a CSA header, but the
    # code works.

    # I added comments and observations, but they're inferences. I might
    # be wrong. YMMV.

    # Some observations --
    # - If you need to debug this code, a hexdump of the tag data will be
    #   your best friend.
    # - The data in the tag is a list of elements, each of which contains
    #   zero or more subelements. The subelements can't be further divided
    #   and are either empty or contain a string.
    # - Everything begins on four byte boundaries.
    # - This code will break on big endian data. I don't know if this data
    #   can be big endian, and if that's possible I don't know what flag to
    #   read to indicate that. However, it's easy to pass an endianness flag
    #   to _get_chunks() should the need to parse big endian data arise.
    # - Delimiters are thrown in here and there; they are 0x4d = 77 which is
    #   ASCII 'M' and 0xcd = 205 which has no ASCII representation.
    # - Strings in the data are C-style NULL terminated.

    # I sometimes read delimiters as strings and sometimes as longs.
    DELIMITERS = ("M", "\xcd", 0x4d, 0xcd)

    # This dictionary of elements is what this function returns
    elements = {}

    # I march through the tag data byte by byte (actually a minimum of four
    # bytes at a time), and current points to my current position in the tag
    # data.
    current = 0

    # The data starts with "SV10" followed by 0x04, 0x03, 0x02, 0x01.
    # It's meaningless to me, so after reading it, I discard it.
    size, chunks = _get_chunks(tag, current, "4s4s")
    current += size

    assert chunks[0] == "SV10"
    assert chunks[1] == "\4\3\2\1"

    # get the number of elements in the outer list
    size, chunks = _get_chunks(tag, current, "L")
    current += size
    element_count = chunks[0]

    # Eat a delimiter (should be 0x77)
    size, chunks = _get_chunks(tag, current, "4s")
    current += size
    assert chunks[0] in DELIMITERS

    for i in range(element_count):
        # Each element looks like this:
        # - (64 bytes) Element name, e.g. ImagedNucleus, NumberOfFrames,
        #   VariableFlipAngleFlag, MrProtocol, etc. Only the data up to the
        #   first 0x00 is important. The rest is helpfully populated with
        #   noise that has enough pattern to make it look like something
        #   other than the garbage that it is.
        # - (4 bytes) VM
        # - (4 bytes) VR
        # - (4 bytes) syngo_dt
        # - (4 bytes) # of subelements in this element (often zero)
        # - (4 bytes) a delimiter (0x4d or 0xcd)
        size, chunks = _get_chunks(tag, current,
                                   "64s" + "4s" + "4s" + "4s" + "L" + "4s")
        current += size

        name, vm_, vr_, syngo_dt_, subelement_count, delimiter = chunks
        assert delimiter in DELIMITERS

        # The subelements hold zero or more strings. Those strings are stored
        # temporarily in the values list.
        values = []

        for j in range(subelement_count):
            # Each subelement looks like this:
            # - (4 x 4 = 16 bytes) Call these four bytes A, B, C and D. For
            #   some strange reason, C is always a delimiter, while A, B and
            #   D are always equal to one another. They represent the length
            #   of the associated data string.
            # - (n bytes) String data, the length of which is defined by
            #   A (and A == B == D).
            # - (m bytes) Padding if length is not an even multiple of four.
            size, chunks = _get_chunks(tag, current, "4L")
            current += size

            assert chunks[0] == chunks[1]
            assert chunks[1] == chunks[3]
            assert chunks[2] in DELIMITERS
            length = chunks[0]

            # get a chunk-o-stuff, length indicated by code above.
            # Note that length can be 0.
            size, chunks = _get_chunks(tag, current, "%ds" % length)
            current += size
            if chunks[0]:
                values.append(chunks[0])

            # If we're not at a 4 byte boundary, move.
            # Clever modulus code below swiped from GDCM
            current += (4 - (length % 4)) % 4

        # The value becomes a single string item (possibly "") or a list
        # of strings
        if len(values) == 0:
            values = ""
        if len(values) == 1:
            values = values[0]

        assert name not in elements
        elements[name] = values

    return elements


def _null_truncate(s):
    """Given a string, returns a version truncated at the first '\0' if
    there is one. If not, the original string is returned."""
    i = s.find(chr(0))
    if i != -1:
        s = s[:i]

    return s


def _scrub(item):
    """Given a string, returns a version truncated at the first '\0' and
    stripped of leading/trailing whitespace. If the param is not a string,
    it is returned unchanged."""
    if isinstance(item, basestring):
        return _null_truncate(item).strip()
    else:
        return item


def _get_chunks(tag, index, format, little_endian=True):
    """Given a CSA tag string, an index into that string, and a format
    specifier compatible with Python's struct module, returns a tuple
    of (size, chunks) where size is the number of bytes read and
    chunks are the data items returned by struct.unpack(). Strings in the
    list of chunks have been run through _scrub().
    """
    # The first character of the format string indicates endianness.
    format = ('<' if little_endian else '>') + format
    size = struct.calcsize(format)
    chunks = struct.unpack(format, tag[index:index + size])

    chunks = [_scrub(item) for item in chunks]

    return (size, chunks)


def _parse_protocol_data(protocol_data):
    """Returns a dictionary containing the name/value pairs inside the
    "ASCCONV" section of the MrProtocol or MrPhoenixProtocol elements
    of a Siemens CSA Header tag.
    """
    # Protocol_data is a large string (e.g. 32k) that lists a lot of
    # variables in a JSONish format with which I'm not familiar. Following
    # that there's another chunk of data delimited by the strings you see
    # below.
    # That chunk is a list of name=value pairs, INI file style. We
    # ignore everything outside of the ASCCONV delimiters. Everything inside
    # we parse and return as a dictionary.
    start = protocol_data.find("### ASCCONV BEGIN ###")
    end = protocol_data.find("### ASCCONV END ###")

    assert start != -1
    assert end != -1

    start += len("### ASCCONV BEGIN ###")
    protocol_data = protocol_data[start:end]

    lines = protocol_data.split('\n')

    # The two lines of code below turn the 'lines' list into a list of
    # (name, value) tuples in which name & value have been stripped and
    # all blank lines have been discarded.
    f = lambda pair: (pair[0].strip(), pair[1].strip())
    lines = [f(line.split('=')) for line in lines if line]

    return dict(lines)
