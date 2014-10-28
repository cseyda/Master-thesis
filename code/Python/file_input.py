"""
code for reading different file formats
output are generators with dictionaries
"""

from os.path import splitext

import hashlib
import codecs


def open_file(file_name):
    """
    returns a pair with the md5 of the file, and fitting generator for
    iterating the file.
    """
    #head,ext = splitext(file_name)
    md5 = hashlib.md5()
    #print head, ext
    with open(file_name,'rb') as f: 
        for chunk in iter(lambda: f.read(128*md5.block_size), b''): 
            md5.update(chunk)
    
    #print md5.hexdigest()
    if file_name.endswith(".gt"):
        return (md5.hexdigest(), gt_generator(file_name))
    elif file_name.endswith(".egt"):
        return (md5.hexdigest(), egt_generator(file_name))
    
def gt_generator(file):
    with codecs.open(file, "r", "utf-8") as f:
        for line in f:
            tmp = line.split()
            yield {"loc" : [float(tmp[1]), float(tmp[2])],#[lon, lat]
                "tag" : tmp[3:]}

def egt_generator(file):
    with codecs.open(file, "r", "utf-8") as f:
        for line in f:
            tmp = line.split()
            yield {"loc" : [float(tmp[1]), float(tmp[2])],
                "id" : int(tmp[3]),
                "tag" : tmp[4:]}
