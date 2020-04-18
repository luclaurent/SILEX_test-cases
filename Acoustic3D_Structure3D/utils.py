###########################################################
# USEFUL VARIOUS TOOLS 
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import importlib
import numpy as np
import os

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# function to import module as global
def global_imports(modulename,shortname = None, asfunction = False):
    if shortname is None: 
        shortname = modulename
    if asfunction is False:
        globals()[shortname] = __import__(modulename)
    else:        
        globals()[shortname] = eval(modulename + "." + shortname)

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# load communication to one processor for MUMPS
class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# function to prepare string (if empty)
def prepareStr(strIn):
    if strIn:
        strOut=strIn
    else:
        strOut='None'
    return strOut

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# function to get the list of nodes from a bbx declaration
def getNodesBBX(coordinateIn,bbxIn):
    #two kind of bbx declaration (list or 2D array)
    bbxTmp=np.array(bbxIn)
    if len(bbxTmp.shape)==2:
        bbxF=bbxTmp.flatten()
    if len(bbxTmp.shape)==1:
        bbxF=bbxTmp
    #extract boundaries
    xmin=bbxF[0]
    xmax=bbxF[1]
    ymin=bbxF[2]
    ymax=bbxF[3]
    flag3d=False
    if len(bbxF)==6:
        flag3d=True
        zmin=bbxF[4]
        zmax=bbxF[5]
    # find nodes
    epsM=10*np.finfo(float).eps
    if flag3d:
        nodeslist=np.where(
            (coordinateIn[:,0]>xmin-epsM)*
            (coordinateIn[:,0]<xmax+epsM)*
            (coordinateIn[:,1]>ymin-epsM)*
            (coordinateIn[:,1]<ymax+epsM)*
            (coordinateIn[:,2]>zmin-epsM)*
            (coordinateIn[:,2]<zmax+epsM))
    else:
        nodeslist=np.where(
            (coordinateIn[:,0]>xmin-epsM)*
            (coordinateIn[:,0]<xmax+epsM)*
            (coordinateIn[:,1]>ymin-epsM)*
            (coordinateIn[:,1]<ymax+epsM))
    #
    return np.array(nodeslist).flatten()

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# check if file exists
def checkFile(filename = None,typeCheck=None):
    flagStatus = filename
    if filename is not None:
        flagStatus = os.path.exists(filename)
        if typeCheck is 'file':
            flagStatus = os.path.isfile(filename)
        if typeCheck is 'dir':
            flagStatus = os.path.isdir(filename)
    return flagStatus

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# tools to obtain the size of a file
def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        num
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0


def file_size(file_path):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        return convert_bytes(file_info.st_size)

    