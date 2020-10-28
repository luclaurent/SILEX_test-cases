###########################################################
# USEFUL VARIOUS TOOLS 
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

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
    if len(bbxF)==6 and coordinateIn.shape[1] == 3:
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
# compute gradient of norm of a complex number
def computeNormGradComplexField(field,gradfield):
    return (np.real(gradfield)+np.imag(gradfield))/np.absolute(field)


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

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# distribution of frequencies per processor
def computeFreqPerProc(nbStep, nbProc, freqInit, freqEnd):
    # compute integer number of freq per proc and remaining steps
    nbFreqProc = nbStep // nbProc
    nbFreqProcRemain = nbStep % nbProc
    # compute frequencies steps
    varCase = 1
    if nbFreqProcRemain == 0:
        varCase = 0
    listFreq = np.zeros((nbFreqProc+varCase, nbProc))
    listAllFreq = np.linspace(freqInit, freqEnd, nbStep)
    # print(np.linspace(freqInit,freqEnd,nbStep))
    # build array of frequencies
    itF = 0
    for itP in range(nbProc):
        for itC in range(nbFreqProc+varCase):
            if itC*nbProc+itP < nbStep:
                listFreq[itC, itP] = listAllFreq[itF]
                itF += 1
    # print(listFreq)
    return listFreq    
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Useful tools to replace None in comlex dictionary
class Null(object):
    def __repr__(self):
        return 'Null'

NULL = Null()

def replace_none(data):
    for k, v in data.items() if isinstance(data, dict) else enumerate(data):
        if v is None:
            data[k] = NULL
        elif isinstance(v, (dict, list)):
            replace_none(v)

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# add extension to filename if necessary
def addExt(filename,ext=None):
    outputFilename = filename
    if outputFilename is not None:        
        if ext is not None:
            if ext[0] is not '.':
                ext = '.'+ext
            if outputFilename[len(outputFilename)-len(ext):len(outputFilename)] != ext:
                outputFilename += ext
    return outputFilename    

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# class used to redirect Fortran stdout/stderr
class RedirectFortran:
    def __init__(self, outputfile=os.devnull):
        self.outfiles = outputfile

    def start(self):
        # open outputfile
        self.out = os.open(self.outfiles,os.O_RDWR|os.O_CREAT)
        # save the current file descriptor
        self.save = os.dup(1)
        # put outfile on 1
        os.dup2(self.out, 1)

    def stop(self, funExport = None):
        # restore the standard output file descriptor
        os.dup2(self.save, 1)
        # close the output file
        os.close(self.out)
        # read it
        file = open(self.outfiles,'r')
        #
        if not funExport:
            for line in file: 
                print(line.replace('\n',''))
        else:
            for line in file: 
                funExport(line.replace('\n',''))
        # delete it
        if self.outfiles != os.devnull:
            os.remove(self.outfiles)
        
