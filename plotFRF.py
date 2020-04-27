#!/usr/bin/env python3

###########################################################
# class to plot the FRF from file or dictionary
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import OpenSSL.tsafe
import matplotlib
import matplotlib.pyplot as pl
import matplotlib.cm as mplcm
import matplotlib.colors as colors
from cycler import cycler
import scipy.io
import pickle
import numpy as np
import os


def changeExt(fileIn=None,newExt=None):
    """
    replace extension of filename
    """
    fileOut = fileIn
    if fileIn is not None and newExt is not None:
        spliTxt=os.path.splitext(fileIn)
        fileOut=spliTxt[0]+'.'+newExt
    return fileOut

class plotFRF:
    def __init__(self,fileIn=None,presRef=20e-6,dictIn=None,fileOut=None):
        matplotlib.interactive(True)

        self.plotCl = None
        self.ready = False
        self.freq = None
        self.frf = None
        self.frfGrad = None
        self.paraName = None
        self.gradName = None
        self.presRef = presRef
        self.paraVal = None
        self.lgd = None
        self.fig = None

        #load data from file
        if fileIn is not None:
            self.ready = self.prepDataFile(filename=fileIn)
        #
        #load data from dictionary
        if dictIn is not None:
            self.ready = self.prepDataDict(dictData=dictIn)
        #        
        #plot
        if self.ready:
            self.plot()
        # save to files
        if fileOut:
            self.exportData(changeExt(fileOut,'csv'))
            if self.ready:
                self.exportPlot(changeExt(fileOut,'tex'))
                self.exportPlot(changeExt(fileOut,'pdf'))

        #restore non intercative mode and replot
        matplotlib.interactive(False)
        if self.ready:
            self.plotCl.show()

    
    def normFRF(self,data,dataF=None,typeData=None):
        """
        function to normlized FRF (compute in dB) and derivatives
        """
        dataOut=data
        if self.presRef is not None:
            if typeData == 'f':
                dataOut=10*np.log10(np.array(data)/(self.presRef**2))
            if typeData == 'g':
                dataOut = np.array(data)*10*1/(np.array(dataF)**2)*np.log10(1/(self.presRef**2))*1/np.log(10)
        return dataOut

    def getUnit(self):
        """
        function to get the unit for plotting
        """
        if self.presRef is None:
            txtUnit = str(1)
        else:
            txtUnit = 'dB'
        return txtUnit

    def checkGrad(self,data):
        """
        function to check if the gradients data are avialable
        """
        gradData = False
        if data is not None:
            if len(data)>0:
                if len(data[0])>0:
                    gradData = True
        return gradData

    def prepDataFile(self,filename=None):
        """
        extract dictionary from save file
        """
        if filename is not None:
            #get extension of the file
            ext = os.path.splitext(filename)
            #
            if ext[1] == '.mat':
                data = scipy.io.loadmat(filename)
            if ext[1] == '.pck':
                f = open(filename,'rb')
                data = pickle.load(f)
            #read the dictionary and return
            return self.prepDataDict(data)


    def prepDataDict(self,dictData=None):
        """
        extract data from dictionary
        """
        if 'presRef' in dictData.keys():
            if dictData['presRef'] is not None:
                self.presRef = dictData['presRef']
        if 'frequencies' in dictData.keys():
            if dictData['frequencies'] is not None:
                if len(dictData['frequencies'])>0:
                    self.freq = dictData['frequencies']
        if 'FRF'  in dictData.keys():
            if dictData['FRF'] is not None:
                if len(dictData['FRF'])>0:
                    self.frf = dictData['FRF']
        if 'FRFgrad' in dictData.keys():
            if dictData['FRFgrad'] is not None:
                if self.checkGrad(dictData['FRFgrad']):
                    self.frfGrad = dictData['FRFgrad']
        if 'paraVal' in dictData.keys():
            if dictData['paraVal'] is not None:
                if len(dictData['paraVal'])>0:
                    # if isinstance(dictData['paraVal'],list):
                    self.paraVal = dictData['paraVal']
                    # else:
                    #     self.paraVal = [dictData['paraVal']]
        if 'name' in dictData.keys():
            if dictData['name'] is not None:
                if len(dictData['name'])>0:
                    self.paraName = dictData['name']
        if 'nameGrad' in dictData.keys():
            if dictData['nameGrad'] is not None:
                if len(dictData['nameGrad'])>0:
                    self.gradName = dictData['nameGrad']

        boolReady = self.freq is not None and self.frf is not None
        return boolReady 
    
    def loadVariables(self,freq=None,frf=None,frfGrad=None,paraName=None,gradName=None,paraVal=None):
        """
        load variables in the object
        """
        if freq is not None:
            self.freq = freq
        if frf is not None:
            self.frf = frf
        if self.checkGrad(frfGrad):
            self.frfGrad = frfGrad
        if paraName is not None:
            self.paraName = paraName
        if gradName is not None:
            self.gradName = gradName
        if paraVal is not None:
            if isinstance(paraVal,list):
                self.paraVal = paraVal
            else:
                self.paraVal = [paraVal]
    
    def getListLabels(self):
        """
        create list labels for the plots
        """
        listLabels = list()
        if self.paraVal is not None and self.paraName is not None:
            pVal=self.paraVal
            for itR in pVal:
                listLabels.append(" ".join(str(x)+"="+str(y) for (x,y) in zip(self.paraName,itR)))
        if self.paraVal is not None and self.paraName is None:
            pVal=self.paraVal
            for itR in pVal:
                listLabels.append(" ".join(str(x) for x in itR))
        if self.paraVal is None and self.paraName is None:
            for it in len(self.frf):
                listLabels.append("Case: "+str(it))
        return listLabels


    def plot(self,freq=None,frf=None,frfGrad=None,paraName=None,gradName=None,paraVal=None):
        """
        plot the figure
        """
        self.loadVariables(freq,frf,frfGrad,paraName,gradName,paraVal)
        #generate table of colors
        NUM_COLORS = len(self.frf)
        cm = pl.get_cmap('gist_rainbow')
        cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
        scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
        default_cycler = (cycler(color=[scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
        #
        self.plotCl = pl
        #
        nbGrad = 0
        nbRows = 1
        nbCols = 1
        if self.checkGrad(self.frfGrad):
            nbGrad = len(self.frfGrad[0])
            nbCols = 3 
            nbRows = int(np.ceil((nbGrad+1)/(nbCols-1)))
        #
        itSub = 0
        itCol = 0
        #open subplots
        self.fig,ax2d = self.plotCl.subplots(nbRows, nbCols, constrained_layout=False,squeeze=False,figsize=(6, 5), tight_layout=True)
        axs = ax2d   #.flatten()    # 1D linear indexing
        axs[itSub,itCol].set_prop_cycle(default_cycler)
        #
        listLabels = self.getListLabels()
        for it,val in enumerate(self.frf):
            l = axs[itSub,itCol].plot(self.freq.flatten(),self.normFRF(data = val,typeData = 'f').flatten(),label = listLabels[it])
            #
            axs[itSub,itCol].set_xlabel('Frequency (Hz)')
            axs[itSub,itCol].set_ylabel('H(f) [%s]'%self.getUnit())
            # axs[itSub,itCol].set_title('Mean quadratic pressure')
            axs[itSub,itCol].grid(True)
        #
        if nbCols>1:
            for itP in range(nbGrad):
                itSub += 1
                if itSub>=nbRows:
                    itCol += 1
                    itSub = 0
                axs[itSub,itCol].set_prop_cycle(default_cycler)
                for it,gradList in enumerate(self.frfGrad):
                    axs[itSub,itCol].plot(self.freq.flatten(),self.normFRF(gradList[itP,:].flatten(),self.frf[it].flatten(),'g'))#,label = listLabels[it])
                    #
                    axs[itSub,itCol].set_xlabel('Frequency (Hz)')
                    if self.paraName is not None:
                        axs[itSub,itCol].set_ylabel('dH(f)/d%s [%s/para]'%(self.gradName[itP],self.getUnit()))
                        # axs[itSub,itCol].set_title('Grad. wrt %s of mean quad. press.'%(self.gradName[itP]))
                    else:
                        axs[itSub,itCol].set_ylabel('dH(f)/dX%s [%s/para]'%(itP))
                        # axs[itSub,itCol].set_title('Grad. wrt X%s of mean quad. press.'%(itP))
                    axs[itSub,itCol].grid(True)
                    # axs[itSub,itCol].legend(loc=4)

        itSub+=1 
        #remove remaining subplot
        while itSub < nbRows:            
            axs[itSub,itCol].axis('off')
            itSub+=1
        for itSub in range(nbRows):
            axs[itSub,2].axis('off')

        self.lgd = self.fig.legend(loc = 'center right',bbox_to_anchor = (1.,0.5),ncol = 1,prop={'size': 6})
        # self.plotCl.figlegend( l, bbox_to_anchor = (1.,0.5), loc = 'center left' , ncol = 1)#, borderaxespad=0.1, labelspacing=0.,  prop={'size': 13} )
        self.plotCl.tight_layout()
        self.plotCl.subplots_adjust(right=0.2)
        self.plotCl.show()


    def exportData(self,fileOut=None):
        """
        export data as CSV
        """
        if self.freq is None or self.frf is None:
            import csv
            #prepare header
            headers = list()
            data = list()
            headers.append('#Freq')
            for it,val in enumerate(frf):
                headers.append('FRF_'+str(it)+'_'+self.getUnit())
            if self.checkGrad(self.frfGrad):
                for it,val in enumerate(frfGrad):
                    for itg,vall in enumerate(val):
                        if self.paraName is None:
                            headers.append('FRF_'+str(it)+'_G'+self.paraName[itg]+'_'+self.getUnit())
                        else:
                            headers.append('FRF_'+str(it)+'_G'+str(itg)+'_'+self.getUnit())
            #prepare data
            if self.checkGrad(self.frfGrad):
                data = np.vstack((self.freq,self.frf,self.frfGrad))
            else:
                data = np.vstack((self.freq,self.frf))
            data = data.T
            #start writing the file
            with open(fileOut) as fp:
                fp.write(';'.join(headers)+'\n')
                np.savetext(fp, data, '%g', ';')
        # write data about parameter sets
        if self.paraVal is not None:
            import csv
            ext = os.path.splitext(fileOut)
            fileOut = ext[0]+'_para'+ext[1]
            with open(fileOut,'w') as fp:
                fp.write(';'.join(self.paraName)+'\n')
                np.savetxt(fp, self.paraVal, '%g', ';')


    def exportPlot(self,fileOut=None):
        """
        export plot as Tikz
        """
        ext = os.path.splitext(fileOut)
        if self.plotCl is not None:
            if ext[1] == '.tex':
                try:
                    import tikzplotlib as tk 
                    tk.save(fileOut,figure=self.fig)
                except ModuleNotFoundError:
                    print('Tikzplotlib is not available')
            if ext[1] == '.png' or ext[1] =='.jpg' or ext[1] =='.jpeg' or ext[1] =='.pdf':
                if self.plotCl is not None:
                    self.fig.savefig(fileOut)#,bbox_extra_artists=(self.lgd), bbox_inches='tight')
        else:
            print('Unable to export picture')


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################

def usage():
    print("Usage: ", sys.argv[0], "file1 [file2 or [-ioh [+arg]]")
    print("\t file1: input file (.mat or .pck) mandatory")
    print("\t file2: optional output file (to export csv, plot, etc...)")
    print("\t -i + file1 : input file1")
    print("\t -o + file2 : input file2")

# function for dealing with options
def manageOpt(parser):
    """
    Manage options/arguments and run
    """
    fileIn = None
    fileOut = None
    #
    # parser.add_argument("file1", nargs=1, help="Input data file (mandatory)",default=None)
    # parser.add_argument("file2", nargs=1, help="Output file (optional)",default=None)
    parser.add_argument("-i", "--input",nargs=1, help="Input data file (mandatory)",required=True)
    parser.add_argument("-o", "--output",nargs=1, help="Output file (optional)",required=False)
    #
    args = parser.parse_args()
    #
    if args.input:
        fileIn = args.input[0]
    if args.output:
        fileOut = args.output[0]    
    #
    fileInOk = os.path.exists(fileIn)
    #
    if fileInOk:
        print('Input file: %s - status: %s'%(fileIn,"OK"))
    else:
        print('Input file: %s - status: %s'%(fileIn,"does not exist"))
    if fileOut:
        fileOutOk = not os.path.exists(fileOut)
        if fileOutOk:
            print('Output file: %s - status: Ok'%fileOut)
        else:
            print('Output file: %s - status: Already exists (will be overwritten)'%fileOut)

    # run computation
    if fileInOk:
        plotFRF(fileIn=fileIn,fileOut=fileOut)


# Run autonomous
if __name__ == '__main__':
    import sys
    import argparse
    #
    parser = argparse.ArgumentParser()
    # run with options
    manageOpt(parser)


