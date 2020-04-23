###########################################################
# class to plot the FRF from file 
# L. LAURENT - 2020 - luc.laurent@lecnam.net
###########################################################

import OpenSSL.tsafe
import matplotlib.pyplot as pl
import matplotlib.cm as mplcm
import matplotlib.colors as colors
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
        self.plot = None
        self.ready = False
        self.freq = None
        self.frf = None
        self.frfGrad = None
        self.paraName = None
        self.gradName = None
        self.presRef = presRef
        self.paraVal = None
        self.manyParametersSets = False

        #load data from file
        if fileIn is not None:
            self.ready = self.prepDataFile(filename=fileIn)
        #
        #load data from dictionary
        if dictIn is not None:
            self.ready = self.prepDataDict(dictIn=dictIn)
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


    def prepDataFile(self,filename=None):
        """
        extract dictionary from save file
        """
        if filename is not None:
            #get extension of the file
            ext = os.path.splitext(filename)
            #
            if ext[1] is '.mat':
                data = scipy.io.loadmat(filename)
            if ext[1] is '.pck':
                f = open(filename,'rb')
                data = pickle.load(f)
            #read the dictionary and return
            return self.prepDataDict(data)


    def prepDataDict(self,dictData=None):
        """
        extract data from dictionary
        """
        if 'frequencies' in dictData.keys():
            if dictData['frequencies'] is not None:
                if len(dictData['frequencies'])>0:
                    self.freq = dictData['frequencies']
        if 'FRF'  in dictData.keys():
            if dictData['FRF'] is not None:
                self.manyParametersSets = False
                if len(dictData['FRF'])>0:
                    if isinstance(dictData['FRF'],'list'):
                        self.frf = dictData['FRF']
                        self.manyParametersSets = True
                    else:
                        self.frf = [dictData['FRF']]
        if 'FRFgrad' in dictData.keys():
            if dictData['FRFgrad'] is not None:
                if len(dictData['FRFgrad'])>0:
                    self.frfGrad = dictData['FRFgrad']
        if 'paraVal' in dictData.keys():
            if dictData['paraVal'] is not None:
                if len(dictData['paraVal'])>0:
                    if isinstance(dictData['paraVal'],'list'):
                        self.paraVal = dictData['paraVal']
                    else:
                        self.paraVal = [dictData['paraVal']]
        if 'paraName' in dictData.keys():
            if dictData['paraName'] is not None:
                if len(dictData['paraName'])>0:
                    self.paraName = dictData['paraName']
        if 'gradName' in dictData.keys():
            if dictData['gradName'] is not None:
                if len(dictData['gradName'])>0:
                    self.gradName = dictData['gradName']

        boolReady = self.freq is not None and self.frf is not None
        return boolReady 
    
    def loadVariables(self,freq=None,frf=None,frfGrad=None,paraName=None,gradName=None,paraVal=None):
        """
        load variables in the object
        """
        if freq is not None:
            self.freq = freq
        if frf is not None:
            if isinstance(frf,'list'):
                self.manyParametersSets = True
                self.frf = frf
            else:
                self.frf = [frf]
        if frfGrad is not None:
            self.frfGrad = frfGrad
        if paraName is not None:
            self.paraName = paraName
        if gradName is not None:
            self.gradName = gradName
        if paraVal is not None:
            if isinstance(paraVal,'list'):
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
                listLabels.append(" ".join(str(x)+"="+str(y) for (x,y) in zip(self.paraName,pVal)))
        if self.paraVal is not None and self.paraName is None:
            pVal=self.paraVal
            for itR in pVal:
                listLabels.append(" ".join(str(x) for x in pVal))
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
        NUM_COLORS = len(frf)
        cm = pl.get_cmap('gist_rainbow')
        cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
        scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
        #
        self.plot = pl
        self.plot.figure(1)
        #
        nbGrad = 0
        nbRows = 1
        nbCols = 1
        if self.frfGrad is not None:
            nbGrad = self.frfGrad[0].shape[1]
            nbCols = 2 
            nbRows = np.ceil((nbGrad+1)/nbCols)
        #
        itSub = nbRows*100+nbCols*10+1
        ax = self.plot.subplot(itSub)
        ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
        #
        listLabels = self.getListLabels()
        for it,val in enumerate(self.frf):
            ax.plot(self.freq,val,label = listLabels[it])
            #
            ax.xlabel('Frequency (Hz)')
            ax.ylabel('Mean quadratic pressure (dB)')
            ax.grid(True)
            ax.legend(loc=4)
        #
        if nbCols>1:
            for itP in range(nbGrad):
                itSub += itP
                ax = self.plot.subplot(itSub)
                ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
                for it,gradList in enumerate(self.frfGrad):
                    self.plot.plot(self.freq,gradList[:,itP])#,label = listLabels[it]
                    #
                    ax.xlabel('Frequency (Hz)')
                    if self.paraName is not None:
                        ax.ylabel('Grad. wrt %s of mean quad. press. (dB)'%self.paraName(itP))
                    else:
                        ax.ylabel('Grad. wrt %s of mean quad. press. (dB)'%itP)
                    ax.grid(True)
                    ax.legend(loc=4)


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
                headers.append('FRF'+str(it))
            if self.frfGrad is not None:
                for it,val in enumerate(frfGrad):
                    for itg,vall in enumerate(val):
                        if self.paraName is None:
                            headers.append('FRF_'+str(it)+'_G'+self.paraName[itg])
                        else:
                            headers.append('FRF_'+str(it)+'_G'+str(itg))
            #prepare data
            if self.frfGrad is not None:
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
            fileOut = os.path.join(fileOut.replace(ext,''),'_para',ext)
            with open(fileOut) as fp:
                fp.write(';'.join(self.paraName)+'\n')
                np.savetext(fp, self.paraVal, '%g', ';')


    def exportPlot(self,fileOut=None):
        """
        export plot as Tikz
        """
        ext = os.path.splitext(fileOut)
        if ext[1] == '.tex':
            try:
                import tikzplotlib as tk 
                tk.save(fileOut,figure=)
            except ModuleNotFoundError:
                print('Tikzplotlib is not available')
        if ext[1] == '.png' or ext[1] =='.jpg' or ext[1] =='.jpeg' or ext[1] =='.pdf':
            if self.plot is not None:
                self.plot.savefig(fileOut)


