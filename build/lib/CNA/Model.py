# Import standard Python and NumPy modules.
import sys
import numpy as np
import matplotlib.pyplot as plt


#import glob and os
import glob
import os
import datetime

#import svc from sklearn
from sklearn import svm
from sklearn.externals import joblib

from CNA import Utilities

class Model_Maker():
    """
    Robert:
        
        This class of functions will allow the user to train a Support Vector Machine
        (SVM) model based on pre-loaded CNA Patterns.
        
        One may also retrain an existing model with new data/structures if one finds
        that the current model is insufficient.
        
        In general, the model SHOULD be independent of elemental composition
        
        A S S U M I N G
        
        that the user has set a reasonable value for R_Cut when loading new data in.
        
        At present, the Parameters argument for the function will contain information on :
            
            MasterKey - List of tuples - Basically, should the default be used or does
            the user wish to use their own? 
            One is advised to leave this argument as None and use the default.
            
            XYZ_Path - String - This will be the path from the working directory in which
            the code has been run to where one may find a directory of *.xyz files.
            Leaving this as None or empty will mean that no new *.xyz files will be used 
            to (re)train the SVM model.
    """
    
    def __init__(self, Parameters = None):

        self.script_path = os.path.dirname(os.path.realpath(__file__))+'/'
        
        self.cwd = os.getcwd()
        

        os.chdir(self.script_path)
        os.chdir('../')
        
        if os.path.exists('Training_Info.txt'):
            with open('Training_Info.txt', 'a') as f:
                
                f.write(
                        '\nTraining the support vector machine model at %s\n'
                        %( 
                            datetime.datetime.now().strftime("%a %d %b %Y %H:%M:%S") 
                            ) 
                        )
                
                f.close()
        else:
            #from Utilities
            with open('Training_Info.txt', 'a')as f:
                f.write(Utilities.Logo().Sapphire())
                f.write('\n')
                f.write("This is a log created by %s on %s to detail the training " 
                        "which has gone into creating a support vector machine model "
                        "\nfor the characterisation of metallic nanoparticles.\n")
                f.write('\nMeaning of values from model prediction:\n')
                f.write('0: Face Centered Cubic structure (To, Co, Oc)\n')
                f.write('1: Icosahedral structure (Ih)\n')
                f.write('2: Decahedral structure (InoDh, mDh (x,x,x))\n')
                f.write('3: Amorphous structure (Am)\n')                
                f.write('\nTraining the support vector machine model at %s\n'
                        %datetime.datetime.now().strftime("%a %d %b %Y %H:%M:%S"))
                f.close()
      
        self.Pattern_Dict = np.load('CNA_npz/pattern_dictionary.npz', allow_pickle=True)
        #This is a path that is intended to be constant throughout the development of Sapphire.
    
        os.chdir(self.cwd)

        if Parameters is not None:
            self.Parameters = Parameters
            with open('Training_Info.txt', 'a')as f:
                f.write('Reading user defined input.\n')
                f.close()

            if self.Parameters['MasterKey'] is None:
                self.MasterKey = Utilities.CNA_Masterkey().Key()
                with open('Training_Info.txt', 'a')as f:
                    f.write('Using default CNA signature masterkey.\n')
                    f.close()
                
            else:
                self.MasterKey = self.Parameters['MasterKey']
                with open('Training_Info.txt', 'a')as f:
                    f.write("Using user's own CNA signature masterkey.\n")
                    f.write("Do so at your own risk.\n")
                    f.write("\n%s\n" %self.MasterKey)
                    f.close()
            
            if self.Parameters['XYZ_Path'] is None:
                self.NewData = False
                with open('Training_Info.txt', 'a')as f:
                    f.write("Will not be reading in new *.xyz files for the model.\n")
                    f.close()
                
            else:
                self.XYZ = self.Parameters['XYZ_Path']
                with open('Training_Info.txt', 'a')as f:
                    f.write("Will be reading in new XYZ files from %s.\n" %self.XYZ)
                    f.close()                
            
        else:
            self.MasterKey = Utilities.CNA_Masterkey().Key()
            with open('Training_Info.txt', 'a')as f:
                f.write('Using default CNA signature masterkey.\n')
                f.close()
            
        self.Bulk_Pattern_Rows=[]
        
        with open('Training_Info.txt', 'a')as f:
            f.write('\nColumn Number of bulk patterns:\n')
            f.close()
            
        if [((5,5,5),12)] in list(self.Pattern_Dict['masterkey']):
            self.icosahedral=list(self.Pattern_Dict['masterkey']).index([((5,5,5),12)])
            with open('Training_Info.txt', 'a')as f:
                f.write('Icosahedral:\t %s' %self.icosahedral)
                f.close()
            self.Bulk_Pattern_Rows.append(self.icosahedral)
            
        if [((4,2,2),10),((5,5,5),2)] in list(self.Pattern_Dict['masterkey']):
            self.twinning_lines=list(self.Pattern_Dict['masterkey']).index([((4,2,2),10),((5,5,5),2)])
            with open('Training_Info.txt', 'a')as f:            
                f.write('HCP lines:\t %s' %self.twinning_lines)
                f.close()
            self.Bulk_Pattern_Rows.append(self.twinning_lines)
            
        if [((4,2,1),6),((4,2,2),6)] in list(self.Pattern_Dict['masterkey']):
            self.twinning_planes=list(self.Pattern_Dict['masterkey']).index([((4,2,1),6),((4,2,2),6)])
            with open('Training_Info.txt', 'a')as f:
                print('HCP planes:\t %s'  %self.twinning_planes)
                f.close()
            self.Bulk_Pattern_Rows.append(self.twinning_planes)
            
        if [((4,2,1),12)] in list(self.Pattern_Dict['masterkey']):
            self.FCC=list(self.Pattern_Dict['masterkey']).index([((4,2,1),12)])
            print('FCC:\t\t %s' %self.FCC)
            self.Bulk_Pattern_Rows.append(self.FCC)
            
        with open('Training_Info.txt', 'a')as f:
            f.writelines(self.Bulk_Pattern_Rows)
            f.close()
        
        #ACCESSING THE COORDINATION NUMBER OF ALL PATTERNS
        k=[]
        for i in range(len(self.Pattern_Dict['masterkey'])):
            #print(Pattern_Dict['masterkey'][i])
            c = 0
            for j in range(len(self.Pattern_Dict['masterkey'][i])):
                c += self.Pattern_Dict['masterkey'][i][j][1]
            k.append(c)
            
        self.Coordination_Number=np.asarray(k)
               
        
##############################################################################
        

    def Train(self, Filepath = '', ModelName = 'Model.pkl'):
        
        #RUNNING THROUGH ALL FILES IN THE XYZ DIRECTORY
        self.Structure_List = []
        self.Known_Filename = []
        self.Unknown_Filename = []
        self.Known_Pattern = []
        self.Unknown_Pattern = []
        self.Known_NAtom = []
        self.Unknown_NAtom = []
        self.folder_path = self.System['Training_Data']

        for filename in Pattern_Dict.files:
            if(filename == 'masterkey'):
                continue
        
            Unknown=False
            #CHECKING THE FILENAME AND CATEGORIZING THE XYZ FILES
            if any(x in filename for x in ['To','Co','Oc']):   #FCC
                self.Structure_List.append(0)
            elif any(x in filename for x in ['Ih']):           #IH
                self.Structure_List.append(1)
            elif any(x in filename for x in ['Dh']):           #DH
                self.Structure_List.append(2)
            elif any(x in filename for x in ['Am']):           #AM
                self.Structure_List.append(3)
            else:                                                               #UNKNOWN
                Unknown=True
        
            #APPENDING THE PATTERN DATA TO KNOWN OR UNKNOWN LISTS
            if(Unknown == True):
                continue
        
            if(Unknown == False):
                self.Known_Pattern.append(np.sum(self.Pattern_Dict[filename],axis=0))
                self.Known_NAtom.append(len(self.Pattern_Dict[filename]))
                self.Known_Filename.append(filename)
        
        self.Structures_Arrays = np.asarray(self.Structure_List)
        self.Known_Pattern_Arrays = np.asarray(self.Known_Pattern)
        self.Unknown_Pattern_Arrays = np.asarray(self.Unknown_Pattern)
        
        self.FCC_indices = np.nonzero(self.Structures_Arrays==0)
        self.Ih_indices = np.nonzero(self.Structures_Arrays==1)
        self.Dh_indices = np.nonzero(self.Structures_Arrays==2)
        self.Am_indices = np.nonzero(self.Structures_Arrays==3)
        
        self.Bulk_Array = np.take(self.Known_Pattern_Arrays, self.Bulk_Pattern_Rows, axis=1)
        
        #MAKING (5,5,5),12 A BOOLEAN PATTERN, AKA IF ITS FOUND OR NOT ONLY
        self.Bulk_Array[np.nonzero(self.Bulk_Array[:,0]),0] = 1
        self.X=self.Bulk_Array.astype(dtype='double', order='C')
        
        
        #DIVIDING THE SECOND TO FOURTH PARAMETERS BY THE NUMBER OF ATOMS
        for i in range(3):
            self.X[:,i+1] /= np.asarray(self.Known_NAtom)
        
        self.scaling = np.maximum(np.max(self.X, axis=0))
        
        self.X /= self.scaling
        self.y = self.Structure_List
        #SCALING DOWN THE (5,5,5),12 PATTERN VALUES, AS THEIR RANGE AFFECTS TOO MUCH
        #THE Ih JUDGEMENT OTHERWISE
        self.X[:,0] /= 10
        
        
        with open('Training_Info.txt', 'a')as f:
            f.write('Coordination Number of each file:\n')
            f.write(Coordination_Number)
            
            f.write('\nNumber of Unrecognized Atoms in each file:\n')
            f.writelines(Known_NAtom-np.sum(Known_Pattern_Arrays,axis=1))
        
            f.write('\nProcessed Known Data:\n')
            for i in range(len(self.X)):
                f.write('%s \t %s \n' %(self.X[i,:], self.y[i]))
            f.write('Processed Unknown Data:\n')
            f.write(self.prediction)
            f.close()
        
        self.clf = svm.SVC()
        self.clf.fit(self.X, self.y)
        
        os.chdir(self.script_path)
        os.chdir('../')
        joblib.dump(self.clf, ModelName)
        os.chdir(self.cwd)
        
        with open('Training_Info.txt', 'a')as f:
            f.write('Support vectors:\n')
            f.write(clf.support_vectors_)
            f.write('\nIndices of support vectors:\n')
            f.write(clf.support_)
            f.write('\nNumber of support vectors for each class:\n')
            f.write(clf.n_support_)