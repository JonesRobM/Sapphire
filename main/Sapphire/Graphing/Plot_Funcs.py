import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
import multiprocessing as mp
from functools import partial
import os
import linecache
import sys
import traceback 
from inspect import getmembers, isfunction
import inspect

plt.rcParams.update({'figure.max_open_warning': 0})

def distance(a, b): 
    
    dx = abs(a[0] - b[0])
     
    dy = abs(a[1] - b[1])
     
    dz = abs(a[2] - b[2])
 
    return np.sqrt(dx**2 + dy**2 + dz**2)

def Gauss(Data,Band, mon=False, Space = None):

    if Space is None:
        Space = np.linspace(2, 8, 100)
        Data = [elem for elem in Data if 1 < elem < 9] 
    A=[]; Data=np.asarray(Data) 
    
    
    if len(Data) > 10:
    
        for i in range(len(Data)):
            A.append(norm.pdf(Space, Data[i],Band))
        Density = np.asarray(np.sum(A, axis=0))
        Density = Density/np.trapz(Density, Space) #For normalisation purposes
        if mon == False:
            Min = (np.diff(np.sign(np.diff(Density))) > 0).nonzero()[0] + 1 # local min
            R_Cut = Space[Min][np.where(Space[Min]>3)][0]
            return Space, Density, R_Cut
        elif mon == True:
            return Space, Density
    else:
        return None

class Plot_Funcs():
    
    def __init__(self, MetaData = None, Errors = None, Quantities = None, System = None):
        
        if System == None:
             self.System = None
             self.Base = ''
             self.Images = ''
             self.single_file = True
             
        else:
            self.System = System
            try:
                self.Base = System['base_dir']
            except KeyError:
                self.Base = ''
                
            try:
                self.Images = System['plot_dir']
                self.ensure_dir(self.Base + self.Images)
            except KeyError:
                self.Images = ''
                
        if MetaData is None:
            sys.exit("\nNo metadata provided for analysis.\nNow exiting.\n")
        else:
            self.Meta = MetaData
            
        if Errors is None:
            self.Errors = False
            with open(self.Base+'Plotting_Info.txt', "a+") as f:
                f.write("\nNo errors have been provided.\nHence, no errors will be plotted.\n")
        else:
            self.Err = Errors
            self.Errors = True
        
        if Quantities is None:
            sys.exit("\nNo quantities requested.\nNow exiting.\n")
        else:
            self.Quantities = Quantities
        

            
        self.functions_list = [o for o in getmembers(Plot_Funcs) if isfunction(o[1])]
        self.Functions = {}
        
        """This provides a dictionary with the function names as keys and the 
        function itself, plus arguments following.
        The reason for the arguments is so that user defined input arguments 
        may be identified and fed in correctly."""
        
        for x in self.functions_list:
            self.Functions[x[0]] = inspect.getfullargspec(x[1])[0][1:]
        self.Q_Keys = self.Quantities.keys()
        self.Meta_Keys = self.Meta.keys()
        self.Plot_Dict = {}
        for obj in self.Q_Keys:
            for item in self.functions_list:
                if obj.lower() in item[0].lower():
                    self.Plot_Dict[item[0]] = [item[1]]
                    
    def ensure_dir(self, file_path=''):
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
    def Make_Plots(self):
        """
        
        Robert:
            This is the function that calls all of the desired functions for
            creating plots.
            The list of function names and arguments are already pre-defined and
            so this function simply parses through the user input.
            
            Still need to make a robust sanitation of user input but that may come later.
                  

        Returns
        -------
        None.

        """
        
        for x in self.Q_Keys:
            if x in self.Functions:
                ArgsList = []
                for y in self.Functions[x]:
                    try:
                        ArgsList.append(self.Quantities[x][y])
                    except KeyError:
                        ArgsList.append(
                            inspect.getargspec(self.Plot_Dict[x][0])[-1][self.Functions[x].index(y)]
                            )
                    
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\nThe arguments for function %s are %s.\n"%(x,ArgsList))
                getattr(self, x)(*ArgsList)
                
                    
        
    def Collect_CNA(self, Sig):
        try:
            Index = self.Meta["masterkey"].index( Sig )
        
            return [ self.Meta['cna_sigs'][x][Index] for x in range(len(self.Meta['cna_sigs'])) ]
        except KeyError:
            with open(self.Base+'Plotting_Info.txt', "a") as f:
                f.write("\nData not found in metadata\n")
            return None

    def Collect_CNA_error(self, Sig):
        try:
            Index = self.Meta["masterkey"].index( Sig )
        
            return [ self.Err['cna_sigs'][x][Index] for x in range(len(self.Err['cna_sigs'])) ]
        except KeyError:
            with open(self.Base+'Plotting_Info.txt', "a") as f:
                f.write("\nData not found in metadata\n")
            return None
        
    def autolabel(self, rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize = 18)
            
    def agcn_heat(self, Name = 'agcn_Heat.png'):
        Bins = np.linspace(3,12,41)
        Heat = []
        try:
            for frame in range( len(self.Meta['agcn']) ):
                a,b = np.histogram( self.Meta['agcn'][frame], bins = Bins )
                Heat.append(a)
                
            YTicks = np.array( [ "{:.1f}".format(x) for x in np.linspace(3,12,20) ] )
            try:
                XTicks = np.array( [ "{:.0f}".format(t) for t in np.linspace( self.Meta['SimTime'][0], self.Meta['SimTime'][-1] ,25) ], dtype = int )
            except KeyError:
                XTicks = np.array( [ "{:.0f}".format(t) for t in np.linspace( self.Meta['Start'], self.Meta['End'] ,25) ], dtype = int )
            Heat = ( np.asanyarray(Heat) ).transpose()
            
            
            ax = sns.heatmap(Heat, cmap = 'hot')
            ax.set_xlabel("Frame", fontsize = 14)
            ax.set_ylabel("AGCN", fontsize =14)
            ax.set_xticklabels(XTicks)
            ax.set_yticklabels(YTicks)
            plt.savefig(self.Base+self.Images+'/'+Name, dpi = 100, bbox_inches='tight')
            plt.close()
        except KeyError:
            print("\nThis quantity does not exist in the metadata.\n")
            return None           
    

            
    def prdf_plot(self, Names = None, Frames = [], He = False, Ho = None, Errors = False):
        
        Frames = list(Frames)
        if self.Errors is True:
            Errors = True
        
        
        """
        Name: str 'pdf' 'rdf'
        
        Frames: list frames to be reviewed
        He: bool Whether to look for hetero quantities default is False
        Homo: List of atomic species to be considered as homo pairs only - default is empty list
        

        Parameters
        ----------
        Name : TYPE
            DESCRIPTION.
        Frames : TYPE, optional
            DESCRIPTION. The default is [].
        He : TYPE, optional
            DESCRIPTION. The default is None.
        Homo : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        
        for Name in Names:
        
            for frame in Frames:
            
                fig, ax = plt.subplots()
                fig.set_size_inches((9,3))
                try:
                    ax.plot(self.Meta[Name][Frames.index(frame)][0], self.Meta[Name][Frames.index(frame)][1], 
                            color='k', linestyle = 'solid', linewidth = 4, label = "Full system")
                
                    if Errors is True:
                        ax.fill_between(self.Meta[Name][Frames.index(frame)][0], 
                                        self.Meta[Name][Frames.index(frame)][1] + self.Err[Name][Frames.index(frame)][1],
                                        self.Meta[Name][Frames.index(frame)][1] - self.Err[Name][Frames.index(frame)][1], 
                                        color='k', alpha = 0.25)
                    fig.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                           ncol=2, mode="expand", borderaxespad=0. ,fontsize = 14)
                except KeyError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\n%s was not found in the provided metadata.\n"%Name)
                    
                if He is False:
                    pass
                else:
                    try:
                        ax.plot(self.Meta['He'+Name.upper()][Frames.index(frame)][0], self.Meta['He'+Name.upper()][Frames.index(frame)][1], 
                                color='r', linestyle = 'dashed', linewidth = 4, label = "Pair different only")
                        if Errors is True:
                            ax.fill_between(self.Meta['He'+Name.upper()][Frames.index(frame)][0], 
                                            self.Meta['He'+Name.upper()][Frames.index(frame)][1] + self.Err['He'+Name.upper()][Frames.index(frame)][1],
                                            self.Meta['He'+Name.upper()][Frames.index(frame)][1] - self.Err['He'+Name.upper()][Frames.index(frame)][1], 
                                            color='r', alpha = 0.25)
                    except KeyError:
                        with open(self.Base+'Plotting_Info.txt', "a") as f:
                            f.write("\n%s was not found in the metadata.\n"%Name)
                
                if Ho is None:
                    pass
                elif type(Ho) is list:
                    for ele in Ho:
                        ax.plot(self.Meta['Ho'+Name.upper()+ele][Frames.index(frame)][0], 
                                self.Meta['Ho'+Name.upper()+ele][Frames.index(frame)][1],
                                linestyle = 'dashdot', linewidth = 4, 
                                label ="%s only"%ele)
                        if Errors is True:
                            ax.fill_between(self.Meta['Ho'+Name.upper()+ele][Frames.index(frame)][0], 
                                            self.Meta['Ho'+Name.upper()+ele][Frames.index(frame)][1] + self.Err['Ho'+Name.upper()+ele][Frames.index(frame)][1],
                                            self.Meta['Ho'+Name.upper()+ele][Frames.index(frame)][1] - self.Err['Ho'+Name.upper()+ele][Frames.index(frame)][1], 
                                            alpha = 0.25)
                        
                else:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nError in input arguments for the prdf_plot.\n")
                
                ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
                ax.set_xlabel(r"Distance (Angstrom)", fontsize = 12)
                ax.set_ylabel(Name, fontsize = 12)
                
                fig.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                           ncol=2, mode="expand", borderaxespad=0. ,fontsize = 14)
                
                try:
                    FrameTemp = "{:.1f}".format(self.Meta['Temp'][int(frame)])
                    ax.text(1.4*np.amin(self.Meta[Name][Frames.index(frame)][0]), 0.9*np.amax(self.Meta[Name][Frames.index(frame)][1]), 
                            "Time: %sps\nTemp: %sK" %(self.Meta['SimTime'][int(frame)], 
                                                      FrameTemp), fontsize=13)
                except KeyError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\n%s threw an error when using the prdf_plot function.\n"%frame)
                
                if He is False:
                    if Ho is None:
                        plt.savefig(self.Base + self.Images + '/' + Name.upper() + str(frame)+'.png', 
                                    dpi = 100, bbox_inches='tight')
                    elif type(Ho) is list:
                        plt.savefig(self.Base + self.Images + '/' + Name.upper() + '_Ho_' + ''.join(map(str, Ho)) +'_' + str(frame)+'.png', 
                                    dpi = 100, bbox_inches='tight')
                    else:
                        plt.savefig(self.Base + self.Images + '/' + Name.upper() + str(frame)+'.png', 
                                    dpi = 100, bbox_inches='tight')
                else:
                    if Ho is None:
                        plt.savefig(self.Base + self.Images + '/' + Name.upper() +'_He_' + str(frame)+'.png', 
                                    dpi = 100, bbox_inches='tight')
                    elif type(Ho) is list:
                        plt.savefig(self.Base + self.Images + '/' + Name.upper() +'_He_' + '_Ho_' + ''.join(map(str, Ho)) +'_' + str(frame)+'.png', 
                                    dpi = 100, bbox_inches='tight')
                    else:
                        plt.savefig(self.Base + self.Images + '/' + Name.upper() +'_He_' + str(frame)+'.png', 
                                    dpi = 100, bbox_inches='tight')
                plt.close()
                        
    
    def plot_stats(self, Stats = [], Species = None, Quants = [], Temp = False, Errors = False, Frames = None):

        if self.Errors is True:
            Errors =  True

        if Frames is None:       
            try:
                TimeAxis = range(int(self.Meta['Start']), 
                                 int(self.Meta['SimTime'][-1]), 
                                 int(int(self.Meta['Skip']) * int(self.Meta['SimTime'][-1]) / int(self.Meta['End'])))
            except KeyError:
                TimeAxis = range(int(self.Meta['Start']), 
                                 int(self.Meta['End']), 
                                 int(self.Meta['Step']))
        else:
            TimeAxis = Frames

        for Stat in Stats:
            fig,ax = plt.subplots()
            fig.set_size_inches((9,3))
            for Quant in Quants:
                try:                
                    ax.plot(TimeAxis, 
                            self.Meta[Stat+Quant.lower()], 
                            label = Quant.lower())
                    if Errors is True:
                        ax.fill_between(TimeAxis, 
                                        self.Meta[Stat+Quant.lower()] - self.Err[Stat+Quant.lower()],
                                        self.Meta[Stat+Quant.lower()] + self.Err[Stat+Quant.lower()],
                                        alpha = 0.25)
                except KeyError:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nNo %s found in metadata.\n"%(Stat+Quant.lower()))
                
            try:
                ax2=ax.twinx()
                
                ax2.scatter(TimeAxis,
                            self.Meta['R_Cut'], 
                            linewidths=4, label = 'R_Cut', color='g')
                if Species is not None:
                    for x in Species:
                        ax2.scatter(TimeAxis,
                                    self.Meta['Cut' + x], 
                                    linewidths=4, label = 'R_Cut' + x)
                
                if Errors is True:
                    ax2.errorbar(TimeAxis, self.Meta['R_Cut'],
                                 self.Err['R_Cut'], color='g',
                                 capsize = 5, capthick = 3)
                    if Species is not None:
                        for x in Species:
                            ax2.errorbar(TimeAxis,self.Meta['Cut' + x], 
                                         self.Err['Cut' + x],
                                         capsize = 5, capthick = 3)
            except KeyError:
                pass
                        
            fig.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                       ncol=3, mode="expand", borderaxespad=0. ,fontsize = 12)
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel(Stat.upper())
            
            if Temp is True:
                ax3 = ax.twiny()
                
                ax1Ticks = ax.get_xticks()
                ax3Ticks = ax1Ticks
                
                ax3.set_xticks(ax2Ticks)
                ax3.set_xbound(ax.get_xbound())
                ax3.set_xticklabels(tick_function(ax2Ticks))
                ax3.set_xlabel('Temperature (K)')
    
            plt.savefig(self.Base + self.Images + '/' + str(Stat) + '.png' , dpi = 100, bbox_inches='tight')
            plt.close()
    
    
    
    def tick_function(self, X):
        try:
           inc = (max(self.Meta['Temp']) - min(self.Meta['Temp']))/( 10*len(self.Meta['Temp']) )
           V = min(self.Meta['Temp']) + X*inc
           return ["%.3f" % z for z in V]
        except KeyError:
            return None
         
    def com_plot_bi(self, Dists = None, Species = None, Frames = [0], Errors =  False):
        if self.Errors is True:
            Errors = True
        
        if Dists is None:
            with open(self.Base+'Plotting_Info.txt', "a") as f:
                f.write("\nNo distributions requested.\n")
            return None
        elif type(Dists) is list:
            for Dist in Dists:
                                        
                if Dist is "MidCoMDist":
                    D = "Cluster Centre"
                elif Dist is "CoMDist":
                    D = "Sub-cluster Centre"
                else:
                    raise KeyError("Invalid distribution.\n")
                
                if Species is None:
                    with open(self.Base+'Plotting_Info.txt', "a") as f:
                        f.write("\nNo chemical species requested.\n")
                elif type(Species) is list:
                
                    for Specie in Species:
                    
                        for frame in Frames:
                                                  
                            try:
                                fig,ax = plt.subplots()
                                fig.set_size_inches(9,3)
                                ax.plot(self.Meta['CoMSpace'], self.Meta[Dist + Specie][frame], color= 'k', linewidth = 4)
                                if Errors is True:
                                    ax.fill_between(self.Meta['CoMSpace'], 
                                                    self.Meta[Dist + Specie][frame] + self.Err[Dist + Specie][frame], 
                                                    self.Meta[Dist + Specie][frame] - self.Err[Dist + Specie][frame],
                                                    color = 'k', alpha = 0.25)
                                
                                ax.set_xlabel('Distance (Angstrom)')
                                
                                ax.set_ylabel('Probability')
                                
                                try:
                                    ax.text(self.Meta['CoMSpace'][5], 0.65*max(self.Meta[Dist + Specie][frame]), "%s to %s\nTime: %sps\nTemp: %sK" 
                                            %(Specie, D, self.Meta['SimTime'][frame], "{:.1f}".format(self.Meta['Temp'][frame])))
                                except KeyError:
                                    pass
                                
                                plt.savefig(self.Base + self.Images + '/' + Dist+Specie+str(frame) + '.png', 
                                            dpi = 100, bbox_inches='tight')
                                plt.close()
                                
                            except KeyError:
                                with open(self.Base+'Plotting_Info.txt', "a") as f:
                                    f.write("\nThere was an error trying to plot %s.\n" %(Dist+Specie))
                                pass
                            
    def cna_plot(self, Name = 'CNA_Time', Frames = [], Errors = False):
        if self.Errors is True:
            Errors = True
                
        for Frame in Frames:
            try:
                X_CNA = [ str(a) for a in self.Meta['masterkey'] ] # Create a set of ticks for the x-axis
                fig = plt.figure(figsize = (9,3) )
                if Errors is True:
                    ax = plt.bar( X_CNA, self.Meta['cna_sigs'][Frame], yerr = self.Err['cna_sigs'][Frame], tick_label = X_CNA )
                else:
                    ax = plt.bar( X_CNA, self.Meta['cna_sigs'][Frame], tick_label = X_CNA)
                plt.xlabel("CNA Signature", fontsize = 14)
                plt.ylabel("Probability", fontsize = 14)
                plt.xticks(rotation=90,fontsize = 14)
                try:
                    plt.text( X_CNA[-7], 0.8*np.amax(self.Meta['cna_sigs'][Frame]), 
                            'Time: %sps\nTemp: %sK' %(self.Meta["SimTime"][Frame], 
                                                      "{:.1f}".format(self.Meta['Temp'][Frame])), fontsize = 14 )
                except KeyError:
                    pass
                plt.savefig(self.Base+self.Images+'/'+Name+str(Frame)+'.png', dpi = 100, bbox_inches = 'tight')
                plt.close()
            except KeyError:
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write("\nThis quantitiy, cna, does not exist in the metadata.\n")
                return None        
        
    
    def agcn_histo(self, Frames = [], Errors = False):

        for Frame in Frames:
            fig, ax = plt.subplots()
            fig.set_size_inches(9,3)
            y,binEdges = np.histogram(self.Meta['agcn'][Frame], bins = 40)
            bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
            ax.bar(bincenters, y, color='r')
            try:
                ax.text(bincenters[4], 0.7*np.amax(y), "Time : %sps\nTemp : %sK"%(self.Meta['SimTime'][Frame], "{:.1f}".format(self.Meta['Temp'][Frame])) )
                plt.savefig(self.Base + self.Images + '/'+ 'AGCNDist'+str(self.Meta['SimTime'][Frame])+'.png', dpi = 100, bbox_inches='tight')
            except KeyError:
                plt.savefig(self.Base + self.Images + '/'+ 'AGCNDist.png', dpi = 100, bbox_inches='tight')
            plt.close()
    
    def com_full_plot(self, Frames = [], Errors =  False):
        if self.Errors is True:
            Errors = True
        
        for Frame in Frames:
            fig, ax = plt.subplots()
            fig.set_size_inches(9,3)
            ax.plot(self.Meta['CoMSpace'], self.Meta['CoMDist'][Frame], color='k')
            if Errors is True:
                ax.fill_between(self.Meta['CoMSpace'] ,
                                self.Meta['CoMDist'][Frame] + self.Err['CoMDist'][Frame],
                                self.Meta['CoMDist'][Frame] - self.Err['CoMDist'][Frame],
                                color='k', alpha = 0.25)
            ax.set_xlabel('Distance (Angstrom)')
            ax.set_ylabel('RDF')
            try:
                ax.text(self.Meta['CoMSpace'][5], 0.65*max(self.Meta['CoMDist'][Frame]), "Full System\nTime: %sps\nTemp: %sK" %(self.Meta['SimTime'][Frame], "{:.1f}".format(self.Meta['Temp'][Frame])))
                plt.savefig(self.Base + self.Images + '/'+ 'FullCoM'+str(self.Meta['SimTime'][Frame])+'.png', 
                            dpi = 100, bbox_inches='tight')
            except KeyError:
                plt.savefig(self.Base + self.Images + '/'+ 'FullCoM.png', dpi = 100, bbox_inches='tight')
            plt.close()
        
        
    def Mass(self,r):
        return ((4/3) * np.pi * r**3 )
    
    def cum_com(self, Frames):
        fig,ax = plt.subplots()
        fig.set_size_inches(9,3)
        for Frame in Frames:
            Int = [ np.trapz(self.Meta['CoMDist'][Frame][:x], self.Meta['CoMSpace'][:x]) for x in range(100) ]
            try:
                ax.plot(self.Meta['CoMSpace'], Int, label = '%sps' %(self.Meta['SimTime'][Frame]))
            except KeyError:
                ax.plot(self.Meta['CoMSpace'], Int, label = str(Frame))
        ax.plot(self.Meta['CoMSpace'], self.Mass(self.Meta['CoMSpace'])/max(self.Mass(self.Meta['CoMSpace'])), label = 'Spherical mass distribution', linestyle = 'dashed')
        ax.set_xlabel('Distance from centre (Angstrom)')
        ax.set_ylabel('M(r) / M(R)')
        fig.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                   ncol=3, mode="expand", borderaxespad=0. ,fontsize = 12)
        try:
            plt.savefig(self.Base + self.Images + '/'+ 'Cum_CoM'+str(self.Meta['SimTime'][Frame])+'.png', 
                            dpi = 100, bbox_inches='tight')
        except KeyError:
            plt.savefig(self.Base + self.Images + '/'+ 'Cum_CoM.png', 
                            dpi = 100, bbox_inches='tight')            
        plt.close()
        
    def cna_traj(self, Sigs = [], Errors = False):
        if self.Errors is True:
            Errors = True
        try:
            Time = self.Meta['SimTime']
        except KeyError:
            Time = range(len(self.Meta['cna_sigs']))
        fig,ax = plt.subplots()
        fig.set_size_inches(9,3)
        for x in Sigs:
            try:
                ax.plot(Time, self.Collect_CNA(x), label = x)
                if Errors is True:
                    ax.fill_between(Time, 
                                    np.asarray(self.Collect_CNA(x)) + np.asarray(self.Collect_CNA_error(x)),
                                    np.asarray(self.Collect_CNA(x)) - np.asarray(self.Collect_CNA_error(x)),
                                    alpha = 0.25)
            except ValueError:
                print(x, type(x))
                with open(self.Base+'Plotting_Info.txt', "a") as f:
                    f.write(f"\nSignature, '{0}', not in metadata.\n".format(x))
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Probability')
        fig.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                   ncol=3, mode="expand", borderaxespad=0. ,fontsize = 12)   
        plt.savefig(self.Base + self.Images + '/'+ 'CNA_Traj'+'.png', 
                        dpi = 100, bbox_inches='tight')
        plt.close()
       
        
    def h_c(self, Errors = False):
        if self.Errors is True:
            Errors = True
            
        Time = self.Meta['SimTime']
        fig,ax = plt.subplots()
        fig.set_size_inches(9,3)

        ax.plot(Time, self.Meta['h'], label = 'Collectivity')
        ax.plot(Time, self.Meta['c'], label = 'Concertedness')
        if Errors is True:
            ax.fill_between(Time[1:], 
                            self.Meta['h']+self.Err['h'],
                            self.Meta['h']-self.Err['h'],
                            alpha = 0.25)
            ax.fill_between(Time[2:-1], 
                            self.Meta['c']+self.Err['c'],
                            self.Meta['c']-self.Err['c'],
                            alpha = 0.25)
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel(' H / C')
        fig.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                   ncol=3, mode="expand", borderaxespad=0. ,fontsize = 12)   
        plt.savefig(self.Base + self.Images + '/'+ 'HC_Stats'+'.png', 
                        dpi = 100, bbox_inches='tight')
        plt.close()
     
        

def pair_plot(Data, System):
    try:
        HeAdj = Data['HeAdj']
        NewHe = []
    except KeyError:
        sys.exit()
    for x in range(len(HeAdj)):
        try:
            NewHe.append(sum(HeAdj[x][1]))
        except TypeError:
            pass

    fig,ax = plt.subplots()
    fig.set_size_inches(9,3)
    ax.plot(Data['SimTime'], [sum(Data['HoAdjPt'][x]) for x in range(len(Data['HoAdjPt']))], 'orange', label='Pt only')
    ax2 = ax.twinx()
    ax2.plot(Data['SimTime'], [sum(Data['HoAdjAu'][x]) for x in range(len(Data['HoAdjAu']))] , 'blue', label = 'Au only')
    ax3 = ax.twiny()
    ax3.plot(NewHe, label = 'Hetero pairs only', color='red')
    ax2.axes.yaxis.set_visible(False)
    ax3.axes.xaxis.set_visible(False)
    labels = [item.get_text() for item in ax.get_yticklabels()]
    
    empty_string_labels = ['']*len(labels)
    ax.set_yticklabels(empty_string_labels)
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Number of pairs')
    fig.legend(bbox_to_anchor=(0, 1.0, 1., 0), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0. ,fontsize = 12)
    plt.savefig(System['base_dir']+System['plot_dir'] + '/Pairs.png', dpi = 100, bbox_inches='tight')      
        
    
def All_CNA_Traj(System, Pipeline, outfile):
    CNA = []
    for x in System['iter_dir']:
        for y in [(4,2,2), (4,2,1), (3,1,1)]:
            Index = Pipeline.BigMeta[x]['cna'][0][0].index(y)
            Temp = [ Pipeline.BigMeta[x]['cna'][i][1][Index] for i in range(len(Pipeline.BigMeta[x]['cna'])) ]
            CNA.append(Temp)
    
    x = Pipeline.BigMeta[System['iter_dir'][0]]['SimTime']
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
    fig.set_size_inches(9,3)
    (ax1, ax2), (ax3, ax4) = axs
    
    ax1.plot(x, CNA[0], label = '(4 2 2)')
    ax1.plot(x, CNA[1], label = '(4 2 1)')
    ax1.plot(x, CNA[2], label = '(3 1 1)')
    
    ax2.plot(x, CNA[3])
    ax2.plot(x, CNA[4])
    ax2.plot(x, CNA[5])
    
    ax3.plot(x, CNA[6])
    ax3.plot(x, CNA[7])
    ax3.plot(x, CNA[8])
    
    ax4.plot(x, CNA[9])
    ax4.plot(x, CNA[10])
    ax4.plot(x, CNA[11])
    
    for ax in axs.flat:
        ax.label_outer()
        
        ax.set_ylim(0, 0.7)
    fig.legend( loc='upper center', ncol=3, fontsize = 10)
    plt.savefig(outfile, dpi = 100, bbox_inches='tight') 
        
    
        
    """
      
    ##########################################################################
    
     The following are old functions with little utility but may be
     reintroduced if there is demand for such things.
     
    def AGCN_Excess():
        Excess = []
        for i in range( len( AverageMeta['agcn'] ) ):
            Temp = [ a>12 for a in AverageMeta['agcn'][i] ]
            Excess.append(np.sum(Temp))
        return Excess
    
    def Strange_CNA():
        Indices = [ 14, 15, 24, 25, 38 ] #37 and on to the end are all odd
        CNA = AverageMeta['cna'] # All of the heights
        Strange_Dict = {}
        for Index in Indices:
            Strange_Dict[AverageMeta['masterkey'][Index]] = np.zeros((len(CNA)), dtype = np.float64)
            for Key in AverageMeta['masterkey'][Indices[-1]:]:
                Strange_Dict[Key] = np.zeros((len(CNA)), dtype = np.float64)
        Key = list(Strange_Dict.keys())
        Mast = AverageMeta['masterkey']
        for frame in range(len(CNA)):
            for Sig in CNA[frame]:
                for obj in Key:
                    if list(CNA[frame]).index(Sig) == Mast.index(obj):
                        if Sig > 0:
                            Strange_Dict[obj][frame] = 1
        Bar_Heights = []                    
        for Item in Strange_Dict:
            Bar_Heights.append( np.sum(Strange_Dict[Item]) )
        return (Strange_Dict.keys(), Bar_Heights)
      
  
    fig, ax = plt.subplots()
    fig.set_size_inches((21,7))
    ax.plot(New, label = '(4,5,5)', color='k')
    Ticks = range(0,1500,50)
    for tick in Ticks:
        ax.vlines(tick, ymin=0, ymax = 1.1*np.amax(New), color='r', linestyle = '--')
    ax2 = ax.twinx()
    ax2.scatter(Ticks, AverageMeta['R_Cut'], linewidths = 6, color='g')
    ax.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = 20)
    ax2.set_ylabel("Nearest neighbour cutoff (Angstrom)", fontsize = 20)
    ax.set_xlabel("Time (ps)", fontsize = 20)
    ax.set_ylabel("Probability", fontsize = 20)

    fig, ax = plt.subplots()
    fig.set_size_inches((21,7))
    rect =ax.bar(X_Key, B, tick_label = X_Key)
    plt.xticks(rotation = 90)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("CNA signature", fontsize = 20)
    ax.set_ylabel("Number of frames", fontsize=20)
    autolabel(rect)  
    
    
    def cna_plotter(Frame):
        X_CNA = [ str(a) for a in AverageMeta['masterkey'][:36] ] # Create a set of ticks for the x-axis
        fig = plt.figure(figsize = (9,3) )
        ax = plt.bar( X_CNA, AverageMeta['cna'][Frame][:36], tick_label = X_CNA )
      
        plt.xlabel("CNA Signature", fontsize = 12)
        plt.ylabel("Probability", fontsize = 12)
        plt.xticks(rotation=90,fontsize = 14)
        plt.text( X_CNA[20], 0.8*np.amax(AverageMeta['cna'][Frame]), 'Time: %sps\nTemp: %sK' %(AverageMeta["SimTime"][Frame], AverageMeta['Temp'][Frame]) )
        plt.savefig(path + 'Images/'+ 'CNA'+str(Frame)+'.png', dpi = 100, bbox_inches='tight')    


    ##########################################################################

    """