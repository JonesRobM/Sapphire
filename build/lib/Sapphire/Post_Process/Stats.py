import numpy as np
from scipy.stats import ks_2samp
import ruptures as rpt

def KB_Dist(P,Q):
    
    """ Robert
    Calculates the Kullback-Liebler divergence between two distributions.
    
    P: The "initial" distribution against which one wishes to measure the mutual
    entropy of the distribution
    
    Q:
    
    At the moment, there is no actual provision to protect against zero division errors.
    One possible solution could be to define a local varaible, epsilon, which is added to 
    every point in P and prevents it from being zero at any point. 
    
    Note that these two distributions must have identical dimensions or the script
    will not run. 
    
    A reasonable work-around is to define both from an identical linspace.
    """
    
    
    return np.sum(np.where(P != 0, P * np.log(P / Q), 0))

def JSD_Dist(P,Q):

    K=0
    Epsilon=0.000001
    Q+=Epsilon
    P+=Epsilon
    for x in range(len(Q)):
        K-=0.5*(P[x]*np.log(2*Q[x]/(Q[x]+P[x])) + Q[x]*np.log(2*P[x]/(P[x]+Q[x])))
    return np.sqrt(K)



class Dist_Stats():
    
    """ Jones
    
    This class of functions is a group of statistical techniques
    that I am experimenting with as a means of identifying "significant"
    changes in distributions of random variables such as:
        
        Radial Distribution Function (RDF)
        
        Pair Distance Distribution Function (PDDF)
        
        NOT PDF as that means Probability Distribution Function 
        {Not to be confused with the Probability Mass Function}
        
        CNA signature distribution.
        
    Note that these tools do not a-priori require you to have normalised distributions.
    Where it is necessary that they are (PDDF, CNA Sigs), the functional form written
    ensures that they already are.
    
    Isn't life nice like that? :D
    
    In this realisation of the code, each analysis code is to be called for each time frame.
    See the example script.
    
    """
            
    def __init__(self, PStats, KL, JSD):
        self.PStats = PStats
        self.KL = KL
        self.JSD = JSD
        
        
    def PStat(Ref_Dist, Test_Dist):
        
        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            PStats: The Kolmogorov Smirnov Test statistic.
            Generally speaking, this being <0.05 is sufficing grounds to
            reject the null hypothesis that two sets of observations are drawn
            from the same distribution
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
            
            
        PStats = (ks_2samp(Ref_Dist,Test_Dist)[1]) #Performs KS testing and returns p statistic
        #if frame+1 == len(Dist):
           #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return PStats
    
    def Kullback(Ref_Dist, Test_Dist):

        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            KL: The Kullback Liebler divergence:
                This is also known as the mutual information between two distributions.
                It may loosely (and dangerously) interpreted as the similarity between
                two distributions. 
                
                I care about plotting this as I suspect strong delineations in the growth
                of mutual entropy as the system undergoes a phase transition.
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
            
        KL = (KB_Dist(Ref_Dist,Test_Dist))
        #if frame+1 == len(Dist):
            #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return KL
        
    def JSD(Ref_Dist, Test_Dist):
        
        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            J: Jenson-Shannon Distance which is a symmetric form the the KL distance above.
            I do not yet understand fully why this should be a superior function to KL but 
            it's another telling discriptor.
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
        
        J = (JSD_Dist(Ref_Dist,Test_Dist))
        #if frame+1 == len(Dist):
            #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return J

def autocorr1(x,lags):
    '''numpy.corrcoef, partial'''

    corr=[1. if l==0 else numpy.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    return numpy.array(corr)

def autocorr2(x,lags):
    '''manualy compute, non partial'''

    mean=numpy.mean(x)
    var=numpy.var(x)
    xp=x-mean
    corr=[1. if l==0 else numpy.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]

    return numpy.array(corr)

def autocorr3(x,lags):
    '''fft, pad 0s, non partial'''

    n=len(x)
    # pad 0s to 2n-1
    ext_size=2*n-1
    # nearest power of 2
    fsize=2**numpy.ceil(numpy.log2(ext_size)).astype('int')

    xp=x-numpy.mean(x)
    var=numpy.var(x)

    # do fft and ifft
    cf=numpy.fft.fft(xp,fsize)
    sf=cf.conjugate()*cf
    corr=numpy.fft.ifft(sf).real
    corr=corr/var/n

    return corr[:len(lags)]

def autocorr4(x,lags):
    '''fft, don't pad 0s, non partial'''
    mean=x.mean()
    var=numpy.var(x)
    xp=x-mean

    cf=numpy.fft.fft(xp)
    sf=cf.conjugate()*cf
    corr=numpy.fft.ifft(sf).real/var/len(x)

    return corr[:len(lags)]

def autocorr5(x,lags):
    '''numpy.correlate, non partial'''
    mean=x.mean()
    var=numpy.var(x)
    xp=x-mean
    corr=numpy.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)

    return corr[:len(lags)]

def Changes(Data, model = 'rbf', lag = 10):
    algo = rpt.Pelt(model=model).fit(Data)
    result = algo.predict(pen=10)
    return result

def R(AdjT, AdjDeltaT):
    TempT = spa.csr_matrix.todense(AdjT)
    TempDeltaT = spa.csr_matrix.todense(AdjDeltaT)
    Temp = TempT-TempDeltaT
    return [ bool(x) for x in Temp.sum(1) ]

def Collectivity(R):
    return float(sum(R)/len(R))

def Concertedness(H1, H2):
    return abs(H2-H1)
