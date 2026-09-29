import numpy as np
from scipy.stats import ks_2samp


def _as_distribution(D):
    """A distribution sampled on a grid, as a probability vector (non-negative, sums to 1).

    Sapphire's PDFs are densities: they integrate to 1 over r, but their values sum to
    1/dr (~33 on the default 200-point grid). Divergences are defined between probability
    vectors, so each input is normalised to unit sum first; that also makes the result
    independent of the grid spacing.
    """
    D = np.clip(np.asarray(D, dtype=float), 0.0, None)
    total = D.sum()
    if not total > 0:
        raise ValueError("cannot compare an empty distribution")
    return D / total


class KB_Dist():

    def __init__(self,P,Q):
        
        
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
        
        self.P = P; self.Q = Q
        
        return None
        
    def calculate(self):
        
        P = _as_distribution(self.P); Q = _as_distribution(self.Q)
        mask = (P > 0) & (Q > 0)          # terms with P = 0 contribute 0; Q = 0 where P > 0 is +inf in theory
        with np.errstate(divide='ignore', invalid='ignore'):
            return float(np.sum(P[mask] * np.log(P[mask] / Q[mask])))

class JSD_Dist():
    
    def __init__(self, P, Q):
        
        self.P = P; self.Q = Q
        
        return None
        
    def calculate(self):
        """Jensen-Shannon distance, base 2: sqrt(JSD) with JSD = (KL(P||M) + KL(Q||M)) / 2.

        Bounded in [0, 1] (0 for identical, 1 for disjoint distributions) and symmetric.

        Until 1.3.x this computed -1/2 sum[P log(2Q/(P+Q)) + Q log(2P/(P+Q))] -- P and Q
        swapped inside the logarithms -- which equals half the symmetric (Jeffreys) KL minus
        the JSD: unbounded, not a Jensen-Shannon quantity, and applied to unnormalised
        densities. Values from earlier versions are not comparable with these.
        """
        P = _as_distribution(self.P); Q = _as_distribution(self.Q)
        M = 0.5 * (P + Q)

        def kl_to_m(A):
            nz = A > 0                    # 0 log 0 = 0; M > 0 wherever A > 0
            return np.sum(A[nz] * np.log2(A[nz] / M[nz]))

        jsd = 0.5 * kl_to_m(P) + 0.5 * kl_to_m(Q)
        return float(np.sqrt(max(jsd, 0.0)))  # clip rounding below zero for identical inputs



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
            
        KL = KB_Dist(Ref_Dist, Test_Dist).calculate()
        #if frame+1 == len(Dist):
            #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return KL
        
    def JSD(Ref_Dist, Test_Dist):
        
        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            J: Jensen-Shannon distance (base 2, in [0, 1]) between the two distributions after
            normalising each to unit sum. Unlike KL it is symmetric, always finite, and a
            true metric, so it is well suited to tracking a structure away from its start.
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
        
        J = JSD_Dist(Ref_Dist, Test_Dist).calculate()
        return J

class autocorr():
    
    def __init__(self):

        return None    
    
    def autocorr1(x,lags):
        '''np.corrcoef, partial'''
    
        corr=[1. if l==0 else np.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
        return np.array(corr)
    
    def autocorr2(x,lags):
        '''manualy compute, non partial'''
    
        mean=np.mean(x)
        var=np.var(x)
        xp=x-mean
        corr=[1. if l==0 else np.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]
    
        return np.array(corr)
    
    def autocorr3(x,lags):
        '''fft, pad 0s, non partial'''
    
        n=len(x)
        # pad 0s to 2n-1
        ext_size=2*n-1
        # nearest power of 2
        fsize=2**np.ceil(np.log2(ext_size)).astype('int')
    
        xp=x-np.mean(x)
        var=np.var(x)
    
        # do fft and ifft
        cf=np.fft.fft(xp,fsize)
        sf=cf.conjugate()*cf
        corr=np.fft.ifft(sf).real
        corr=corr/var/n
    
        return corr[:len(lags)]
    
    def autocorr4(x,lags):
        '''fft, don't pad 0s, non partial'''
        mean=x.mean()
        var=np.var(x)
        xp=x-mean
    
        cf=np.fft.fft(xp)
        sf=cf.conjugate()*cf
        corr=np.fft.ifft(sf).real/var/len(x)
    
        return corr[:len(lags)]
    
    def autocorr5(x,lags):
        '''np.correlate, non partial'''
        mean=x.mean()
        var=np.var(x)
        xp=x-mean
        corr=np.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)
    
        return corr[:len(lags)]

class ChangePoints():
    
    def __init__(self, Data, model = 'rbf', lag = 10):
        
        self.Data = Data
        self.model = model
        self.lag = lag
        
        return None

    def calculates(self):
        import ruptures as rpt  # optional: the 'changepoint' extra
        algo = rpt.Pelt(model=self.model).fit(self.Data)
        result = algo.predict(pen=self.lag)
        return result

class Mobility():
    
    def __init__(self, All_Adjacencies):
        self.Adj = All_Adjacencies
        return None

    def R(AdjT, AdjDeltaT):
        """Per-atom flag: did atom i's neighbour list change between the two frames?

        Accepts dense arrays or scipy sparse matrices. (The original summed the signed
        difference, so losing one neighbour and gaining another cancelled to 'unchanged'.)
        """
        a = np.asarray(AdjT.todense() if hasattr(AdjT, 'todense') else AdjT)
        b = np.asarray(AdjDeltaT.todense() if hasattr(AdjDeltaT, 'todense') else AdjDeltaT)
        return (a != b).any(axis=1)

    def Collectivity(R):
        """Fraction of atoms whose neighbourhood changed: 0 = frozen, 1 = every atom rearranged."""
        R = np.asarray(R)
        return float(R.sum() / len(R)) if len(R) else 0.0
    
    def Concertedness(H1, H2):
        return abs(H2-H1)


    def calculate(self):
        return None