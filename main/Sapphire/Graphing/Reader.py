import numpy as np
_trapz = getattr(np, 'trapezoid', None) or np.trapz  # noqa: NPY201
from scipy.stats import norm
import os

from Sapphire.Utilities.log import get_logger

log = get_logger('Sapphire.Graphing.Reader')

"""
Supported=[
        'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn',
        'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
           ]
"""

def get_heights_asap3(CNA, Masterkey, frame):
    CNA[frame] = list(CNA[frame])
    CNA[frame][1] = list(CNA[frame][1])
    def getkey(item):
        return item[0]
    
    FullSample=[]; Heights=[]
    Temp1=CNA[frame][0]
    for x in Masterkey:
        if x not in Temp1:    
            CNA[frame][0][x] = 0
    Sample=[]
    for j in Masterkey:
        Sample.append((j, CNA[frame][0][j]))
    FullSample.append(Sample)
    FullSample.sort()
    A,B=zip(*FullSample[0])
    Heights.append(B/np.sum(B))
    Temp = FullSample[0]
    FullCNA = [ item[0] for i, item in enumerate(Temp) ]
    
    Sample = (FullCNA, Heights[0])
    Altern = []
    for x in range(len(FullCNA)):
        Altern.append( ( FullCNA[x],Heights[0][x] ) )
    Sample = sorted(Altern, key = getkey)
    
    return ( [ Sample[x][0] for x in range(len(Sample)) ], [ Sample[x][1] for x in range(len(Sample)) ])

def Get_Heights_Ovito(CNAs, Masterkey, Norm = False):
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        CNAs: The dictionary containng the time ordered CNA signatures 
        and the number of observed occurances.
        
        
        MasterKey: The output from calling the Master function.
        This is to do pairwise comparrison for creating full 
        distributions without having to know what the craic is.
            
        Norm: Default - False
        Whether or not the user wishes to normalise the distribution of 
        CNA signatures for each frame in order to perform meaningful
        statistical analysis.
        
        
    Returns:
        
        Heights: np.array(Frames/Skip, len(MasterKey)) The array containing 
        the (if desired) normalised
        distribution of CNA signature occurances. 
        
    """
    

    Heights=np.zeros((len(CNAs),len(Masterkey)))
    

    for frame in range(len(CNAs)):
        Temp = CNAs[frame].keys()
        for x in Masterkey:
            if x not in Temp:
        
                CNAs[frame][x] = 0

            Heights[frame][Masterkey.index(x)] = CNAs[frame][x]
            
            if Norm == True:
                Heights[frame] = Heights[frame]/sum(Heights[frame])
            
    return Heights

def ensure_dir(file_path=''):
    directory = os.path.dirname(file_path+'/')
    if not os.path.exists(directory):
        os.makedirs(directory)
        log.info("Made a new directory.")


class Read_Meta:
    """Build the metadata layout :class:`Plot_Funcs` expects from one or more Sapphire run
    directories, via :class:`Sapphire.IO.Reader.Reader`.

    ``System`` keys (all optional except ``base_dir``):

    * ``base_dir``   – run directory (or parent of the ``iter_dir`` runs); figures go to
      ``base_dir/plot_dir``.
    * ``iter_dir``   – list of run sub-directories to average over (errors = std across runs).
    * ``plot_dir``   – output folder for images, default ``Images/``.
    * ``frame_dt``   – ps per output frame; gives ``SimTime``. Default: frame index.
    * ``temperature``– optional per-frame temperatures (K) for annotations.
    * ``com_band``   – KDE bandwidth for centre-of-mass distance distributions (default 0.3 Å).

    Legacy keys produced: ``pdf/rdf`` as ``[(space, heights), ...]`` per frame, ``HePDF``,
    ``HoPDF<X>``, ``R_Cut``, ``Cut<X>``, ``cna_sigs`` (normalised), ``masterkey`` (tuples),
    ``agcn``, ``CoMSpace``, ``CoMDist``, ``CoMDist<X>``, ``MidCoMDist<X>``, ``h``, ``c``,
    ``Start/End/Step/Skip``, ``SimTime``, ``Temp``.
    """

    def __init__(self, System=None):
        from Sapphire.IO.Reader import Reader
        self.System = dict(System or {})
        self.Base = self.System.get('base_dir', '')
        self.Images = os.path.join(self.Base, self.System.get('plot_dir', 'Images/'))
        self.Iter = self.System.get('iter_dir') or False
        os.makedirs(self.Images, exist_ok=True)
        runs = [os.path.join(self.Base, str(d), '') for d in self.Iter] if self.Iter else [self.Base]
        self.BigMeta = {str(run): self._adapt(Reader(run).load_all()) for run in runs}
        self.Keys = sorted(set().union(*(m.keys() for m in self.BigMeta.values())))
        self.AverageMeta, self.Errors = {}, {}

    # ------------------------------------------------------------------ adapter
    def _kde(self, samples, space):
        band = self.System.get('com_band', 0.3)  # Å; centre-of-mass distributions have ~N samples, not N^2
        out = np.zeros_like(space)
        s = np.asarray(samples, dtype=float)
        if s.size:
            out = np.sum(norm.pdf(space[None, :], s[:, None], band), axis=0)
            area = _trapz(out, space)
            if area > 0:
                out = out / area
        return out

    def _adapt(self, m):
        """Map Reader keys (see IO/OutputInfo*) onto the layout Plot_Funcs was written for."""
        out = {}
        species = sorted({k[len('hocut'):] for k in m if k.startswith('hocut')})
        n_frames = max((len(v) for k, v in m.items() if k in ('nn', 'pdf', 'rcut', 'cna_sigs') and hasattr(v, '__len__')), default=0)

        def dist(space_key, height_key, name):
            if space_key in m and height_key in m:
                sp, h = np.atleast_2d(m[space_key]), np.atleast_2d(m[height_key])
                out[name] = [(sp[min(i, len(sp) - 1)], h[i]) for i in range(len(h))]

        dist('pdfspace', 'pdf', 'pdf'); dist('rdfspace', 'rdf', 'rdf')
        dist('hepdfspace', 'hepdf', 'HePDF'); dist('herdfspace', 'herdf', 'HeRDF')
        for x in species:
            dist('hopdfspace' + x, 'hopdf' + x, 'HoPDF' + x); dist('hordfspace' + x, 'hordf' + x, 'HoRDF' + x)
            if 'hocut' + x in m:
                out['Cut' + x] = np.asarray(m['hocut' + x], dtype=float)
        if 'rcut' in m:
            out['R_Cut'] = np.asarray(m['rcut'], dtype=float)
        if 'agcn' in m:
            out['agcn'] = np.asarray(m['agcn'], dtype=float)
        if 'nn' in m:
            out['nn'] = np.asarray(m['nn'], dtype=float)
        if 'cna_sigs' in m:
            # Rows grow as new signatures join the masterkey: right-pad (absent == 0 counts).
            rows = [np.asarray(r, dtype=float) for r in m['cna_sigs']]
            width = max([len(r) for r in rows] + [len(m.get('masterkey', []))])
            sigs = np.array([np.pad(r, (0, width - len(r))) for r in rows])
            tot = sigs.sum(axis=1, keepdims=True); tot[tot == 0] = 1
            out['cna_sigs'] = sigs / tot
            out['masterkey'] = [tuple(int(c) for c in s) for s in m.get('masterkey', [])]
        # centre-of-mass distance distributions on a common grid
        com_keys = [k for k in ('comdist',) if k in m] + [k for k in m if k.startswith(('hocomdist', 'homidcomdist'))]
        if com_keys:
            rmax = max(float(np.max(np.asarray(m[k], dtype=float))) for k in com_keys)
            space = np.linspace(0, 1.1 * rmax, 100)
            out['CoMSpace'] = space
            for k in com_keys:
                name = {'comdist': 'CoMDist'}.get(k) or k.replace('hocomdist', 'CoMDist').replace('homidcomdist', 'MidCoMDist')
                out[name] = np.array([self._kde(frame, space) for frame in np.asarray(m[k], dtype=float)])
        if 'collect' in m:
            out['h'] = np.asarray(m['collect'], dtype=float)
        if 'concert' in m:
            out['c'] = np.asarray(m['concert'], dtype=float)
        for k in ('mix', 'gyration', 'stat_radius', 'surf_area'):
            if k in m:
                out[k] = np.asarray(m[k], dtype=float)
        for k, v in m.items():  # statistics and any other simple per-frame series
            if k not in out and isinstance(v, np.ndarray) and v.ndim == 1 and n_frames - 2 <= len(v) <= n_frames:
                out[k] = v.astype(float)
        # bookkeeping expected by the plotters
        out['Start'], out['End'], out['Step'], out['Skip'] = 0, n_frames, 1, 1
        dt = self.System.get('frame_dt')
        out['SimTime'] = np.arange(n_frames) * dt if dt else np.arange(n_frames)
        if self.System.get('temperature') is not None:
            out['Temp'] = np.asarray(self.System['temperature'], dtype=float)
        return out

    # ------------------------------------------------------------------ averaging
    def Average(self):
        """Mean (and std as error) over the run directories; a single run passes through."""
        runs = list(self.BigMeta.values())
        for key in self.Keys:
            vals = [r[key] for r in runs if key in r]
            if not vals:
                continue
            if key in ('masterkey', 'Start', 'End', 'Step', 'Skip', 'SimTime', 'Temp') or len(vals) == 1:
                self.AverageMeta[key] = vals[0]
                continue
            try:
                if isinstance(vals[0], list):  # list of (space, heights) tuples
                    n = min(len(v) for v in vals)
                    self.AverageMeta[key] = [(vals[0][i][0], np.mean([v[i][1] for v in vals], axis=0)) for i in range(n)]
                    self.Errors[key] = [(vals[0][i][0], np.std([v[i][1] for v in vals], axis=0)) for i in range(n)]
                else:
                    n = min(len(v) for v in vals)
                    stack = np.array([np.asarray(v)[:n] for v in vals], dtype=float)
                    self.AverageMeta[key] = stack.mean(axis=0)
                    self.Errors[key] = stack.std(axis=0)
            except (ValueError, TypeError) as e:
                log.warning("could not average %s across runs: %s", key, e)
                self.AverageMeta[key] = vals[0]
        return self.AverageMeta, (self.Errors if len(runs) > 1 else None)
