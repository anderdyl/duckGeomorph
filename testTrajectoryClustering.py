import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.tracking.streamline import Streamlines
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle
from dipy.data import get_fnames
#from dipy.viz import window, actor, colormap
fname = get_fnames('fornix')
fornix = load_tractogram(fname, 'same', bbox_valid_check=False)
streamlines = fornix.streamlines