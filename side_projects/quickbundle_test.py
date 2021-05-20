
#Tests using the quick bundles function, get centroids, view the clusters, etc

import numpy as np
from dipy.io.streamline import load_tractogram
from dipy.tracking.streamline import Streamlines
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle
from dipy.data import get_fnames
from dipy.viz import window, actor


fname = get_fnames('fornix')
fornix = load_tractogram(fname, 'same', bbox_valid_check=False)
streamlines = fornix.streamlines

outpath = "/Users/alex/jacques/whiston_results/";
#outpath = "C:\\Users\\Jacques Stout\\Documents\\Work\\results\\"
qb = QuickBundles(threshold=10.)
clusters = qb.cluster(streamlines)

print("Nb. clusters:", len(clusters))
print("Cluster sizes:", map(len, clusters))
print("Small clusters:", clusters < 10)
print("Streamlines indices of the first cluster:\n", clusters[0].indices)
print("Centroid of the last cluster:\n", clusters[-1].centroid)

#Cluster sizes: [64, 191, 47, 1]

#Small clusters: array([False, False, False, True], dtype=bool)

interactive = True

scene = window.Scene()
scene.SetBackground(1, 1, 1)
scene.add(actor.streamtube(streamlines, window.colors.white))
window.record(scene, out_path=outpath + 'fornix_initial.png', size=(600, 600))
if interactive:
    window.show(scene)

colormap = actor.create_colormap(np.arange(len(clusters)))

scene.clear()
scene.SetBackground(1, 1, 1)
scene.add(actor.streamtube(streamlines, window.colors.white, opacity=0.05))
scene.add(actor.streamtube(clusters.centroids, colormap, linewidth=0.4))
window.record(scene, out_path= outpath + 'fornix_centroids.png', size=(600, 600))
if interactive:
    window.show(scene)

colormap_full = np.ones((len(streamlines), 3))
for cluster, color in zip(clusters, colormap):
    colormap_full[cluster.indices] = color

scene.clear()
scene.SetBackground(1, 1, 1)
scene.add(actor.streamtube(streamlines, colormap_full))
window.record(scene, out_path=outpath + 'fornix_clusters.png', size=(600, 600))
if interactive:
    window.show(scene)

save_pickle(outpath + 'QB.pkl', clusters)