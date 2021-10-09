from dipy.io.streamline import load_trk as load_trk_spe
import os

TRK_folder = '/Users/alex/jacques/AMD_TRK_testing/TRK_MDT'
for trk_path in os.listdir(TRK_folder):
    if trk_path.endswith('.trk'):
        try:
            load_trk(trk_path)
        except ValueError:
            load_trk_spe(trk_path)

