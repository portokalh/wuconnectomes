
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
import cv2

def save_slice_volume(volume, save_path):
    '''
    the function save volume data to slices in the specific directory
    :param volume: input volume data
    :param save_path:
    :return:
    '''
    shape = volume.shape
    # translate intensity to 0-255
    v_max = np.max(volume)
    v_min = np.min(volume)
    volume_norm = (volume - v_min) / (v_max - v_min)
    volume_norm = (volume_norm * 255).astype("int")
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    for i in range(shape[-1]):
        abs_path = os.path.join(save_path, str(i)+".png")
        cv.imwrite(abs_path, volume_norm[..., i])

def show_slices(slices):
    """ Function to display row of image slices """
    fig, axes = plt.subplots(1, len(slices))
    for i, slice in enumerate(slices):
        axes[i].imshow(slice.T, cmap="gray", origin="lower")

img = nib.load("/Users/alex/brain_data/atlases/chass_symmetric3/chass_symmetric3_FA.nii.gz")
img_data = img.get_fdata()
img_affine = img.affine
save_path = "/Users/alex/jacques/testslice.png"
slice_0 = img_data[64, :, :]
slice_1 = img_data[:, 64, :]
slice_2 = img_data[:, :, 64]
nft_img = nib.Nifti1Image(slice_0, img_affine)
nib.save(nft_img, save_path)

#show_slices([slice_0,slice_1,slice_2])

"""
fig, axes = plt.subplots(1, 1)
axes[0].imshow(slice.T, cmap="gray", origin="lower")
fig.savefig("/Users/alex/jacques/testslice.png")
nft_img = nib.Nifti1Image(slice_img, img_affine)
nib.save(nft_img, "/Users/alex/jacques/testslice.png")
"""
