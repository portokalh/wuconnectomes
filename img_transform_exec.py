import os
import warnings
from file_tools import splitpath, mkcdir
import nibabel as nib
import numpy as np
from dipy.io.image import load_nifti, save_nifti

def header_superpose(target_path, origin_path, outpath=None, transpose=None):
    target_nii=nib.load(target_path)
    origin_nii=nib.load(origin_path)
    if np.shape(target_nii._data)[0:3] != np.shape(origin_nii._data)[0:3]:
        raise TypeError('not implemented')
    else:
        target_affine=target_nii._affine
        #added this to add a little translocation onto the atlas whenever necessary (like space_transpose). A bit convoluted
        # but didnt want to save a nifti file over itself multiple times.
        if transpose is not None:
            target_affine[:3,3] = transpose
        if np.any(target_affine != origin_nii._affine):
            new_nii = nib.Nifti1Image(origin_nii._data, target_affine, origin_nii._header)
            if outpath is None:
                outpath = origin_path
            nib.save(new_nii, outpath)

def space_transpose(origin_path, transpose=[0,0,0], outpath=None):
    if outpath is None:
        outpath = origin_path
    origin_nii=nib.load(origin_path)
    newaffine = origin_nii._affine
    newaffine[:3,3] = transpose


def get_transpose(img):
    from nibabel import load
    chass_sym_nii = load(img)
    transpose = chass_sym_nii[:3, 3]
    return transpose


def img_transform_exec(img, current_vorder, desired_vorder, output_path=None, write_transform=0, recenter=1):

    is_RGB = 0;
    is_vector = 0;
    is_tensor = 0;

    if not os.path.exists(img):
        raise('Nifti img file does not exists')


    for char in current_vorder:
        if char.islower():
            warnings.warn("Use uppercase for current order")
            current_vorder=current_vorder.upper()
    for char in desired_vorder:
        if char.islower():
            warnings.warn("Use uppercase for desired order")
            current_vorder=desired_vorder.upper()

    ordervals='RLAPSI'
    if not ordervals.find(current_vorder[0]) and not ordervals.find(current_vorder[1]) and not ordervals.find(current_vorder[2]):
        raise TypeError('Please use only R L A P S or I for current voxel order')

    if not ordervals.find(desired_vorder[0]) and not ordervals.find(desired_vorder[1]) and not ordervals.find(desired_vorder[2]):
        raise TypeError('Please use only R L A P S or I for desired voxel order')

    dirname, filename, ext = splitpath(img)
    filename=filename+"."+ext
    if output_path is None:
        output_name=filename.replace('.nii','_'+desired_vorder+'.nii')
        output_path = os.path.join(dirname,output_name)
    if not output_path.find('.'):
        mkcdir(output_path)
        output_name=filename.replace('.nii','_'+desired_vorder+'.nii')
        output_path = os.path.join(output_path,output_name)

    out_dir, _, _ = splitpath(output_path)
    affine_out = os.path.join(out_dir, current_vorder + '_to_' + desired_vorder + '_affine.pickle')

    overwrite = True;
    if os.path.isfile(output_path) and not overwrite and (not write_transform or (write_transform and os.path.exists(affine_out))):
        warnings.warn('Existing output:%s, not regenerating', output_path);

    orig_string = 'RLAPSI';
    flip_string = 'LRPAIS';
    orig_current_vorder = current_vorder;

    try:
        nii = nib.load(img)
    except:
        raise('Could not load img at '+img)
    nii_data = nii.get_data()
    hdr = nii.header

    dims = nii.shape
    if np.size(dims) > 6:
        raise('Image has > 5 dimensions')
    elif np.size(dims) < 3:
        raise('Image has < 3 dimensions')

    new_data = nii_data

    if desired_vorder!=current_vorder:
        if ((np.size(dims) > 4) and (dims(5) == 3)):
            is_vector = 1;
        elif ((np.size(dims) > 5) and (dims(5) == 6)):
            is_tensor = 1;

    #x_row = [1, 0, 0];
    #y_row = [0, 1, 0];
    #z_row = [0, 0, 1];
    affine = nii._affine
    x_row = affine[0,:]
    y_row = affine[1,:]
    z_row = affine[2,:]

    xpos=desired_vorder.find(current_vorder[0])
    if xpos == -1:
        print('Flipping first dimension')
        val=0
        new_data = np.flip(new_data, 0)
        orig_ind=orig_string.find(current_vorder[0])
        current_vorder = current_vorder[0:val] + flip_string[orig_ind] + current_vorder[val+1:]
        if is_vector:
            new_data[:,:,:,0,1]=-new_data[:,:,:,0,1]
        x_row = [-1 * val for val in x_row]

    ypos=desired_vorder.find(current_vorder[1])
    if ypos == -1:
        print('Flipping second dimension')
        val=1
        new_data = np.flip(new_data, 1)
        orig_ind=orig_string.find(current_vorder[1])
        current_vorder = current_vorder[0:val] + flip_string[orig_ind] + current_vorder[val+1:]
        if is_vector:
            new_data[:,:,:,0,2]=-new_data[:,:,:,0,2]
        y_row = [-1 * val for val in y_row]

    zpos=desired_vorder.find(current_vorder[2])
    if zpos == -1:
        print('Flipping third dimension')
        val=2
        new_data = np.flip(new_data, 2)
        orig_ind=orig_string.find(current_vorder[2])
        current_vorder = current_vorder[0:val] + flip_string[val] + current_vorder[val+1:]
        if is_vector:
            new_data[:,:,:,0,2]=-new_data[:,:,:,0,2]
        z_row = [-1 * val for val in z_row]

    xpos=current_vorder.find(desired_vorder[0])
    ypos=current_vorder.find(desired_vorder[1])
    zpos=current_vorder.find(desired_vorder[2])

    print(['Dimension order is:' + str(xpos) + ' ' + str(ypos) + ' ' + str(zpos)] )
    if not os.path.isfile(output_path) or overwrite:
        if np.size(dims) == 5:
            if is_tensor:
                new_data = new_data.tranpose(xpos, ypos, zpos, 4, 5)
            else:
                if is_vector:# =>> honestly looking at the original code, this doesnt really make sense to me, so deactivated for now. Will raise warning in case it happens
                    warnings.warn('is vector not properly implemented')
                    #    new[:,:,:,1,:] = new[:,:,:,1].transpose(xpos, ypos, zpos)
                    #new=new(:,:,:,[xpos, ypos, zpos]);
                new_data.transpose(xpos, ypos, zpos, 4, 5)
        elif np.size(dims) == 4:
            if is_RGB:
                ('is rgb not properly implemented')
                #new=new(:,:,:,[xpos, ypos, zpos]);
            new_data = new_data.transpose(xpos, zpos, ypos, 4)
        elif np.size(dims) == 3:
            new_data = new_data.transpose(xpos, ypos, zpos)

    if not os.path.isfile(affine_out) and write_transform:
        intermediate_affine_matrix = [x_row , y_row, z_row];
        iam = intermediate_affine_matrix;
        affine_matrix_for_points = [iam[xpos,:], iam[ypos,:], iam[zpos,:]]
        affine_matrix_for_images = np.inv(affine_matrix_for_points)
        am4i = affine_matrix_for_images
        affine_matrix_string = [am4i[1,:] + am4i[2,:] + am4i[3,:] + '0 0 0']
        affine_fixed_string = ['0', '0', '0',];
        try:
            #write_affine_xform_for_ants(affine_out,affine_matrix_string,affine_fixed_string);
            #needs to be implemented if this is a desired functionality
            print("nope, not implemented")
        except:
            print("nope, not implemented")

    if recenter:
        origin=np.round([val/2 for val in np.shape(new_data)])
        origin=origin[0:3]
        trueorigin=origin*[x_row[0],y_row[1],z_row[2]]
        trueorigin[2]=trueorigin[2]*(-1)
    else:
        origin=affine[0:3,3]
        trueorigin=origin*[x_row[0],y_row[1],z_row[2]]
        trueorigin[2]=trueorigin[2]*(-1)

    newaffine=np.zeros([4,4])
    newhdr=hdr
    if not desired_vorder==orig_current_vorder:

        newaffine[0:3,0:3]=affine[0:3,0:3]
        newaffine[3,:]=[0,0,0,1]
        newaffine[:3,3]=trueorigin
        #newaffine[0,:]=x_row
        #newaffine[1,:]=y_row
        #newaffine[2,:]=z_row
        #newhdr.srow_x=[newaffine[0,0:3]]
        #newhdr.srow_y=[newaffine[1,0:3]]
        #newhdr.srow_z=[newaffine[2,0:3]]
        #newhdr.pixdim=hdr.pixdim
    else:
        newaffine=affine

    new_nii=nib.Nifti1Image(new_data, newaffine, newhdr)
    output_path = str(output_path)
    nib.save(new_nii, output_path)
    """
    #newnii.hdr.dime.intent_code = nii.hdr.dime.intent_code
    new_nii=nib.Nifti1Image(new, newaffine, newhdr)
    output_pathold=output_path
    output_path = str(output_path)
    output_path=output_path.replace('.nii', '_test.nii')
    nib.save(new_nii, output_path)
    output_path2=output_path.replace('.nii','_2.nii')
    new_nii_oldaffine=nib.Nifti1Image(new, affine, hdr)
    nib.save(new_nii_oldaffine,output_path2)
    #newnii.hdr.dime.intent_code = nii.hdr.dime.intent_code
    nib_test=nib.load(output_path)
    nib_test2=nib.load(output_path)
    nib_test_target=nib.load(output_pathold)
    hdr_test=nib_test._header
    hdr_test2=nib_test2._header
    hdr_target = nib_test_target._header
    print('hi')
    """