
import numpy as np
from dipy.io.image import load_nifti
import warnings
from dipy.viz import window, actor
from time import sleep
import os, glob
from dipy.segment.clustering import ClusterCentroid
from dipy.tracking.streamline import Streamlines
import vtk
import nibabel as nib


def close_window(iren):
    render_window = iren.GetRenderWindow()
    render_window.Finalize()
    iren.TerminateApp()
    del render_window, iren


def win_callback(obj, event):
    global size
    if size != obj.GetSize():
        size_old = size
        size = obj.GetSize()
        size_change = [size[0] - size_old[0], 0]
        panel.re_align(size_change)


def show_bundles(bundles, colors=None, show=True, fname=None, str_tube=False, ref=None):
    ren = window.Renderer()
    ren.SetBackground(1., 1, 1)
    if str_tube:
        bundle_actor = actor.streamtube(bundles, colors, linewidth=0.5)
        ren.add(bundle_actor)
    else:
        for (i, bundle) in enumerate(bundles):
            color = colors[i]
            #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05

            lines_actor = actor.line(bundle, color, linewidth=2.5)
            # lines_actor.RotateX(-90)
            # lines_actor.RotateZ(90)
            ren.add(lines_actor)

    if ref is not None:
        if os.path.exists(ref):
            data, affine = load_nifti(ref)
            #data, affine = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
            ref_actor = actor.slicer(data, affine)
            ren.add(ref_actor)
        else:
            txt = f'Was asked to find reference file for background at {ref} but path did not exist'
            warnings.warn(txt)

    if show:
        window.show(ren)
    if fname is not None:
        sleep(1)
        window.record(ren, n_frames=1, out_path=fname, size=(900, 900))


def launch_interactive_view(scene):
    show_m = window.ShowManager(scene, size=(1200, 900))
    show_m.initialize()
    show_m.add_window_callback(win_callback)
    show_m.render()
    show_m.start()

def setup_view_legacy(trk_object, colors=None, world_coords=False, show=True, fname=None, str_tube=False, ref=None,
               objectvals = None, colorbar=False, record = None, scene = None, interactive = True):

    from dipy.viz import actor, window, ui

    if isinstance(ref, str):
        if os.path.exists(ref):
            data, affine = load_nifti(ref)
            shape = data.shape
            # data, affine = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
        else:
            txt = f'Was asked to find reference file for background at {ref} but path did not exist'
            warnings.warn(txt)
    if isinstance(ref, np.ndarray):
        data = ref
        shape = np.shape(ref)
    if isinstance(ref, nib.Nifti1Image):
        data = np.asarray(ref.dataobj)
        shape = np.shape(data)
        affine = ref._affine
    if np.size(shape)==4:
        data = np.squeeze(data[:,:,:,0])
        shape = shape[:3]
    if not world_coords:
        image_actor_z = actor.slicer(data, affine=np.eye(4))
    else:
        image_actor_z = actor.slicer(data, affine)

    if scene is None:
        scene = window.Scene()
    #if objectvals is None:
    #    objectvals = np.random(np.shape(trk_object)[0])

    slicer_opacity = 0.9
    image_actor_z.opacity(slicer_opacity)

    #adding slicer sliders

    image_actor_x = image_actor_z.copy()
    x_midpoint = int(np.round(shape[0] / 2))
    image_actor_x.display_extent(x_midpoint,
                                 x_midpoint, 0,
                                 shape[1],
                                 0,
                                 shape[1])

    image_actor_y = image_actor_z.copy()
    y_midpoint = int(np.round(shape[1] / 2))
    image_actor_y.display_extent(0,
                                 shape[0],
                                 y_midpoint,
                                 y_midpoint,
                                 0,
                                 shape[1])


    scene.add(image_actor_z)
    scene.add(image_actor_x)
    scene.add(image_actor_y)

    if colorbar:
        bar3 = actor.scalar_bar(colors)
        scene.add(bar3)


    if isinstance(trk_object[0], ClusterCentroid):
        bundles = trk_object
        if str_tube:
            object_actor = actor.streamtube(bundles, colors, linewidth=0.5)
            #ren.add(bundle_actor)
        else:
            for (i, bundle) in enumerate(bundles):
                color = colors[i]
                #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05
                #color = (0.0, 1.0, 0.0)
                object_actor = actor.line(bundle, color, linewidth=1.0)
                # lines_actor.RotateX(-90)
                # lines_actor.RotateZ(90)
                #scene.add(object_actor)
    elif isinstance(trk_object[0][0], ClusterCentroid):
        for group in np.arange(np.shape(trk_object)[0]):
            color = colors[group]
            bundles = trk_object[group]
            for (i, bundle) in enumerate(bundles):
                #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05
                # color = (0.0, 1.0, 0.0)
                object_actor = actor.line(bundle, color, linewidth=1.0)
                # lines_actor.RotateX(-90)
                # lines_actor.RotateZ(90)
                #scene.add(object_actor)
    elif isinstance(trk_object, Streamlines):
        if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
            object_actor = actor.line(trk_object,objectvals, linewidth=0.1,
                               lookup_colormap=colors)
        else:
            object_actor = actor.line(trk_object)

        #scene.add(object_actor)

    else:
        raise Exception('Unindentified object')

    scene.add(object_actor)



    """
    if os.path.exists(ref):
        data, affine = load_nifti(ref)
        shape = data.shape
        # data, affine = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
    else:
        txt = f'Was asked to find reference file for background at {ref} but path did not exist'
        warnings.warn(txt)

    if not world_coords:
        image_actor_z = actor.slicer(data, affine=np.eye(4))
    else:
        image_actor_z = actor.slicer(data, affine)

    slicer_opacity = 0.6
    image_actor_z.opacity(slicer_opacity)

    #adding slicer sliders

    image_actor_x = image_actor_z.copy()
    x_midpoint = int(np.round(shape[0] / 2))
    image_actor_x.display_extent(x_midpoint,
                                 x_midpoint, 0,
                                 shape[1],
                                 0,
                                 shape[2])

    image_actor_y = image_actor_z.copy()
    y_midpoint = int(np.round(shape[1] / 2))
    image_actor_y.display_extent(0,
                                 shape[0],
                                 y_midpoint,
                                 y_midpoint,
                                 0,
                                 shape[2])


    scene.add(object_actor)
    scene.add(image_actor_z)
    scene.add(image_actor_x)
    scene.add(image_actor_y)

    if colorbar:
        bar3 = actor.scalar_bar(colors)
        scene.add(bar3)
    """
    show_m = window.ShowManager(scene, size=(2000, 900))
    show_m.initialize()

    line_slider_z = ui.LineSlider2D(min_value=0,
                                    max_value=shape[2] - 1,
                                    initial_value=shape[2] / 2,
                                    text_template="{value:.0f}",
                                    length=140)

    line_slider_x = ui.LineSlider2D(min_value=0,
                                    max_value=shape[0] - 1,
                                    initial_value=shape[0] / 2,
                                    text_template="{value:.0f}",
                                    length=140)

    line_slider_y = ui.LineSlider2D(min_value=0,
                                    max_value=shape[1] - 1,
                                    initial_value=shape[1] / 2,
                                    text_template="{value:.0f}",
                                    length=140)

    opacity_slider = ui.LineSlider2D(min_value=0.0,
                                     max_value=1.0,
                                     initial_value=slicer_opacity,
                                     length=140)

    def change_slice_z(slider):
        z = int(np.round(slider.value))
        image_actor_z.display_extent(0, shape[0] - 1, 0, shape[1] - 1, z, z)

    def change_slice_x(slider):
        x = int(np.round(slider.value))
        image_actor_x.display_extent(x, x, 0, shape[1] - 1, 0, shape[1] - 1)

    def change_slice_y(slider):
        y = int(np.round(slider.value))
        image_actor_y.display_extent(0, shape[0] - 1, y, y, 0, shape[1] - 1)

    def change_opacity(slider):
        slicer_opacity = slider.value
        image_actor_z.opacity(slicer_opacity)
        image_actor_x.opacity(slicer_opacity)
        image_actor_y.opacity(slicer_opacity)

    line_slider_z.on_change = change_slice_z
    line_slider_x.on_change = change_slice_x
    line_slider_y.on_change = change_slice_y
    opacity_slider.on_change = change_opacity

    def build_label(text):
        label = ui.TextBlock2D()
        label.message = text
        label.font_size = 18
        label.font_family = 'Arial'
        label.justification = 'left'
        label.bold = False
        label.italic = False
        label.shadow = False
        label.background_color = (0, 0, 0)
        label.color = (1, 1, 1)

        return label

    line_slider_label_z = build_label(text="Z Slice")
    line_slider_label_x = build_label(text="X Slice")
    line_slider_label_y = build_label(text="Y Slice")
    opacity_slider_label = build_label(text="Opacity")

    panel = ui.Panel2D(size=(300, 200),
                       color=(1, 1, 1),
                       opacity=0.1,
                       align="right")
    panel.center = (1030, 120)

    panel.add_element(line_slider_label_x, (0.1, 0.75))
    panel.add_element(line_slider_x, (0.38, 0.75))
    panel.add_element(line_slider_label_y, (0.1, 0.55))
    panel.add_element(line_slider_y, (0.38, 0.55))
    panel.add_element(line_slider_label_z, (0.1, 0.35))
    panel.add_element(line_slider_z, (0.38, 0.35))
    panel.add_element(opacity_slider_label, (0.1, 0.15))
    panel.add_element(opacity_slider, (0.38, 0.15))

    scene.add(panel)

    global size
    size = scene.GetSize()

    scene.zoom(1.5)
    scene.reset_clipping_range()

    if interactive:

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()

    if record is not None:
        if os.path.exists(record):
            record_name = os.path.basename(record)
            dir_name = os.path.dirname(record)
            record_name = record_name.replace('.','_1.')
            record = os.path.join(dir_name, record_name)
            if os.path.exists(record):
                while os.path.exists(record):
                    if record_name.split('.')[0].split('_')[-1].isnumeric():
                        val = int((record_name.split('.')[0].split('_'))[-1])
                        newval = val + 1
                        record_name = record_name.replace('_' + str(val)+'.', '_' + str(newval)+'.')
                    else:
                        raise Exception('wtf??')
                    record = os.path.join(dir_name, record_name)
        window.record(scene, out_path=record, size=(1200, 900),
                      reset_camera=False)
        print(f'Saved figure at {record}')

    scene.rm(object_actor)
    return scene


def view_test(scene,testtype='interactive',record_path = None):

    from dipy.viz import actor, window, ui
    import warnings

    if testtype=='interactive':
        show_m = window.ShowManager(scene, size=(2000, 900))
        show_m.initialize()

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()
        del show_m
    elif testtype=='record':
        window.record(scene, out_path=record_path, size=(1200, 900),
                      reset_camera=False)
        print(f'Saved figure at {record_path}')
    else:
        warnings.warn('invalid test')
    return


def setup_view(trk_object, colors=None, world_coords=False, show=True, fname=None, str_tube=False, ref=None,
               objectvals = None, colorbar=False, record = None, scene = None, plane = 'all', interactive = True):

    from dipy.viz import actor, window, ui

    if isinstance(ref, str):
        if os.path.exists(ref):
            data, affine = load_nifti(ref)
            shape = data.shape
            # data, affine = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
        else:
            txt = f'Was asked to find reference file for background at {ref} but path did not exist'
            raise Exception(txt)
    if isinstance(ref, np.ndarray):
        data = ref
        shape = np.shape(ref)
    if isinstance(ref, nib.Nifti1Image):
        data = np.asarray(ref.dataobj)
        shape = np.shape(data)
        affine = ref._affine
    if np.size(shape)==4:
        data = np.squeeze(data[:,:,:,0])
        shape = shape[:3]
    value_range = (0,1500)
    if not world_coords:
        image_actor_z = actor.slicer(data, affine=np.eye(4),value_range = value_range)
    else:
        image_actor_z = actor.slicer(data, affine,value_range = value_range)

    object_actors_toremove = []

    def select_plane(plane):
        show_x = False
        show_y = False
        show_z = False
        if plane == 'all':
            show_x = True
            show_y = True
            show_z = True
        elif plane == 'x':
            show_x=True
        elif plane == 'x':
            show_y=True
        elif plane == 'x':
            show_z=True
        else:
            warnings.warn('Plane was not a recognized parameter, showing all planes by default')
            show_x = True
            show_y = True
            show_z = True
        return show_x, show_y, show_z

    def change_slice_z(slider):
        z = int(np.round(slider.value))
        image_actor_z.display_extent(0, shape[0] - 1, 0, shape[1] - 1, z, z)

    def change_slice_x(slider):
        x = int(np.round(slider.value))
        image_actor_x.display_extent(x, x, 0, shape[1] - 1, 0, shape[1] - 1)

    def change_slice_y(slider):
        y = int(np.round(slider.value))
        image_actor_y.display_extent(0, shape[0] - 1, y, y, 0, shape[1] - 1)

    def change_opacity(slider):
        slicer_opacity = slider.value
        image_actor_z.opacity(slicer_opacity)
        image_actor_x.opacity(slicer_opacity)
        image_actor_y.opacity(slicer_opacity)

    def build_label(text, color=(1,1,1)):
        label = ui.TextBlock2D()
        label.message = text
        label.font_size = 18
        label.font_family = 'Arial'
        label.justification = 'left'
        label.bold = False
        label.italic = False
        label.shadow = False
        label.background_color = (0, 0, 0)
        label.color = color

        return label

    if scene is None:
        scene = window.Scene()
        #scene = window.Renderer()
        # if objectvals is None:
        #    objectvals = np.random(np.shape(trk_object)[0])

        slicer_opacity = 0.7
        image_actor_z.opacity(slicer_opacity)

        # adding slicer sliders

        image_actor_x = image_actor_z.copy()
        x_midpoint = int(np.round(shape[0] / 2))
        image_actor_x.display_extent(x_midpoint,
                                     x_midpoint, 0,
                                     shape[1],
                                     0,
                                     shape[1])

        image_actor_y = image_actor_z.copy()
        y_midpoint = int(np.round(shape[1] / 2))
        image_actor_y.display_extent(0,
                                     shape[0],
                                     y_midpoint,
                                     y_midpoint,
                                     0,
                                     shape[1])

        show_x, show_y, show_z = select_plane(plane)

        if show_z:
            scene.add(image_actor_z)
        if show_x:
            scene.add(image_actor_x)
        if show_y:
            scene.add(image_actor_y)

        if colorbar:
            if colors is not None:
                bar3 = actor.scalar_bar(colors)
                scene.add(bar3)
                object_actors_toremove.append(bar3)
            else:
                if np.shape(objectvals)[0] == np.size(trk_object):
                    for i,color in enumerate(objectvals):
                        colorobject = actor.text_3d(text=f'Bundle {str(i + 1)}', position=(150, -i * 10, 0), color=color, justification = 'right', vertical_justification = 'middle')
                        scene.add(colorobject)
                        object_actors_toremove.append(colorobject)
                else:
                    raise Exception('Unrecognized circumstances where the given colors do not correspond to the number of trk objects, to add a case here for this circumstance')

        show_m = window.ShowManager(scene, size=(2000, 900))
        show_m.initialize()

        line_slider_z = ui.LineSlider2D(min_value=0,
                                        max_value=shape[2] - 1,
                                        initial_value=shape[2] / 2,
                                        text_template="{value:.0f}",
                                        length=140)

        line_slider_x = ui.LineSlider2D(min_value=0,
                                        max_value=shape[0] - 1,
                                        initial_value=shape[0] / 2,
                                        text_template="{value:.0f}",
                                        length=140)

        line_slider_y = ui.LineSlider2D(min_value=0,
                                        max_value=shape[1] - 1,
                                        initial_value=shape[1] / 2,
                                        text_template="{value:.0f}",
                                        length=140)

        opacity_slider = ui.LineSlider2D(min_value=0.0,
                                         max_value=1.0,
                                         initial_value=slicer_opacity,
                                         length=140)

        panel = ui.Panel2D(size=(300, 200),
                           color=(1, 1, 1),
                           opacity=0.1,
                           align="right")
        panel.center = (1030, 120)


        if show_z:
            line_slider_z.on_change = change_slice_z
            line_slider_label_z = build_label(text="Z Slice")
            panel.add_element(line_slider_label_z, (0.1, 0.35))
            panel.add_element(line_slider_z, (0.38, 0.35))
        if show_x:
            line_slider_x.on_change = change_slice_x
            line_slider_label_x = build_label(text="X Slice")
            panel.add_element(line_slider_label_x, (0.1, 0.75))
            panel.add_element(line_slider_x, (0.38, 0.75))
        if show_y:
            line_slider_y.on_change = change_slice_y
            line_slider_label_y = build_label(text="Y Slice")
            panel.add_element(line_slider_label_y, (0.1, 0.55))
            panel.add_element(line_slider_y, (0.38, 0.55))

        opacity_slider.on_change = change_opacity
        opacity_slider_label = build_label(text="Opacity")

        panel.add_element(opacity_slider_label, (0.1, 0.15))
        panel.add_element(opacity_slider, (0.38, 0.15))

        scene.add(panel)

        global size
        size = scene.GetSize()

        scene.zoom(1.5)
        scene.reset_clipping_range()
    else:
        if colorbar:
            if colors is not None:
                bar3 = actor.scalar_bar(colors)
                scene.add(bar3)
                object_actors_toremove.append(bar3)
            else:
                if np.shape(objectvals)[0] == np.size(trk_object):
                    for i,color in enumerate(objectvals):
                        colorobject = actor.text_3d(text=f'Bundle {str(i + 1)}', position=(150, -i * 10, 0), color=color, justification = 'right', vertical_justification = 'middle')
                        scene.add(colorobject)
                        object_actors_toremove.append(colorobject)
                else:
                    raise Exception('Unrecognized circumstances where the given colors do not correspond to the number of trk objects, to add a case here for this circumstance')
    """
    test=1
    if test==1:
        bundle = trk_object[0]
        streamline = trk_object[0][0]
        str_vals = objectvals[0][0]
        #object_actor = actor.line(streamline,str_vals, linewidth=0.2,
        #                          lookup_colormap=colors)
        #object_actor = actor.line(bundle,window.colors.green, linewidth=0.2)

        colors_points = []
        for s in range(len(bundle)):
            stream = bundle[s]
            for idx in range(len(stream)):
                colors_points.append(objectvals[0][s][idx])

        colors_new = []
        for s in range(len(bundle)):
            stream = bundle[s]
            colors_new.append(list(objectvals[0][s]))

        object_actor = actor.line(bundle, colors_points, linewidth=0.2, lookup_colormap=colors)
        scene.add(object_actor)

    if test==1:
        bundles = trk_object
        if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
            for bundle, vals in zip(trk_object, objectvals):
                colors_points = []
                for s in range(len(bundle)):
                    stream = bundle[s]
                    for idx in range(len(stream)):
                        colors_points.append(vals[s][idx])
                object_actor = actor.line(bundle, colors_points, linewidth=0.2, lookup_colormap=colors)
                scene.add(object_actor)
    """

    if isinstance(trk_object, ClusterCentroid):
        if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
            object_actor = actor.line(trk_object, objectvals, linewidth=0.1,
                                      lookup_colormap=colors)
            scene.add(object_actor)
            object_actors_toremove.append(object_actor)
        else:
            object_actor = actor.line(trk_object)
            scene.add(object_actor)
            object_actors_toremove.append(object_actor)
    elif isinstance(trk_object[0], ClusterCentroid):
        bundles = trk_object
        if str_tube:
            object_actor = actor.streamtube(bundles, colors, linewidth=0.5)
            # ren.add(bundle_actor)
        if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
            if np.size(np.shape(objectvals[0]))<=1: #if object vals just contains a list of values where the number of values is the same as the number of bundles, then do this (color per bundle)
                for (i, bundle) in enumerate(bundles):
                    object_actor = actor.line(bundle, objectvals[i], linewidth=0.1,
                                              lookup_colormap=colors)
                    scene.add(object_actor)
                    object_actors_toremove.append(object_actor)
            else:
                for bundle, vals in zip(trk_object, objectvals):
                    colors_points = []
                    for s in range(len(bundle)):
                        stream = bundle[s]
                        for idx in range(len(stream)):
                            colors_points.append(vals[s][idx])
                    object_actor = actor.line(bundle, colors_points, linewidth=0.2, lookup_colormap=colors)
                    scene.add(object_actor)
                    object_actors_toremove.append(object_actor)
        else:
            for (i, bundle) in enumerate(bundles):
                color = colors[i]
                #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05
                # color = (0.0, 1.0, 0.0)
                object_actor = actor.line(bundle, color, linewidth=1.0)
                # lines_actor.RotateX(-90)
                # lines_actor.RotateZ(90)
                scene.add(object_actor)
                object_actors_toremove.append(object_actor)
    elif isinstance(trk_object[0][0], ClusterCentroid):
        for group in np.arange(np.shape(trk_object)[0]):
            color = colors[group]
            bundles = trk_object[group]
            for (i, bundle) in enumerate(bundles):
                #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05
                # color = (0.0, 1.0, 0.0)
                object_actor = actor.line(bundle, color, linewidth=1.0)
                # lines_actor.RotateX(-90)
                # lines_actor.RotateZ(90)
                scene.add(object_actor)
                object_actors_toremove.append(object_actor)
    elif isinstance(trk_object, Streamlines):
        if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
            object_actor = actor.line(trk_object, objectvals, linewidth=0.1,
                                      lookup_colormap=colors)
        else:
            object_actor = actor.line(trk_object)
        scene.add(object_actor)
        object_actors_toremove.append(object_actor)

    elif isinstance(trk_object[0], Streamlines):
        #if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
        if objectvals[0] is not None:
            if np.size(np.shape(objectvals[0]))>1:      #treat objet vals as point coloring
                for bundle, vals in zip(trk_object, objectvals):
                    colors_points = []
                    for s in range(len(bundle)):
                        stream = bundle[s]
                        for idx in range(len(stream)):
                            colors_points.append(vals[s][idx])
                    object_actor = actor.line(bundle, colors_points, linewidth=0.2, lookup_colormap=colors)
                    scene.add(object_actor)
                    object_actors_toremove.append(object_actor)

            elif np.size(np.shape(objectvals[0]))<=1:### treat objectvals as line coloring
                for bundle, vals in zip(trk_object, objectvals):
                    object_actor = actor.line(bundle, vals, linewidth=0.1,
                                                  lookup_colormap=colors)
                    scene.add(object_actor)
                    object_actors_toremove.append(object_actor)
        else:
            object_actor = actor.line(trk_object)
            scene.add(object_actor)
            object_actors_toremove.append(object_actor)

    else:
        raise Exception('Unindentified object')

    #scene.add(object_actor)

    if interactive:
        if not 'show_m' in locals():
            show_m = window.ShowManager(scene, size=(2000, 900))
            show_m.initialize()

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()
        del show_m

    if record is not None:
        if os.path.exists(record):
            record_name = os.path.basename(record)
            dir_name = os.path.dirname(record)
            record_name = record_name.replace('.','_1.')
            record = os.path.join(dir_name, record_name)
            if os.path.exists(record):
                while os.path.exists(record):
                    if record_name.split('.')[0].split('_')[-1].isnumeric():
                        val = int((record_name.split('.')[0].split('_'))[-1])
                        newval = val + 1
                        record_name = record_name.replace('_' + str(val)+'.', '_' + str(newval)+'.')
                    else:
                        raise Exception('wtf??')
                    record = os.path.join(dir_name, record_name)
        window.record(scene, out_path=record, size=(1200, 900),
                      reset_camera=False)
        print(f'Saved figure at {record}')

    #scene.rm(object_actor)
    #scene.clear()
    for object_actor in object_actors_toremove:
        scene.rm(object_actor)
    return scene


def setup_view_colortest(trk_object, colors=None, world_coords=False, show=True, fname=None, str_tube=False, ref=None,
               objectvals = None, colorbar=False, record = None, scene = None, plane = 'all', interactive = True):

    from dipy.viz import actor, window, ui

    object_actors_toremove = []

    def select_plane(plane):
        show_x = False
        show_y = False
        show_z = False
        if plane == 'all':
            show_x = True
            show_y = True
            show_z = True
        elif plane == 'x':
            show_x=True
        elif plane == 'x':
            show_y=True
        elif plane == 'x':
            show_z=True
        else:
            warnings.warn('Plane was not a recognized parameter, showing all planes by default')
            show_x = True
            show_y = True
            show_z = True
        return show_x, show_y, show_z

    if scene is None:
        scene = window.Scene()

        if colorbar:
            if colors is not None:
                bar3 = actor.scalar_bar(colors)
                scene.add(bar3)
                object_actors_toremove.append(bar3)
            else:
                objectvals = [objectvals[0],objectvals[1]]
                for i,color in enumerate(objectvals):
                    #axis1 = actor.axes(colorx=color[0],colory=color[1],colorz=color[2])
                    #colorobject = actor.texture(color)
                    #colorobject = billboard(np.array([50,50,50],[10,10,10]))
                    #colorobject = actor.label(text=f'Bundle {str(i+1)}', pos=(0, i*1, 0), scale=(0.1, 0.1, 0.1), color=color)
                    #colorobject = actor.text_3d(text=f'Bundle {str(i + 1)}', position=(0, -i * 10, 0), color=color)
                    line_slider_label_z = build_label(text="Z Slice")
                    panel.add_element(line_slider_label_z, (0.1, 0.35))
                    scene.add(colorobject)
                    #object_actors_toremove.append(colorobject)

        show_m = window.ShowManager(scene, size=(2000, 900))
        show_m.initialize()

        panel = ui.Panel2D(size=(300, 200),
                           color=(1, 1, 1),
                           opacity=0.1,
                           align="right")
        panel.center = (1030, 120)

        scene.add(panel)

        global size
        size = scene.GetSize()

        scene.zoom(1.5)
        scene.reset_clipping_range()
    else:
        if colorbar:
            bar3 = actor.scalar_bar(colors)
            scene.add(bar3)
            object_actors_toremove.append(bar3)

    if interactive:
        if not 'show_m' in locals():
            show_m = window.ShowManager(scene, size=(2000, 900))
            show_m.initialize()

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()
        del show_m

    for object_actor in object_actors_toremove:
        scene.rm(object_actor)
    return scene



def setup_view_scene_experiment(trk_object, colors=None, world_coords=False, show=True, fname=None, str_tube=False, ref=None,
               objectvals = None, colorbar=False, record = None, scene = None, interactive = True):

    from dipy.viz import actor, window, ui

    if isinstance(ref, str):
        if os.path.exists(ref):
            data, affine = load_nifti(ref)
            shape = data.shape
            # data, affine = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
        else:
            txt = f'Was asked to find reference file for background at {ref} but path did not exist'
            warnings.warn(txt)
    if isinstance(ref, np.ndarray):
        data = ref
        shape = np.shape(ref)
    if isinstance(ref, nib.Nifti1Image):
        data = np.asarray(ref.dataobj)
        shape = np.shape(data)
        affine = ref._affine
    if np.size(shape)==4:
        data = np.squeeze(data[:,:,:,0])
        shape = shape[:3]
    if not world_coords:
        image_actor_z = actor.slicer(data, affine=np.eye(4))
    else:
        image_actor_z = actor.slicer(data, affine)

    if scene is None:
        scene = window.Scene()
        #if objectvals is None:
        #    objectvals = np.random(np.shape(trk_object)[0])

        slicer_opacity = 0.6
        image_actor_z.opacity(slicer_opacity)

        #adding slicer sliders

        image_actor_x = image_actor_z.copy()
        x_midpoint = int(np.round(shape[0] / 2))
        image_actor_x.display_extent(x_midpoint,
                                     x_midpoint, 0,
                                     shape[1],
                                     0,
                                     shape[1])

        image_actor_y = image_actor_z.copy()
        y_midpoint = int(np.round(shape[1] / 2))
        image_actor_y.display_extent(0,
                                     shape[0],
                                     y_midpoint,
                                     y_midpoint,
                                     0,
                                     shape[1])


        scene.add(image_actor_z)
        scene.add(image_actor_x)
        scene.add(image_actor_y)

        if colorbar:
            bar3 = actor.scalar_bar(colors)
            scene.add(bar3)

        show_m = window.ShowManager(scene, size=(2000, 900))
        show_m.initialize()

        line_slider_z = ui.LineSlider2D(min_value=0,
                                        max_value=shape[2] - 1,
                                        initial_value=shape[2] / 2,
                                        text_template="{value:.0f}",
                                        length=140)

        line_slider_x = ui.LineSlider2D(min_value=0,
                                        max_value=shape[0] - 1,
                                        initial_value=shape[0] / 2,
                                        text_template="{value:.0f}",
                                        length=140)

        line_slider_y = ui.LineSlider2D(min_value=0,
                                        max_value=shape[1] - 1,
                                        initial_value=shape[1] / 2,
                                        text_template="{value:.0f}",
                                        length=140)

        opacity_slider = ui.LineSlider2D(min_value=0.0,
                                         max_value=1.0,
                                         initial_value=slicer_opacity,
                                         length=140)

        def change_slice_z(slider):
            z = int(np.round(slider.value))
            image_actor_z.display_extent(0, shape[0] - 1, 0, shape[1] - 1, z, z)

        def change_slice_x(slider):
            x = int(np.round(slider.value))
            image_actor_x.display_extent(x, x, 0, shape[1] - 1, 0, shape[1] - 1)

        def change_slice_y(slider):
            y = int(np.round(slider.value))
            image_actor_y.display_extent(0, shape[0] - 1, y, y, 0, shape[1] - 1)

        def change_opacity(slider):
            slicer_opacity = slider.value
            image_actor_z.opacity(slicer_opacity)
            image_actor_x.opacity(slicer_opacity)
            image_actor_y.opacity(slicer_opacity)

        line_slider_z.on_change = change_slice_z
        line_slider_x.on_change = change_slice_x
        line_slider_y.on_change = change_slice_y
        opacity_slider.on_change = change_opacity

        def build_label(text):
            label = ui.TextBlock2D()
            label.message = text
            label.font_size = 18
            label.font_family = 'Arial'
            label.justification = 'left'
            label.bold = False
            label.italic = False
            label.shadow = False
            label.background_color = (0, 0, 0)
            label.color = (1, 1, 1)

            return label

        line_slider_label_z = build_label(text="Z Slice")
        line_slider_label_x = build_label(text="X Slice")
        line_slider_label_y = build_label(text="Y Slice")
        opacity_slider_label = build_label(text="Opacity")

        panel = ui.Panel2D(size=(300, 200),
                           color=(1, 1, 1),
                           opacity=0.1,
                           align="right")
        panel.center = (1030, 120)

        panel.add_element(line_slider_label_x, (0.1, 0.75))
        panel.add_element(line_slider_x, (0.38, 0.75))
        panel.add_element(line_slider_label_y, (0.1, 0.55))
        panel.add_element(line_slider_y, (0.38, 0.55))
        panel.add_element(line_slider_label_z, (0.1, 0.35))
        panel.add_element(line_slider_z, (0.38, 0.35))
        panel.add_element(opacity_slider_label, (0.1, 0.15))
        panel.add_element(opacity_slider, (0.38, 0.15))

        scene.add(panel)

    ### put back here




    # scene.add(panel)
    ### end of panel segment

    global size
    size = scene.GetSize()

    scene.zoom(1.5)
    scene.reset_clipping_range()

    if isinstance(trk_object[0], ClusterCentroid):
        bundles = trk_object
        if str_tube:
            object_actor = actor.streamtube(bundles, colors, linewidth=0.5)
            #ren.add(bundle_actor)
        else:
            for (i, bundle) in enumerate(bundles):
                color = colors[i]
                #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05
                #color = (0.0, 1.0, 0.0)
                object_actor = actor.line(bundle, color, linewidth=1.0)
                # lines_actor.RotateX(-90)
                # lines_actor.RotateZ(90)
                scene.add(object_actor)
    elif isinstance(trk_object[0][0], ClusterCentroid):
        for group in np.arange(np.shape(trk_object)[0]):
            color = colors[group]
            bundles = trk_object[group]
            for (i, bundle) in enumerate(bundles):
                #         lines_actor = actor.streamtube(bundle, color, linewidth=0.05
                # color = (0.0, 1.0, 0.0)
                object_actor = actor.line(bundle, color, linewidth=1.0)
                # lines_actor.RotateX(-90)
                # lines_actor.RotateZ(90)
                scene.add(object_actor)
    elif isinstance(trk_object, Streamlines):
        if isinstance(colors, vtk.vtkLookupTable) and objectvals[0] is not None:
            object_actor = actor.line(trk_object,objectvals, linewidth=0.1,
                               lookup_colormap=colors)
        else:
            object_actor = actor.line(trk_object)

        scene.add(object_actor)

    else:
        raise Exception('Unindentified object')

    scene.add(object_actor)

    """
    if os.path.exists(ref):
        data, affine = load_nifti(ref)
        shape = data.shape
        # data, affine = load_nifti('/Volumes/Data/Badea/Lab/mouse/VBM_19IntractEP01_IITmean_RPI-work/dwi/SyN_0p5_3_0p5_dwi/dwiMDT_NoNameYet_n7_i6/median_images/MDT_fa.nii.gz')
    else:
        txt = f'Was asked to find reference file for background at {ref} but path did not exist'
        warnings.warn(txt)

    if not world_coords:
        image_actor_z = actor.slicer(data, affine=np.eye(4))
    else:
        image_actor_z = actor.slicer(data, affine)

    slicer_opacity = 0.6
    image_actor_z.opacity(slicer_opacity)

    #adding slicer sliders

    image_actor_x = image_actor_z.copy()
    x_midpoint = int(np.round(shape[0] / 2))
    image_actor_x.display_extent(x_midpoint,
                                 x_midpoint, 0,
                                 shape[1],
                                 0,
                                 shape[2])

    image_actor_y = image_actor_z.copy()
    y_midpoint = int(np.round(shape[1] / 2))
    image_actor_y.display_extent(0,
                                 shape[0],
                                 y_midpoint,
                                 y_midpoint,
                                 0,
                                 shape[2])


    scene.add(object_actor)
    scene.add(image_actor_z)
    scene.add(image_actor_x)
    scene.add(image_actor_y)

    if colorbar:
        bar3 = actor.scalar_bar(colors)
        scene.add(bar3)
    """



    if interactive:

        show_m.add_window_callback(win_callback)
        show_m.render()
        show_m.start()

    if record is not None:
        if os.path.exists(record):
            record_name = os.path.basename(record)
            dir_name = os.path.dirname(record)
            record_name = record_name.replace('.','_1.')
            record = os.path.join(dir_name, record_name)
            if os.path.exists(record):
                while os.path.exists(record):
                    if record_name.split('.')[0].split('_')[-1].isnumeric():
                        val = int((record_name.split('.')[0].split('_'))[-1])
                        newval = val + 1
                        record_name = record_name.replace('_' + str(val), '_' + str(newval))
                    else:
                        raise Exception('wtf??')
                    record = os.path.join(dir_name, record_name)
        window.record(scene, out_path=record, size=(1200, 900),
                      reset_camera=False)
        print(f'Saved figure at {record}')

    scene.clear()
    return scene