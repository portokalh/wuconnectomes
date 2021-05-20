import numpy as np
from dipy.viz import actor, window, ui

from dipy.tracking.streamline import Streamlines
from dipy.data.fetcher import fetch_bundles_2_subjects, read_bundles_2_subjects


def win_callback(obj, event):
    global size
    if size != obj.GetSize():
        size_old = size
        size = obj.GetSize()
        size_change = [size[0] - size_old[0], 0]
        panel.re_align(size_change)

def change_slice_z(slider):
    z = int(np.round(slider.value))
    image_actor_z.display_extent(0, shape[0] - 1, 0, shape[1] - 1, z, z)


def change_slice_x(slider):
    x = int(np.round(slider.value))
    image_actor_x.display_extent(x, x, 0, shape[1] - 1, 0, shape[2] - 1)


def change_slice_y(slider):
    y = int(np.round(slider.value))
    image_actor_y.display_extent(0, shape[0] - 1, y, y, 0, shape[2] - 1)


def change_opacity(slider):
    slicer_opacity = slider.value
    image_actor_z.opacity(slicer_opacity)
    image_actor_x.opacity(slicer_opacity)
    image_actor_y.opacity(slicer_opacity)


fetch_bundles_2_subjects()

res = read_bundles_2_subjects('subj_1', ['t1', 'fa'],
                              ['af.left', 'cst.right', 'cc_1'])

streamlines = Streamlines(res['af.left'])
streamlines.extend(res['cst.right'])
streamlines.extend(res['cc_1'])

data = res['fa']
shape = data.shape
affine = res['affine']

world_coords = True

if not world_coords:
    from dipy.tracking.streamline import transform_streamlines
    streamlines = transform_streamlines(streamlines, np.linalg.inv(affine))

scene = window.Scene()
stream_actor = actor.line(streamlines)

if not world_coords:
    image_actor_z = actor.slicer(data, affine=np.eye(4))
else:
    image_actor_z = actor.slicer(data, affine)

slicer_opacity = 0.6
image_actor_z.opacity(slicer_opacity)

image_actor_x = image_actor_z.copy()
x_midpoint = int(np.round(shape[0] / 2))
image_actor_x.display_extent(x_midpoint,
                             x_midpoint, 0,
                             shape[1] - 1,
                             0,
                             shape[2] - 1)

image_actor_y = image_actor_z.copy()
y_midpoint = int(np.round(shape[1] / 2))
image_actor_y.display_extent(0, shape[0] - 1,
                             y_midpoint,
                             y_midpoint,
                             0,
                             shape[2] - 1)

show_m = window.ShowManager(scene, size=(1200, 900))
show_m.initialize()

interactive = True

scene.add(stream_actor)
scene.add(image_actor_z)
scene.add(image_actor_x)
scene.add(image_actor_y)

global size
size = scene.GetSize()

#show_m.add_window_callback(win_callback)
#show_m.render()
#show_m.start()

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



line_slider_z.on_change = fetch_bundles_2_subjects
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

"""
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
"""
"""
size = scene.GetSize()



show_m.initialize()

interactive = True

scene.zoom(1.5)
scene.reset_clipping_range()
"""

if interactive:

    show_m.add_window_callback(win_callback)
    show_m.render()
    show_m.start()

else:

    window.record(scene, out_path='bundles_and_3_slices.png', size=(1200, 900),
                  reset_camera=False)

del show_m
