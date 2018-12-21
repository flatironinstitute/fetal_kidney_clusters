import numpy as np
import os
from scipy import signal

interior_sequence = [
    (0, 0, 1),
    (1, 0, 1),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 1, 1)
]

tetrahedral_vertices = []
last = interior_sequence[-1]
for vertex in interior_sequence:
    vertices = [(1,1,1), vertex, last, (0,0,0)]
    tetrahedral_vertices.append(vertices)
    last = vertex
tetrahedral_vertices = np.array(tetrahedral_vertices)

def threshold_outline(array, threshhold=0.5, strides=1):
    if strides > 1:
        array = array[::strides, ::strides, ::strides]
    (nrows, ncols, npiles) = array.shape
    outlines = np.zeros(array.shape)
    maxneighbors = np.zeros(array.shape)
    minneighbors = np.zeros(array.shape) + 1
    for di in (0,1):
        for dj in (0,1):
            for dk in (0,1):
                minneighbors[:-1,:-1,:-1] = np.minimum(
                    minneighbors[:-1,:-1,:-1],
                    array[di: nrows-1+di, dj: ncols-1+dj, dk:npiles-1+dk]
                )
                maxneighbors[:-1,:-1,:-1] = np.maximum(
                    maxneighbors[:-1,:-1,:-1],
                    array[di: nrows-1+di, dj: ncols-1+dj, dk:npiles-1+dk]
                )
    outlines = (minneighbors < threshhold) & (maxneighbors > threshhold)
    outlines = outlines.astype(np.int)
    return (outlines, array)

def blur_category(array, extent=10, scale=0.3):
    nx, ny, nz = (extent, extent, extent)
    x = np.linspace(-extent, extent, nx)
    y = np.linspace(-extent, extent, ny)
    z = np.linspace(-extent, extent, nz)
    xv, yv, zv = np.meshgrid(x, y, z)
    kernel3 = 1.0/(1.0 + scale * (xv*xv + yv*yv + zv*zv))
    return signal.fftconvolve(array, kernel3, mode='valid')


def dump_named_json_vector(name, vector, f, comma=True, indent=""):
    f.write(indent)
    f.write('"%s": [' % name)
    inside = False
    for value in vector:
        if inside:
            f.write(",")
        f.write(str(int(value)))
        inside = True
    f.write("]")
    if comma:
        f.write(",")
    f.write("\n")


def nonzero_tetrahedra(array, limit=10000000, strides=1):
    Abuffer = []
    Bbuffer = []
    Cbuffer = []
    Dbuffer = []
    fbuffer = []
    (outlines, array) = threshold_outline(array, threshhold=0.5, strides=strides)
    (i_indices, j_indices, k_indices) = outlines.nonzero()
    assert len(i_indices) > 0
    for count in range(len(i_indices)):
        i = i_indices[count]
        j = j_indices[count]
        k = k_indices[count]
        ijk = np.array([[i,j,k]])
        for tetrahedron in tetrahedral_vertices:
            assert len(tetrahedron) == 4
            ijk_tet = tetrahedron + ijk
            values = array[ijk_tet[:,0], ijk_tet[:,1], ijk_tet[:,2]]
            if values.max() > 0.5 and values.min() < 0.5:
                fbuffer.extend(values)
                Abuffer.extend(ijk_tet[0])
                Bbuffer.extend(ijk_tet[1])
                Cbuffer.extend(ijk_tet[2])
                Dbuffer.extend(ijk_tet[3])
        if limit and len(Abuffer) > limit:
            break
    center = (i_indices.mean(), j_indices.mean(), k_indices.mean())
    return (Abuffer, Bbuffer, Cbuffer, Dbuffer, fbuffer, center)

special_categories = {1,3,8,9,10,11,13}

def dump_tetrahedral_data(to_file, indent, A, B, C, D, f, center,
        hexcolor="0xffff00", comma=True, fscale=20, name="shape", special_categories=special_categories):
    indent2 = "    " + indent
    to_file.write(indent + "{\n")
    to_file.write(indent2 + '"hexcolor": "%s",\n'  % (hexcolor))
    to_file.write(indent2 + '"category": "%s",\n'  % (name))
    to_file.write(indent2 + '"center": %s,\n'  % (repr(list(center))))
    opacity = 0
    if name in special_categories:
        opacity = 1.0
    to_file.write(indent2 + '"opacity": %s,\n'  % (opacity))
    dump_named_json_vector("A", A, to_file, indent=indent2)
    dump_named_json_vector("B", B, to_file, indent=indent2)
    dump_named_json_vector("C", C, to_file, indent=indent2)
    dump_named_json_vector("D", D, to_file, indent=indent2)
    dump_named_json_vector("f", np.array(f) * fscale, to_file, comma=False, indent=indent2)
    to_file.write(indent + "}\n")

def all_categories_json(category_volume_array, to_json_path="all_categories.json", 
        blur=True, fscale=20, strides=2, special_categories=special_categories):
    outfile = open(to_json_path, "w")
    outfile.write("[\n")
    inside = False
    center = None
    for category in range(1, 15):
        print ('writing category', category)
        if inside:
            outfile.write(",\n")
        #outfile.write('"%s": ' % category)
        category_array = np.where(category_volume_array==category, 1, 0)
        if blur:
            category_array = blur_category(category_array)
        (A, B, C, D, f, center) = nonzero_tetrahedra(category_array, strides=strides)
        dump_tetrahedral_data(outfile, "", A, B, C, D, f, center, fscale=fscale, name=category,
            special_categories=special_categories)
        inside = True
    outfile.write("]\n")
    print("wrote ", to_json_path)
    return center

def one_category_html(category_volume_array, category, to_html_path=None, to_json_path=None, 
        hexcolor="0xffff00", fscale=20, blur=True):
    if to_html_path is None:
        to_html_path = "category_" + str(category) + ".html"
    if to_json_path is None:
        to_json_path = "category_" + str(category) + ".json"
    category_array = np.where(category_volume_array==category, 1, 0)
    if blur:
        category_array = blur_category(category_array)
    (A, B, C, D, f, center) = nonzero_tetrahedra(category_array)
    to_json = open(to_json_path, "w")
    dump_tetrahedral_data(to_json, "", A, B, C, D, f, fscale=fscale)
    print ("wrote", to_json_path)
    (i,j,k) = center  # map(int, center)
    print ("center", center)
    html = open("template.html").read()
    html = html.replace("FILE_PATH", to_json_path)
    html = html.replace("HEXCOLOR", hexcolor)
    look_at = repr((i,j,k))[1:-1]
    html = html.replace("CAMERA_LOOK_AT", look_at)
    orbit_center = look_at
    if not (i or j or k):
        orbit_center = "2,2,2"
    html = html.replace("ORBIT_CENTER", orbit_center)
    cutoff = fscale * 0.49
    html = html.replace("CUTOFF", str(cutoff))
    html = html.replace("CAMERA_X", str(-2*i))
    html = html.replace("CAMERA_Y", str(-2*j))
    html = html.replace("CAMERA_Z", str(-2*k))
    open(to_html_path, "w").write(html)
    print((to_json_path, hexcolor, orbit_center, look_at))
    print("wrote", to_html_path)

def read_kidney_slices(path="../"):
    "Read rachels data files into a single array."
    L = []
    for i in range(285):
        S = []
        fn = path + ("c%04d.txt" % i)
        f = open(fn)
        for line in f:
            r = list(map(int, line.strip().split()))
            S.append(r)
        L.append(S)
    # trim off the first element because it doesn't have the right length
    L = L[1:]
    return np.array(L, dtype="int")