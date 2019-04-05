import glob
import os
import numpy as np
import json

DEFAULT_STRIDES = 8

def read_file_array(path):
    f = open(path)
    lines = f.readlines()
    flines = [list(map(float, line.split())) for line in lines]
    nrows = len(flines)
    [ncols] = list(set(len(nums) for nums in flines))
    array = np.array(flines, dtype=np.float)
    return array

def read_all_files(index_to_file):
    index_to_arrays = {}
    for (index, path) in index_to_file.items():
        index_to_arrays[index] = read_file_array(path)
    [(dim2, dim3)] = list(set(A.shape for A in index_to_arrays.values()))
    dim1 = len(index_to_file)
    result = np.zeros((dim1, dim2, dim3), dtype=np.float)
    for index in range(dim1):
        result[index] = index_to_arrays[index]
    return result

def get_protein_to_directories(parent="all_proteins"):
    bsuffix = "_bkgrd"
    lbsuffix = len(bsuffix)
    result = {}
    for subdir in os.listdir(parent):
        if subdir.endswith(bsuffix):
            protein = subdir[:-lbsuffix]
            result[protein] = parent + "/" + subdir
    return result

def protein_index_to_file(protein_directory):
    prefix = protein_directory + "/c"
    suffix = ".txt"
    lp = len(prefix)
    ls = len(suffix)
    glob_pattern = prefix + "*" + suffix
    files = list(sorted(glob.glob(glob_pattern)))
    index_to_file = {}
    assert files, "no files " + repr(glob_pattern)
    for path in files:
        indexstr = path[lp:-ls]
        index = int(indexstr)
        index_to_file[index] = path
    return index_to_file

def proteins_to_json(source="all_proteins", destination="../docs/data", rebin=DEFAULT_STRIDES):
    protein_to_json_path = {}
    p2d = get_protein_to_directories(source)
    for (protein, protein_directory) in p2d.items():
        index_to_file = protein_index_to_file(protein_directory)
        full_array = read_all_files(index_to_file)
        print ("...processing protein", protein, full_array.shape)
        rebin_array = average_down3(full_array, rebin)
        print ("... rebin shape", rebin_array.shape)
        s = np.array(rebin_array.shape)
        center = s * 0.5
        camera = s * 2.0
        array_json = json3darray(rebin_array)
        json_path = destination + "/protein_" + protein + ".json"
        out = open(json_path, "w")
        w = out.write
        w("{\n")
        w('"name": "' + protein + '",\n')
        w('"center": ' + str(list(center)) + ",\n")
        w('"camera": ' + str(list(camera)) + ",\n")
        w('"array": ' + array_json + "\n")
        w("}")
        out.close()
        protein_to_json_path[protein] = json_path
    list_file_path = destination + "/all_proteins.json"
    f = open(list_file_path, "w")
    json.dump(protein_to_json_path, f)
    print ("wrote", list_file_path)
    return protein_to_json_path

# https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array

def rebin_avg3(a, shape):
    (am, an, ap) = a.shape
    (sm, sn, sp) = shape
    sh = sm, am//sm, sn, an//sn, sp, ap//sp
    #return a.reshape(sh)
    return a.reshape(sh).mean(5).mean(3).mean(1)

def average_down3(a, factor):
    ashape = np.array(a.shape, dtype=np.int)
    truncated = ashape - (ashape % factor)
    [t1, t2, t3] = truncated
    a_truncated = a[:t1, :t2, :t3]
    fshape = truncated // factor
    return rebin_avg3(a_truncated, fshape)

def json3darray(A, fmt="%2.2f"):
    L = []
    ap = L.append
    ap("[")
    inside_sheet = False
    for sheet in A:
        inside_line = False
        if inside_sheet:
            ap(",\n\n ")
        ap("[")
        for line in sheet:
            if inside_line:
                ap(",\n  ")
            ap("[")
            inside_value = False
            for value in line:
                if inside_value:
                    ap(",")
                ap(fmt % value)
                inside_value = True
            ap("]")
            inside_line = True
        ap("]")
        inside_sheet = True
    ap("]")
    return "".join(L)


def test():
    proteins_to_json()

if __name__ == "__main__":
    test()
