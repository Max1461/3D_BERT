from pymol.cgo import *
from pymol import cmd
import numpy
import h5py
import os
import itertools

# You have to play with the size value, should be roughly as bog as a gridvoxel
def addSphere(color, feature_vector, size=.5):
    color = [float(i)for i in color]
    return [
       COLOR,	color[0],	color[1],   color[2],
       SPHERE,   feature_vector[0],   feature_vector[1],   feature_vector[2], size,
    ]

def heatmap_color(value, vmin, vmax):	
    blue = numpy.array([0,0,255])
    white = numpy.array([255,255,255])
    red = numpy.array([255,0,0])
    t = (value - vmin) / ((vmax - vmin) / 2)/2
    print(t)
    if t < .5:
        t = t*2
        color = blue*(1-t) + white*t
    else:
        t = (t-.5)*2
        color = red*t + white*(1-t)
    return color/255
    
    
def plotFeatures(feature_vectors, feature_name='very_important_feature'):
    # assumes feature_vector -> [ x, y, z, single_value ]
    values = [fv[-1] for fv in feature_vectors]
    min_value = min(values)
    max_value = max(values)
    colors = [heatmap_color(value, min_value, max_value)for value in values]
    spheres = []	
    for i in range(len(values)):
        spheres += addSphere(colors[i], feature_vectors[i])
    cmd.load_cgo(spheres, feature_name,   1)
    return feature_name
    
def draw_gridbox(nX, nY, nZ, lw=0.01, boxname='grid'):
    gridbox = [
    LINEWIDTH, float(lw),
    BEGIN, LINES,
    COLOR, 1.0, 1.0, 1.0,
    ]
    for i, x in enumerate(nX[:-1]):
        for j, y in enumerate(nY[:-1]):
            for k, z in enumerate(nZ[:-1]):
                dots= [
                    VERTEX, nX[i], nY[j], nZ[k],
                    VERTEX, nX[i], nY[j], nZ[k+1],
                    VERTEX, nX[i], nY[j+1], nZ[k],
                    VERTEX, nX[i], nY[j+1], nZ[k+1],
                    VERTEX, nX[i+1], nY[j], nZ[k],
                    VERTEX, nX[i+1], nY[j], nZ[k+1],
                    VERTEX, nX[i+1], nY[j+1], nZ[k],
                    VERTEX, nX[i+1], nY[j+1], nZ[k+1],
                    VERTEX, nX[i], nY[j], nZ[k],
                    VERTEX, nX[i+1], nY[j], nZ[k],
                    VERTEX, nX[i], nY[j+1], nZ[k],
                    VERTEX, nX[i+1], nY[j+1], nZ[k],
                    VERTEX, nX[i], nY[j+1], nZ[k+1],
                    VERTEX, nX[i+1], nY[j+1], nZ[k+1],
                    VERTEX, nX[i], nY[j], nZ[k+1],
                    VERTEX, nX[i+1], nY[j], nZ[k+1],
                    VERTEX, nX[i], nY[j], nZ[k],
                    VERTEX, nX[i], nY[j+1], nZ[k],
                    VERTEX, nX[i+1], nY[j], nZ[k],
                    VERTEX, nX[i+1], nY[j+1], nZ[k],
                    VERTEX, nX[i], nY[j], nZ[k+1],
                    VERTEX, nX[i], nY[j+1], nZ[k+1],
                    VERTEX, nX[i+1], nY[j], nZ[k+1],
                    VERTEX, nX[i+1], nY[j+1], nZ[k+1],
                ]
                gridbox += dots
    gridbox.append(END)

    cmd.load_cgo(gridbox,boxname)

def index2xyz(index, X, Y, Z):
    # indx=range(X)
    # indy=range(Y)
    # indz=range(Z)
    indexes = [X, Y, Z]
    points = list(itertools.product(*indexes))
    
    return points[index]

def plotMappedFeatures(f, case, feature_name='very_important_feature'):
    values = f[case]['mapped_features']['Feature_ind'][feature_name]['value'][:]
    try:
        min_value = min(values)
        max_value = max(values)
    except ValueError:
        min_value = -1.0
        max_value = 1.0
    colors = [heatmap_color(value, min_value, max_value)for value in values]
    spheres = []
    X = f[case]['grid_points']['x'][:]
    Y = f[case]['grid_points']['y'][:]
    Z = f[case]['grid_points']['z'][:]
    indexes = [X, Y, Z]
    points = list(itertools.product(*indexes))
    for i, index in enumerate(f[case]['mapped_features']['Feature_ind'][feature_name]['index'][:]):
        x, y, z = points[index]
        value = f[case]['mapped_features']['Feature_ind'][feature_name]['value'][i]
        spheres += addSphere(colors[i], [x, y, z, value], size=.1)
    cmd.load_cgo(spheres, f'mapped_{feature_name}', 1)
    return f'mapped_{feature_name}'
    
# def test(f, case):
#     #examples = [[0,0,0, 100], [0,10,0, 0] ,[0,10,10,50],[0,0,10,25]]
    
#     # for each mapped_features index
#     #f['1AK4_10w']['grid_points']['x'][index], f['1AK4_10w']['grid_points']['y'][index], f['1AK4_10w']['grid_points']['z'][index]
#     d = f[case]['features']['charge'][:][:,1:]
#     plotFeatures(d, 'dario_example')

def visualize_hdf5(filepath, case, features=True, mapped_features=False):
    """Draws features and grid from DeepRank HDF5 file

    Args:
        filepath (str): file to HDF5 file
        case (str): Case id used as key in the hdf5 file
        features (bool, optional): Set it to 0 to not draw the features. Defaults to True.
        mapped_features (bool, optional): Whether to draw or not mapped features. 
            It is very slow, so it is better to draw as little as possible. It has to be
            given as list of feature names. Defaults to False.
    """    
    #cmd.feedback('disable', 'all', 'everything')
    f = h5py.File(filepath)
    features_objects = []
    if features != 'False' and features !='0':
        for feature in f[case]['features'].keys():
            d = f[case]['features'][feature][:][:,1:]
            feature_obj = plotFeatures(d, f'{case}_{feature}')
            features_objects.append(feature_obj)
        
    if mapped_features == 'all':
        for feature in f[case]['mapped_features']['Feature_ind'].keys():
            feature_obj = plotMappedFeatures(f, case, feature)
            features_objects.append(feature_obj)
            
    elif type(mapped_features)==str and mapped_features !='all':
        mapped_features = mapped_features.replace('[',
                                                  '').replace(']',
                                                             '').replace(' ','').split(',')
        #if type(mapped_features)==str:
        #    mapped_features = list(mapped_features)
        for feature in f[case]['mapped_features']['Feature_ind'].keys():
            if any(x in feature for x in mapped_features):
                feature_obj = plotMappedFeatures(f, case, feature)
                features_objects.append(feature_obj)
            # else:
            #     cmd.feedback('enable', 'all', 'everything')
            #     print(f'WARNING: mapped feature {feature} not found')
            #     cmd.feedback('disable', 'all', 'everything')
    else:
        pass
    
    x = f[case]['grid_points']['x'][:]
    y = f[case]['grid_points']['y'][:]
    z = f[case]['grid_points']['z'][:]
    boxname=f'{case}_grid'
    draw_gridbox(x, y, z, boxname=boxname)
   

    name = f'{case}.pdb'
    with open(name, 'w') as outfile:
        outfile.writelines([str(x).replace("b'","").replace("'","")+'\n' for x in f[case]['complex'][:]])
    cmd.load(name)
    
    cmd.group(name = f'{case}_features', members =(' ').join(features_objects + [ boxname, case]))
    
    os.system(f'rm {name}')

#visualize_hdf5('/Users/Dario/Documents/DeepXplorer/example.hdf5', '1AK4_10w')
cmd.extend ("visualize_hdf5", visualize_hdf5)