import numpy as np
import pandas as pd
import cv2
import os
import pyOpenLPT as lpt

def print_dict(dict):
    for key in dict:
        print(key+":", dict[key])
    print('')

    
def load_csv(directory, keywords):
    # Get list of files in the directory
    files = os.listdir(directory)
    
    # Filter files that are CSV and contain any of the keywords
    matched_files = [f for f in files if f.endswith('.csv') and any(keyword in f for keyword in keywords)]
    
    # Load matched CSV files into DataFrames
    df_list = {}
    for file in matched_files:
        file_path = os.path.join(directory, file)
        df_list[file] = pd.read_csv(file_path)
    
    return df_list

def load_tracks(folder, load_previous=False, keywords=['ExitTrack', 'LongTrackActive', 'LongTrackInactive']):
    tracks = pd.DataFrame()
    
    # load latest tracks 
    df_list = load_csv(folder, keywords)
    
    for key in df_list:
        df = df_list[key]
        if not tracks.empty:
            df['TrackID'] += tracks['TrackID'].unique().shape[0]
        tracks = pd.concat([tracks, df], axis=0, ignore_index=True)
        
    if load_previous:
        # load previous converge tracks
        subfolder = os.path.join(folder, 'ConvergeTrack')
        df_list = load_csv(subfolder, ['ExitTrack', 'LongTrackInactive'])
        
        for key in df_list:
            df = df_list[key]
            if not tracks.empty:
                df['TrackID'] += tracks['TrackID'].unique().shape[0]
            tracks = pd.concat([tracks, df], axis=0, ignore_index=True)
    
    return tracks


def load_txt(directory, keywords):
    # Get list of files in the directory
    files = os.listdir(directory)
    
    # Filter files that are CSV and contain any of the keywords
    matched_files = [f for f in files if f.endswith('.txt') and any(keyword in f for keyword in keywords)]
    
    # Load matched CSV files into DataFrames
    df_list = {}
    for file in matched_files:
        file_path = os.path.join(directory, file)
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if len(lines) == 0:
                continue
            lines[0] = lines[0][:-1]
            
        data = np.loadtxt(lines, delimiter=',').reshape(-1, 5)
        df = pd.DataFrame(data, columns=['TrackID', 'FrameID', 'WorldX', 'WorldY', 'WorldZ'])
        df_list[file] = df
    
    return df_list

def load_tracks_old(folder, load_previous=False, keywords=['ExitTrack', 'ActiveLongTracks', 'InactiveLongTracks']):
    tracks = pd.DataFrame()
    
    # load latest tracks
    df_list = load_txt(folder, keywords)
    
    for key in df_list:
        df = df_list[key]
        if not tracks.empty:
            df['TrackID'] += tracks['TrackID'].unique().shape[0]
        tracks = pd.concat([tracks, df], axis=0, ignore_index=True)
        
    if load_previous:
        # load previous converge tracks
        subfolder = os.path.join(folder, 'ConvergeTrack')
        df_list = load_txt(subfolder, ['ExitTrack', 'InactiveLongTracks'])
        
        for key in df_list:
            df = df_list[key]
            if not tracks.empty:
                df['TrackID'] += tracks['TrackID'].unique().shape[0]
            tracks = pd.concat([tracks, df], axis=0, ignore_index=True)
    
    tracks['TrackID'] = tracks['TrackID'].astype(int)
    tracks['FrameID'] = tracks['FrameID'].astype(int) - 1
    
    return tracks
    

def getCamParam(file):
    camParam = {
        'Noffh': [],
        'Noffw': [],
        'R': [],
        'T': [],
        'Rinv': [],
        'Tinv': [],
        'f_eff': [],
        'hpix': [],
        'wpix': [],
        'Npixh': [],
        'Npixw': [],
        'k1': [],
        'k1star': [],
        'p1': [],
        'p1star': [],
        'p2': [],
        'p2star': [],
        'err_x': [],
        'err_y': [],
        'err_t': []
    }
    
    content = []
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(('#', '\n')):
                continue
            else:
                content.append(float(line.strip().split('#')[0]))
    
    ncam = int(content[0])
    camParam['ncam'] = ncam
    
    id = 1
    for _ in range(ncam):
        camParam['Noffh'].append(content[id])
        camParam['Noffw'].append(content[id+1])
        camParam['Npixw'].append(content[id+2])
        camParam['Npixh'].append(content[id+3])
        camParam['wpix'].append(content[id+4])
        camParam['hpix'].append(content[id+5])
        camParam['f_eff'].append(content[id+6])
        camParam['k1'].append(content[id+7])
        camParam['k1star'].append(0)
        camParam['R'].append(np.array(content[id+9:id+17+1]).reshape(3,3))
        camParam['T'].append(np.array(content[id+18:id+20+1]).reshape(3,1))
        camParam['Rinv'].append(np.array(content[id+21:id+29+1]).reshape(3,3))
        camParam['Tinv'].append(np.array(content[id+30:id+32+1]).reshape(3,1))
        camParam['p1'].append(0)
        camParam['p1star'].append(0)
        camParam['p2'].append(0)
        camParam['p2star'].append(0)
        camParam['err_x'].append(np.nan)
        camParam['err_y'].append(np.nan)
        camParam['err_t'].append(np.nan)

        id += 33
    return camParam

def project_old(tracks, camfile):
    camParam = getCamParam(camfile)
    ncam = camParam['ncam']
    
    xyz = tracks[['WorldX', 'WorldY', 'WorldZ']].values
    
    for i in range(ncam):
        temp = xyz @ camParam['R'][i].T + camParam['T'][i].T
        temp = temp / temp[:,2].reshape(-1,1) * camParam['f_eff'][i]
        pt2d = temp[:,0:2]
        pt2d = pt2d / (1+camParam['k1'][i]*(pt2d[:,0]**2 + pt2d[:,1]**2)).reshape(-1,1)
        
        pt2d[:,0] = pt2d[:,0] / camParam['wpix'][i] + camParam['Noffw'][i] + camParam['Npixw'][i]/2
        pt2d[:,1] = - pt2d[:,1] / camParam['hpix'][i] - camParam['Noffh'][i] + camParam['Npixh'][i]/2

        tracks['cam'+str(i)+'_x(col)'] = pt2d[:,0]
        tracks['cam'+str(i)+'_y(row)'] = pt2d[:,1]
        
    return tracks


def project_new(pt3d, folder, suffix=''):
    # load cam file 
    pt2d_list = []
    for i in range(4):
        file = os.path.join(folder, 'cam'+str(i+1)+suffix+'.txt')
        cam = lpt.math.Camera(file)
        pt2d = np.vstack([lpt.math.matrix_to_numpy(cam.project(lpt.math.Pt3D(pt3d[j,0], pt3d[j,1], pt3d[j,2]))).reshape(2,) for j in range(pt3d.shape[0])])
        pt2d_list.append(pt2d)
    pt2d_list = np.hstack(pt2d_list)
    return pt2d_list


def cal_img_combine(n_row, n_col, folder, name, frame_start, frame_end, dtype=np.uint8):
    img_combine = np.zeros((n_row, n_col)).astype(dtype)
    for frame_id in range(frame_start,frame_end+1):
        # name = 'img'+'{:05d}'+'.tif'
        file = os.path.join(folder, name.format(frame_id))
        img = cv2.imread(file, cv2.IMREAD_GRAYSCALE)
        judge = img > img_combine
        img_combine[judge] = img[judge]
    return img_combine


def get_lpt_tracks(tracks, n_init=4):
    tracks_lpt = []
    
    tracksID = tracks['TrackID'].unique()
    for trackID in tracksID:
        track = tracks[tracks['TrackID'] == trackID]
        
        if track.shape[0] < n_init:
            continue
        
        track_lpt = lpt.stb.TracerTrack()
        
        for i in range(0,track.shape[0]):
            pt3d = lpt.math.Pt3D(track['WorldX'].values[i], track['WorldY'].values[i], track['WorldZ'].values[i])
            tr3d = lpt.object.Tracer3D(pt3d)
            track_lpt.addNext(tr3d, track['FrameID'].values[i])
        
        tracks_lpt.append(track_lpt)
    
    return tracks_lpt