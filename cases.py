import toml
import pandas as pd
import os
import shutil
import subprocess
from matplotlib import pyplot as plt
from pathlib import Path
import matplotlib.animation as animation
import numpy as np
import re
import h5py


def parse_name(name):
    return float(re.search('=([\d.]+)',name).group(1))

class BaseCase:
    def __init__(self,index,params=None):
        self.index=index

    def get_folder_name(self):
        return os.path.join(self.path,str(self.index))
class Case(BaseCase):
    params_name =['rel',
                  'time',
                  'pr',
                  'gr_v',
                  'omega',
                  'rel_v',
                  'le',
                  'sor',
                  'bm',
                  'l_sed'
                 ]

    params_c_name =['le','sor','bm','l_sed']
    run_name = './hele-shaw'
    path = './'
    def __init__(self,index,params=None,path='./'):
        self.index=index
        self.path=path
        if params:
            self._params = params
        else:
            self._params = self.read_params_from_df()
    def read_params_from_df(self):
        df = pd.read_csv(os.path.join(self.path,'db.csv'),index_col=0)
        parmas_row = df.iloc[self.index,:]
        return {k:float(parmas_row[k]) for k in self.params_name}

    def get_params(self):
        return self._params
    def get_config_params(self):
        conf = {k:v for k,v in self._params.items() if not( k in self.params_c_name)}
        conf['params_c'] = {k:v for k,v in self._params.items() if ( k in self.params_c_name)}
        return conf
    def save_config(self):
        config = self.get_config_params()
        file_name = os.path.join(self.get_folder_name(),'config.toml')

        with open(file_name,'w') as f:
            toml.dump(config,f)
    def create_folder(self):
        path = self.get_folder_name()
        os.makedirs(path,exist_ok=True)
        self.save_config()
        shutil.copy(self.run_name,path)
    def run(self):
        path =  self.get_folder_name()
        subprocess.run([self.run_name],cwd=path)
    def get_df(self):
        path =  os.path.join(self.get_folder_name(),'foo.csv')
        return pd.read_csv(path)
    def get_energy(self):
        paths =  os.path.join(self.get_folder_name(),'storag.h5')

        f = h5py.File(paths,'r')
        paths=list(f['map'].keys())
        paths.sort(key=parse_name)
        res = 0
        for stage in paths[-4000:]:
            try:
                df = f['map'][stage]
                z = np.array(df['psi']).T
                res+=(z*z).mean()
            except:
                pass

        return res/len(paths[-4000:])
    def get_energy2(self):
        paths =  os.path.join(self.get_folder_name(),'storag.h5')

        f = h5py.File(paths,'r')
        paths=list(f['map'].keys())
        paths.sort(key=parse_name)
        res = 0
        for stage in paths[-4000:]:
            try:
                df = f['map'][stage]
                z = np.array(df['psi']).T
                e=(np.diff(z,axis=0)**2)[:-1,1:-1]+(np.diff(z,axis=1)**2)[1:-1,:-1]/(20/40)**2
                res+=np.sum(e)/(e.shape[0]*e.shape[1])
            except:
                pass
        return res/len(paths[-4000:])



    def get_max_psi(self):
        df = self.get_df()
        df = df.dropna()
        df['psi_m']=df['psi_m'].astype(float)
        return df['psi_m'][len(df)//6*4:].mean()
    def get_max_nu(self):
        df = self.get_df()
        df = df.dropna()
        df['nu']=df['nu'].astype(float)
        return df['nu'][len(df)//6*4:].mean()

    def get_animation(self):

        paths = sorted(Path(os.path.join(self.get_folder_name(),'res')).iterdir(), key=os.path.getmtime)
        return get_animation(paths)

class CounturControl(BaseCase):
    field_name = ['psi','T','C']
    _stage = None
    _shape = None
    _meshgrid =None
    path = './'
    def __init__(self,index,params=None):
        self.index=index
    def get_stage(self,index=-1):
        if self._stage is None:
            paths = sorted(Path(os.path.join(self.get_folder_name(),'res')).iterdir(), key=os.path.getmtime)
            self._stage = pd.read_csv(paths[index])
        return self._stage

    def get_meshgrid(self):
        if self._meshgrid is None:
            df= self.get_stage()
            x = df.x.unique()
            y = df.y.unique()
            self._meshgrid = np.meshgrid(x, y)
        return self._meshgrid


    def get_shape(self):
        if self._shape is None:
            df= self.get_stage()
            x = df.x.unique()
            y = df.y.unique()
            NX= len(df.x.unique())
            NY= len(df.y.unique())
            self._shape= (NY,NX)
        return self._shape

    def get_countur_field(self,ax,name):
        shape = self.get_shape()
        df = self.get_stage()
        z = df[name].to_numpy().reshape(shape)

        X,Y = self.get_meshgrid()
        cs = ax.contourf(X, Y, z)

    def get_countur(self):
        fig,ax = plt.subplots(1,3,figsize=(9,3))
        for i,x in enumerate(self.field_name):
            self.get_countur_field(ax[i],x)
        return fig,ax



def get_animation(paths):
    fig,ax = plt.subplots(1,3,figsize=(9,3))
    paths=paths[:]
    # Method to change the contour data points
    def animate(i):
        stage= paths[i]
        df = pd.read_csv(stage)
        NX= len(df.x.unique())
        NY= len(df.y.unique())

        x = df.x.unique()
        y = df.y.unique()
        ax[0].clear()
        ax[1].clear()
        ax[2].clear()

        X, Y = np.meshgrid(x, y)
        z = df['psi'].to_numpy().reshape((NY, NX))
        cs = ax[0].contourf(X, Y, z,cmap=plt.cm.rainbow)
        #z = df['phi'].to_numpy().reshape((NY, NX))
    # cs = ax[1].contourf(X, Y, z,cmap=plt.cm.rainbow)
        ax[1].set_title(stage)
        z = df['T'].to_numpy().reshape((NY, NX))
        cs = ax[1].contourf(X, Y, z,cmap=plt.cm.rainbow,levels=100)
        z = df['C'].to_numpy().reshape((NY, NX))
        cs = ax[2].contourf(X, Y, z,cmap=plt.cm.rainbow,levels=100)
    # Call animate method
    ani = animation.FuncAnimation(fig, animate, len(paths), interval=100, blit=False)
    return ani

# Display the plot
plt.show()
