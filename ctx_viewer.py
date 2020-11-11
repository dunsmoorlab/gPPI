import os
import glob
import cortex
import numpy as np
import nibabel as nib
from preprocess_library import *
from fc_config import *
from glob import glob
import time

#### consider applying these things
# Subject.view.controls
# Subject.view.animate
# Subject.view.setData

class View(object):

    def __init__(self,subj=None,xfm_name='rai2rfi'):

        #to use fsaverage just type 'fs'
        # self.subj = meta(subj)

        #find the correct transform in the pycortex subject folder
        if self.subj.fsub == 'fsaverage': self.xfm_name = 'mni2fsavg'
        else: self.xfm_name = 'rai2rfi'

        self.pycx_id = self.subj.fsub

        proj_dir = data_dir
        self.opdir = os.path.join('graphing','op_decode','importance_maps') + os.sep


        if self.pycx_id not in cortex.db.subjects.keys():
            print('Subject {:s} does NOT have a structure in pycortex.'.format(self.subj.fsub))
            print('Run View.create_pycx_subj() to do this.') 

        # list of images that can be shown in viewer
        # self.avail_imgs = [ os.path.basename(f).split('.')[0]
        #                   for f in glob(os.path.join(self.ref_dir,'*.nii.gz')) ]


    def _get_nii_data(self,nii_name=None):
        if nii_name is None:
            if self.subj.fsub == 'gusbrain': data = np.zeros([36,100,100])
            elif self.subj.fsub == 'fsaverage': data = np.zeros([182,218,182])
            else: data = np.zeros([48,76,76])
        else:
            # fname = os.path.join(self.ref_dir,'{:s}.nii.gz'.format(nii_name))
            # if not os.path.exists(fname):
            #   raise IOError(('fname {:s} not found!'.format(fname)
            #       +'\nPick from these available images:\n{:s}'.format(
            #           ','.join([ x for x in self.available_imgs]))))
            # data = nib.load(fname).get_data().swapaxes(0,2) #,'float32')
            data = os.path.join(data_dir,nii_name)
            # data = os.path.join(data_dir,nii_name)
        return data


    def show(self,nii=None,**kwargs):
        # get data from nii image
        data = self._get_nii_data(nii_name=nii)
        print(data)
        # create pycortex volume structure
        #####use vmin/vmax kwrags######
        self.vol = cortex.Volume(data,self.pycx_id,self.xfm_name,
                                    #           with_curvature=True,
                                    # curvature_contrast=0.65,
         #                            curvature_brightness=0.16,
         #                            curvature_threshold=True,
                                    **kwargs)
        
        # open in browser
        self.view = cortex.webshow(self.vol,autoclose=False,
                                    # with_curvature=True,
                                    # curvature_contrast=0.65,
         #                            curvature_brightness=0.16,
         #                            curvature_threshold=True,
                                    **kwargs)
        # view = cortex.webs`how(vol, port=port, autoclose=False, 
        #                     template=html_template, title='instabrain')
        self.view.animate([{'state':'surface.%s.specularity'%(self.subj.fsub),'idx':0,'value':0}])
        if self.subj.fsub == 'fsaverage':
            time.sleep(3)
            self.view.animate([{'state':'mix','idx':.5,'value':.5}])

        if nii is None or 'mask' in nii:
            self.view.setVminmax(-1,1)


    def update_show(self,nii):
        data = os.path.join(data_dir,nii)
        vol = cortex.Volume(data,self.pycx_id,self.xfm_name)
        mos, _ = cortex.mosaic(vol.data)
        # mos, _ = cortex.mosaic(vol.volume[0],show=False)
        self.view.dataviews.data.data[0]._setData(0,mos)
        # self.view.dataviews.data.data[0]._setData(1,mos)
        # if init:
        #     self.view.setVminmax(np.percentile(data,75),np.percentile(data,99))
        self.view.setFrame(1)

    def reset(self):

        self.view.animate([{'state':'camera.azimuth','idx':0,'value':45},
                            {'state':'camera.altitude','idx':0,'value':75},
                            {'state':'camera.radius','idx':0,'value':400},
                            {'state':'surface.%s.pivot'%(self.subj.fsub),'idx':0,'value':0}])

    def move(self,clip=None):
        # self.view.animate([{'state':'camera.altitude','idx':0,'value':90}])
        # time.sleep(2)
        if clip == 'mask':
            self.view.animate([{'state':'camera.azimuth','idx':2,'value':180}])
            time.sleep(2.1)
            self.view.animate([{'state':'surface.%s.pivot'%(self.subj.fsub),'idx':2,'value':-99}])
        elif clip == 'CS+':
            self.view.animate([{'state':'camera.azimuth','idx':0,'value':0},
                    {'state':'camera.altitude','idx':0,'value':90},
                    {'state':'camera.radius','idx':0,'value':400},
                    {'state':'surface.%s.pivot'%(self.subj.fsub),'idx':0,'value':-180}])

            self.view.animate([{'state':'camera.azimuth','idx':0,'value':0},
                    {'state':'camera.altitude','idx':0,'value':90},
                    {'state':'camera.radius','idx':2,'value':200},
                    {'state':'surface.%s.pivot'%(self.subj.fsub),'idx':2,'value':180}])

        elif clip == 'corr':
            self.view.animate([{'state':'camera.azimuth','idx':0,'value':0},
                    {'state':'camera.altitude','idx':0,'value':90},
                    {'state':'camera.radius','idx':0,'value':400},
                    {'state':'surface.%s.pivot'%(self.subj.fsub),'idx':0,'value':-180}])

            time.sleep(3)
            self.view.animate([{'state':'camera.azimuth','idx':0,'value':0},
                    {'state':'camera.altitude','idx':0,'value':90},
                    {'state':'camera.radius','idx':2,'value':200},
                    {'state':'surface.%s.pivot'%(self.subj.fsub),'idx':2,'value':180}])

            self.view.animate([{'state':'camera.azimuth','idx':0,'value':0},
                    {'state':'camera.altitude','idx':0,'value':90},
                    {'state':'camera.radius','idx':2,'value':200},
                    {'state':'surface.%s.pivot'%(self.subj.fsub),'idx':2,'value':180}])

        elif clip == 'spin':
            self.view.animate([{'state':'camera.azimuth','idx':3,'value':180},
                                #{'state':'camera.azimuth','idx':1.5,'value':360}
                                ])
            time.sleep(.01)
            self.view.animate([{'state':'camera.azimuth','idx':3,'value':360}])

    def create_pycx_subj(self):
        
        if not os.path.exists(self.subj.fs_regmat):
            reg_cmd = 'bbregister --s %s/%sfs --mov %s --init-fsl --bold --reg %s'%(self.subj.fsub, self.subj.fsub, self.subj.refvol_be, self.subj.fs_regmat)
            os.system(reg_cmd)

        cortex.freesurfer.import_subj(freesurfer_subject_dir=self.subj.subj_dir,
                                        subject=self.subj.fs_id,sname=self.pycx_id)
        # NOTE there is also a cortex.freesurfer.import_flat()

        # save a transform into the pycortex subject's dir
        cortex.xfm.Transform.from_freesurfer(
                freesurfer_subject_dir=self.subj.subj_dir,
                fs_register=self.subj.fs_regmat,
                func_nii=self.subj.refvol,
                subject=self.subj.fs_id,
        ).save(subject=self.pycx_id,
               name=self.xfm_name, # places into subj/transforms/xfm_name
               xfmtype='coord')

        print('Created pycortex subject {:s} successfully'.format(self.pycx_id))
        print('Structure placed in {:s}'.format(cortex.database.default_filestore))


        print('Reload pycortex before viewing subject.')