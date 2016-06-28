	# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
albanPRF.py
lots of functions that have been changed for the weibull PRF experiment
"""

import os, sys
from tempfile import mkdtemp
import gc
#from sklearn.externals import joblib
from random import sample as randomsample
from pylab import *
import numpy as np
import scipy as sp
from scipy.stats import spearmanr,spearmanr
from scipy import ndimage
import seaborn as sns
import statsmodels
import socket
from matplotlib import animation
import hrf_estimation as he
import random as random

import multiprocessing as mp

import time as time_module

from nifti import *
from math import *
import shutil

from joblib import Parallel, delayed

from sklearn.linear_model import ARDRegression, BayesianRidge, Ridge, RidgeCV, ElasticNet, ElasticNetCV
from skimage.morphology import disk
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, report_errors


#from joblib import Parallel, delayed


from IPython import embed as shell

# import toolbox funcitonality
from Tools.Operators.ArrayOperator import *
from Tools.Operators.EyeOperator import *
from Tools.Operators.PhysioOperator import *
from Tools.Operators.CommandLineOperator import *
from Tools.Sessions import *
from Tools.Sessions.Session import *
from Tools.Subjects.Subject import *
from Tools.Run import *
from Tools.Projects.Project import *

from ModelExperiment import *
from FitPRFModel import *
from EyeFromfMRISession import *
from PopulationReceptiveFieldMappingSession import *

import pkg_resources
pkg_resources.require("joblib==0.8.4")
import joblib

class WeibullPopulationReceptiveFieldMappingSession(PopulationReceptiveFieldMappingSession):

	"""
	Class for population receptive field mapping sessions analysis focussed on Weibull image statistics
	10 conditions, 8 full directional trials and 4 empty ones per condition. 
	"""

	def __init__(self, ID, date, project, subject, this_project_folder, parallelize = True, loggingLevel = logging.DEBUG,**kwargs):
		super(WeibullPopulationReceptiveFieldMappingSession, self).__init__(ID, date, project, subject, parallelize = parallelize, 
			loggingLevel = loggingLevel,this_project_folder=this_project_folder)

		self.n_pixel_elements_raw = 101
		self.n_pixel_elements_convolved = 31
		self.nr_dummy_scans = 6
		self.this_project_folder = this_project_folder
		self.stim_duration_in_TR = 24
		self.saccade_threshold = 5.0
		self.orientations = [0,45,90,135,180,225,270,315]
		self.hrf_shape = 'doubleGamma'
		self.TR = 1.594
		for k,v in kwargs.items():
			setattr(self, k, v)
		# results_frames = {'real_polar':0,'polar':1,'delta_amplitude':2,'ecc':3,'xo':4,'yo':5,'surround_size':6,'amp_center':7,'imag_polar':8,'amp_surround':9,'sigma_surround':10,'SI':11,'delta_sigma':12,'sigma_center':13,'fwhm':14,'real_eccen':15,'imag_eccen':16}

		#self.task_names = {'PRF':['Color', 'Speed', 'Fix', 'Fix_no_stim'],'Mapper':['fix_no_stim','no_color_no_speed','yes_color_no_speed','no_color_yes_speed','yes_color_yes_speed']}
		self.task_names = {'PRF':['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']}

		self.n_TR = {'PRF':300,'Mapper':510}

		# wat doet die? self.n_TR? was 765; heb er 300 van gemaakt
		#self.run_types = ['PRF','Mapper']
		self.run_types = ['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'] #,'T2_anat','T1']


	def create_moco_check_gifs(self,conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'],postFix=['mcf'],fps=10):

		pl.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
		general_plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/'))
		if not os.path.isdir(general_plotdir): os.mkdir(general_plotdir)
		plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/moco/'))
		if not os.path.isdir(plotdir): os.mkdir(plotdir)	

		for condition in conditions:

			# creat single run animation with shuffled TRs
			for r in [self.runList[i] for i in self.conditionDict[condition]]:


				func_fn = self.runFile(stage = 'processed/mri', run = r, postFix=postFix )
				func_data = NiftiImage(func_fn).data
				
				self.logger.info('creating %s animation for run %s'%(condition,r.ID))
				ims = []
				f=pl.figure(figsize=(4,4))
				slice_dim = np.argmin(np.shape(func_data))
				timepoints = np.arange(len(func_data))
				timepoints = randomsample(timepoints,100)
				# np.random.shuffle(timepoints)
				for t in timepoints:
					s=f.add_subplot(111)
					if slice_dim == 1:
						im=pl.imshow(func_data[t,15,:,:],origin='lowerleft',interpolation='nearest',cmap='gray')
					elif slice_dim == 2:
						im=pl.imshow(func_data[t,:,15,:],origin='lowerleft',interpolation='nearest',cmap='gray')
					elif slice_dim == 3:
						im=pl.imshow(func_data[t,:,:,15],origin='lowerleft',interpolation='nearest',cmap='gray')
					pl.axis('off')
					ims.append([im])
				ani = animation.ArtistAnimation(f, ims)#, interval=5, blit = True, repeat_delay = 1000)
				mywriter = animation.FFMpegWriter(fps = fps)
				self.logger.info('saving to %s_%s_run_%s.mp4'%('_'.join(postFix),condition,r.ID))
				ani.save(os.path.join(plotdir,'%s_%s_run_%s.mp4'%('_'.join(postFix),condition,r.ID)),writer=mywriter,dpi=200,bitrate=200)#,fps=2)#,dpi=100,bitrate=50)
				pl.close()

	def create_B0_distortion_check_gifs(self,conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'],fps=2,n_repetitions=50,func_data_postFix=['mcf','meanvol','NB','flirted2sessionT2']):
	# def create_B0_distortion_check_gifs(self,conditions=['PRF','Mapper'],fps=2,n_repetitions=50,func_data_postFix=['mcf','meanvol','NB','flirted2sessionT2']):


		pl.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
		general_plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/'))
		if not os.path.isdir(general_plotdir): os.mkdir(general_plotdir)
		plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/B0distortion/'))
		if not os.path.isdir(plotdir): os.mkdir(plotdir)	
		sub_plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/B0distortion/%s'%'_'.join(func_data_postFix)))
		if not os.path.isdir(sub_plotdir): os.mkdir(sub_plotdir)	

		for condition in conditions:

			for r in [self.runList[i] for i in self.conditionDict[condition]]:

				target_T2_run_idx = np.where(np.array([self.runList[ri].ID for ri in range(len(self.runList))]) == r.thisSessionT2ID)[0][0]
				func_data = NiftiImage(self.runFile(stage = 'processed/mri', run = r, postFix = func_data_postFix )).data
				T2_data = NiftiImage(self.runFile( stage = 'processed/mri', run = self.runList[target_T2_run_idx],postFix=['NB'])).data

				for direction in range(3):

					for depth in [1,2]:

						f=pl.figure(figsize=(4,4))
						ims = []
						s=f.add_subplot(111)
						for repetition in range(n_repetitions):
							if direction ==0:
								im=pl.imshow(func_data[int(np.shape(func_data)[0]/3*depth),:,:],origin='lowerleft',interpolation='nearest',cmap='gray')
							elif direction ==1:
								im=pl.imshow(func_data[:,int(np.shape(func_data)[1]/3*depth),:],origin='lowerleft',interpolation='nearest',cmap='gray')
							elif direction ==2:
								im=pl.imshow(func_data[:,:,int(np.shape(func_data)[2]/3*depth)],origin='lowerleft',interpolation='nearest',cmap='gray')
							pl.axis('off')
							ims.append([im])

							if direction ==0:
								im=pl.imshow(T2_data[int(np.shape(T2_data)[0]/3*depth),:,:],origin='lowerleft',interpolation='nearest',cmap='gray')
							elif direction ==1:
								im=pl.imshow(T2_data[:,int(np.shape(T2_data)[1]/3*depth),:],origin='lowerleft',interpolation='nearest',cmap='gray')
							elif direction ==2:
								im=pl.imshow(T2_data[:,:,int(np.shape(T2_data)[2]/3*depth)],origin='lowerleft',interpolation='nearest',cmap='gray')
							pl.axis('off')
							ims.append([im])

						ani = animation.ArtistAnimation(f, ims)#, interval=5, blit = True, repeat_delay = 1000)
						mywriter = animation.FFMpegWriter(fps = 3)
						self.logger.info('saving to %s_run_%s_direction_%s_depth_%s.mp4'%(condition,r.ID,direction,depth))
						ani.save(os.path.join(sub_plotdir,'%s_run_%s_direction_%s_depth_%s.mp4'%(condition,r.ID,direction,depth)),writer=mywriter,dpi=200,bitrate=200)#,fps=2)#,dpi=100,bitrate=50)
						pl.close()

	def check_EPI_alignment(self,conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'],fps=10,n_repetitions=10,postFix=['mcf','meanvol','NB','flirted2targetEPI']):
	#def check_EPI_alignment(self,conditions=['PRF','Mapper'],fps=10,n_repetitions=10,postFix=['mcf','meanvol','NB','flirted2targetEPI']):

		pl.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
		general_plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/'))
		if not os.path.isdir(general_plotdir): os.mkdir(general_plotdir)
		plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/between_session_registration/'))
		if not os.path.isdir(plotdir): os.mkdir(plotdir)	
		sub_plotdir = os.path.join(self.stageFolder('processed/mri/figs/check_registration/between_session_registration/%s'%'_'.join(postFix)))
		if not os.path.isdir(sub_plotdir): os.mkdir(sub_plotdir)	

		# load epis
		all_epis = []
		self.logger.info('loading all niftis')
		#for r in [self.runList[i] for i in self.conditionDict[condition]]:
		for r in [self.runList[i] for i in self.scanTypeDict['epi_bold']]:

			filename = self.runFile(stage = 'processed/mri', run = r, postFix=postFix )
			all_epis.append(NiftiImage(filename).data)

		all_epis = np.array(all_epis)

		self.logger.info('creating animation')
		pl.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
		ims = []
		for direction in range(3):
			for depth in [1,2]:
				f=pl.figure(figsize=(4,4))
				ims = []
				s=f.add_subplot(111)
				for repetition in range(n_repetitions):
					for epir in range(len([self.runList[i] for i in self.scanTypeDict['epi_bold']])):
						if direction ==0:
							im=pl.imshow(all_epis[epir,int(np.shape(all_epis)[0]/3*depth),:,:],origin='lowerleft',interpolation='nearest',cmap='gray')
						elif direction ==1:
							im=pl.imshow(all_epis[epir,:,int(np.shape(all_epis)[1]/3*depth),:],origin='lowerleft',interpolation='nearest',cmap='gray')
						elif direction ==2:
							im=pl.imshow(all_epis[epir,:,:,int(np.shape(all_epis)[2]/3*depth)],origin='lowerleft',interpolation='nearest',cmap='gray')
						pl.axis('off')
						ims.append([im])

				ani = animation.ArtistAnimation(f, ims)#, interval=5, blit = True, repeat_delay = 1000)
				mywriter = animation.FFMpegWriter(fps = fps)
				self.logger.info('saving to direction_%s_depth_%s.mp4'%(direction,depth))
				ani.save(os.path.join(sub_plotdir,'direction_%s_depth_%s.mp4'%(direction,depth)),writer=mywriter,dpi=200,bitrate=300)#,fps=2)#,dpi=100,bitrate=50)
				pl.close()
					
	def fnirt_mean_moco_to_mean_target_EPI(self):

		fnirts = []

		for r in [self.runList[i] for i in self.scanTypeDict['epi_bold']]:

			# fnirt mean mcf volume to target epi
			target_run_idx = np.where(np.array([self.runList[ri].ID for ri in range(len(self.runList))]) == self.targetEPIID)[0][0]
			fnirt_target_fn = self.runFile( stage = 'processed/mri', run = self.runList[target_run_idx],postFix = ['mcf','meanvol'])
			input_fn = self.runFile(stage = 'processed/mri', run = r, postFix = ['mcf','meanvol'])
			output_fn = self.runFile(stage = 'processed/mri', run = r, postFix = ['mcf','meanvol','fnirted2targetEPI'])
			initial_flirt_mat_fn = self.runFile(stage = 'processed/mri', run = r, postFix = ['mcf','meanvol','NB','flirted2targetEPI'], extension='.mat')
			coefs_fn = self.runFile(stage = 'processed/mri', run = r, postFix = ['mcf','meanvol','fnirted2targetEPI','coefs'])

			flO = FnirtOperator(inputObject=input_fn, referenceFileName = fnirt_target_fn)
			flO.configure(AffineTransMatrixFileName=initial_flirt_mat_fn,outputFileName=output_fn,coefsFileName=coefs_fn)
			fnirts.append(flO)

		# now execute in parallel
		ppservers = ()
		job_server = pp.Server(ppservers=ppservers, secret='mc')
		self.logger.info("starting pp with", job_server.get_ncpus(), "workers for " + sys._getframe().f_code.co_name)
		ppResults = [job_server.submit(ExecCommandLine,(fnirt.runcmd,),(),('subprocess','tempfile',)) for fnirt in fnirts]
		for fMcf in ppResults:
			fMcf()
		job_server.print_stats()


	def applywarp_to_moco_data(self):

		AWOs = []

		# the fnirt includes the flirt matrix, so we can simply apply the warp
		# to the motion corrected data
		for r in [self.runList[i] for i in self.scanTypeDict['epi_bold']]:


			# get the rawest func file to apply to
			func_fn = self.runFile(stage = 'processed/mri', run = r,postFix = ['mcf'])

			# the fnirt warpfields:
			warpfield_fn = self.runFile(stage = 'processed/mri', run = r, postFix = ['mcf','meanvol','fnirted2targetEPI','coefs'])

			# then, we'll create the applywarp command
			ofn = self.runFile(stage = 'processed/mri', run = r, postFix=['mcf','fnirted'])
			AWO = ApplyWarpOperator(inputObject = func_fn, referenceFileName = func_fn )
			AWO.configure(outputFileName = ofn,warpfieldFileName=warpfield_fn	)
			AWOs.append(AWO)

		# execute the applywarp operations in parallel
		ppservers = ()
		job_server = pp.Server(ppservers=ppservers, secret='mc')
		self.logger.info("starting pp with", job_server.get_ncpus(), "workers for " + sys._getframe().f_code.co_name)
		ppResults = [job_server.submit(ExecCommandLine,(AWO.runcmd,),(),('subprocess','tempfile',)) for AWO in AWOs]
		for fMcf in ppResults:
			fMcf()
		job_server.print_stats()


	def design_matrices_for_averaged_data(self, gamma_hrfType = 'doubleGamma', hrf_params = [1,0,0]#gamma_hrfParameters = {'a1' : 6, 'a2' : 12, 'b1' : 0.9, 'b2' : 0.9, 'c' : 0.35}, 
				,sample_duration = 0.01, ssr = 30,animate_dm=True,
				orientations = [315,270,135,180,225,45,90,0]):
		"""
		This function creates a design matrix with one trial per orientation direction (and a fix_no_stim trial). 
		Connected together, these trials create a model experiment. The arrange data function below will select
		windows from the continuous data and interpolate and median them to match this design matrix.

		A trial in the experiment contained the following timing:
		0-1.5 seconds: instruction (1.5 s)
		1.5-22.5 seconds: bar sweep (21 s)
		22.5-~26 seconds: ITI (3.5 s)

		A trial in the dm is modelled as:
		0-21 seconds: bar_sweep (21 s)
		21-end seconds: ITI (whatever left in TRs)
		"""

		for convolve in [False,True]:

			if convolve:
				n_pixel_elements =  31
			else:
				n_pixel_elements = 101

			self.logger.info('Creating design matrix with convolution is %s'%convolve)

			# create design matrix at higher time resolution to avoid aliasing artifacts in the creation process. 
			# note that the ssr should be the amount of slices
			sample_duration = self.TR/ssr

			# now create a model run for an 'experiment' in which every direction only occured once. 
			# Some variables need to be added to a (new) run object, as the PRFmodelRun demands it in that format.
			class empty_object(object): pass 
			r = empty_object()
			after_period = 3 * self.TR # period to include in dm after stim offset
			stim_duration = self.stim_duration_in_TR * self.TR

			trial_duration = stim_duration + after_period # in chronological order
			TRs_in_trial = int(round(trial_duration/self.TR)) # take floor to be sure not to add TRs contaminated by instruction response
			trial_types = np.hstack([np.tile('PRF_01',8),'Fix_no_stim'])
			r.orientations = np.radians(np.concatenate([orientations,[0]]))  # add extra 0 orientation for fix_no_stim trial
			stim_start_times = [(TRs_in_trial*self.TR*t) for t in np.arange(len(r.orientations))]
			stim_end_times = np.array(stim_start_times) + stim_duration 
			r.trial_times = [[trial_types[t],stim_start_times[t],stim_end_times[t]] for t in range(len(stim_start_times))]
			n_TRs = TRs_in_trial * len(trial_types)
			sample_times = np.linspace(0,TRs_in_trial*self.TR,TRs_in_trial,endpoint=False) # include this for only one trial, for later data selection
			n_samples = n_TRs*ssr

			bar_width = 4.0/9.0 # the screen is 2 wide, so 0.25 is bar width of 0.125*screen diameter
			#for weibullPRF the bar is 2/9th of the width. x2 = 4/9
			mr = PRFModelRun(r, n_samples=n_samples, n_pixel_elements=n_pixel_elements,sample_duration = sample_duration, bar_width = bar_width)
			self.logger.info('simulating model experiment with %d pixel elements and %1.2f s sample_duration'%(n_pixel_elements, self.TR))
			mr.simulate_run_original()

			if not convolve:
				down_sample_index = np.arange(0,np.shape(mr.run_matrix)[0],ssr).astype(int) # this downsamples from first sample in steps of ssr
				raw_design_matrix = mr.run_matrix[down_sample_index]
				self.logger.info('saving design matrix')
				with open(os.path.join(self.stageFolder('processed/mri/PRF/'), 'raw_design_matrix_for_average_data.pickle'), 'w') as f:
					pickle.dump( {'raw_design_matrix':raw_design_matrix,'sample_times':sample_times} , f)
				#np.save(os.path.join(self.stageFolder('processed/mri/PRF/'), 'raw_design_matrix_for_average_data.npy'),raw_design_matrix)

			else:
				self.logger.info('convolving design matrix with hrf')
				# add empty spaces between trials to avoid hrf lag contamination
				samples_in_trial = TRs_in_trial*ssr
				add_empty_samples = 20*ssr
				padded_dm = np.zeros(((samples_in_trial+add_empty_samples)*len(trial_types),n_pixel_elements,n_pixel_elements))
				for triali in range(len(trial_types)):
					padded_dm[(samples_in_trial+add_empty_samples)*triali:(samples_in_trial+add_empty_samples)*triali + samples_in_trial] = mr.run_matrix[samples_in_trial*triali:samples_in_trial*(triali+1)]
				run_design = Design(padded_dm.shape[0], mr.sample_duration, subSamplingRatio = 1) # ssr can be one, as the mr.run_matrix is still upsampled
				rdm = padded_dm.reshape((padded_dm.shape[0], padded_dm.shape[1] * padded_dm.shape[2])).T
				run_design.rawDesignMatrix = np.repeat(rdm, 1, axis=1) # ssr can be one, as the raw dm is already upsampled
				# run_design.convolveWithHRF(hrfType = gamma_hrfType, hrfParameters = gamma_hrfParameters)
				run_design.convolveWithHRF(hrfParameters = hrf_params)
				convolved_dm = np.reshape(np.swapaxes(run_design.designMatrix,0,1),(-1,n_pixel_elements,n_pixel_elements))
				# remove empty space again:
				convolved_design_matrix = np.zeros((samples_in_trial*len(trial_types),n_pixel_elements,n_pixel_elements))
				for triali in range(len(trial_types)):
					convolved_design_matrix[samples_in_trial*triali:samples_in_trial*(triali+1)] = convolved_dm[(samples_in_trial+add_empty_samples)*triali:(samples_in_trial+add_empty_samples)*triali + samples_in_trial]
				# we'll not downsample the convolved dm so we can later select slices
				self.logger.info('saving design matrix')
				with open(os.path.join(self.stageFolder('processed/mri/PRF/'), 'convolved_design_matrix_for_average_data.pickle'), 'w') as f:
					pickle.dump({'convolved_design_matrix':convolved_design_matrix,'sample_times':sample_times}, f)
				#np.save(os.path.join(self.stageFolder('processed/mri/PRF/'), 'convolved_design_matrix_for_average_data.npy'),convolved_design_matrix)

			if animate_dm:
				pl.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
				ims = []
				f=pl.figure(figsize=(4,4))
				if convolve:
					dm = convolved_design_matrix
					timepoints = np.linspace(0,len(dm),243)
					#timepoints = np.linspace(0,len(dm),ssr)
					#lelijke oplossing met hardcoded 243, maar het werkt.
					#timepoints hier zijn er maar 30 van voor de convolved versie. Waarom?
					# voor convolved is die 7290, voor raw is hij 243. 7290 / 243 = 30

				else:
					dm = raw_design_matrix
					timepoints = np.arange(len(dm))
				for t in timepoints[:-1]:
				#original: 
				#for t in timepoints:
				#calling on out of bounds list position otherwise for the convolved one
					s=f.add_subplot(111)
					im=pl.imshow(dm[t,:,:],origin='lowerleft',interpolation='nearest',cmap='gray')
					pl.clim(0,np.max(dm))
					pl.axis('off')
					ims.append([im])
				ani = animation.ArtistAnimation(f, ims)#, interval=5, blit = True, repeat_delay = 1000)
				if not convolve:
					mywriter = animation.FFMpegWriter(fps = 1)
				else:
					mywriter = animation.FFMpegWriter(fps = 1)
					#mywriter = animation.FFMpegWriter(fps = ssr)

				self.logger.info('saving animation')
				plotdir = os.path.join(self.stageFolder('processed/mri/figs/dm_animations'))
				if not os.path.isdir(plotdir): os.mkdir(plotdir)
				if convolve:
					ani.save(os.path.join(plotdir,'convolved_design_matrix_%dx%d.mp4'%(n_pixel_elements, n_pixel_elements)),writer=mywriter,dpi=300,bitrate=200)
				else:
					ani.save(os.path.join(plotdir,'raw_design_matrix_%dx%d.mp4'%(n_pixel_elements, n_pixel_elements)),writer=mywriter,dpi=300,bitrate=200)


	def design_matrices_for_concatenated_data(self, gamma_hrfType = 'doubleGamma', hrf_params = [1,0,0],n_pixel_elements_raw = 101,
		n_pixel_elements_convolved=31,task_conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'],change_type='all_data',run_num=0):
		"""
		This function creates a design matrix per stimulus type
		"""

		#the code below will save it all in 1 PRF folder

		dm_dir = os.path.join(self.stageFolder('processed/mri/PRF/design_matrices/'))
		try:
			os.makedirs(dm_dir)
		except:
			pass

		# create design matrix at higher time resolution to avoid aliasing artifacts in the creation process. 
		for this_condition in task_conditions:

			#this block saves the design matrix in the relevant condition folder.
			# dm_dir = os.path.join(self.stageFolder('processed/mri/%s/design_matrices/')%this_condition)
			# try:
			# 	os.makedirs(dm_dir)
			# except:
			# 	pass


			self.logger.info('Creating design matrix for condition %s'%this_condition)
			self.stimulus_timings(this_condition)

			raw_dms = []
			convolved_dms = []
			for ri,r in enumerate([self.runList[i] for i in self.conditionDict[this_condition]]): 

				this_nii_file = NiftiImage(self.runFile(stage = 'processed/mri', run = r))
				TR = this_nii_file.rtime
				n_TRs = this_nii_file.timepoints
				n_slices = np.min(this_nii_file.extent)

				if TR > 10:
					TR /= 1000
				elif TR < 0.01:
					TR *= 1000		
				

				# select trials
				#if this_condition == 'All':
				which_trials = np.ones(len(r.trial_times)).astype(bool)
				#else: 
				#	which_trials = (r.trial_times[:,0]==np.where(np.array(self.task_names['PRF'])==this_condition)[0][0]) 
				r.orientations = np.radians(r.orientations[which_trials])
				r.trial_times = r.trial_times[which_trials]
				# first, let's create a raw design matrix at TR time resolution with 101 pixel elements
				sample_duration = TR
				bar_width = 4.0/9.0
				mr = PRFModelRun(r, n_samples = n_TRs, n_pixel_elements = n_pixel_elements_raw, sample_duration = sample_duration, bar_width = bar_width)
				self.logger.info('simulating model experiment run %d with %d pixel elements and %1.4f s sample_duration'%(ri,n_pixel_elements_raw, sample_duration))

				#mr.simulate_run(save_images_to_file = dm_dir)
				if this_condition == "PRF_04":
					mr.simulate_run(save_images_to_file = dm_dir)
				else:
					mr.simulate_run()
				
				raw_dms.append(mr.run_matrix)

				#if this_condition == 'All':
				# then, let's create a design matrix at slice time resolution with 31 pixel elements for convolution
				sample_duration = TR/float(n_slices)
				bar_width = 4.0/9.0
				n_samples = n_TRs*n_slices
				mr = PRFModelRun(r, n_samples = n_samples, n_pixel_elements = n_pixel_elements_convolved, sample_duration = sample_duration, bar_width = bar_width)
				self.logger.info('simulating model experiment run %d with %d pixel elements and %1.4f s sample_duration'%(ri,n_pixel_elements_convolved, sample_duration))
				mr.simulate_run()

				self.logger.info('convolving design matrix with hrf')
				run_design = Design(mr.run_matrix.shape[0], mr.sample_duration, subSamplingRatio = 1)
				run_design.rawDesignMatrix = mr.run_matrix.reshape((n_TRs*n_slices,-1)).T
				run_design.convolveWithHRF(hrfParameters = hrf_params)
				convolved_dms.append(run_design.designMatrix.T.reshape(n_TRs*n_slices,n_pixel_elements_convolved,n_pixel_elements_convolved))

			#if this_condition == 'All':
			combined_convolved_dm = np.vstack(convolved_dms)
			del convolved_dms
			self.logger.info('saving convolved design matrix')
			np.save(os.path.join(dm_dir, 'convolved_design_matrix_for_concatenated_data_%s.npy'%(this_condition)),combined_convolved_dm)

			combined_dm = np.vstack(raw_dms)
			del raw_dms
			self.logger.info('saving raw design matrix')
			if change_type != 'all_data':
				np.save(os.path.join(dm_dir, 'raw_design_matrix_for_concatenated_data_%s_%s_%s.npy'%(this_condition,change_type,run_num)),combined_dm)
			else:
				np.save(os.path.join(dm_dir, 'raw_design_matrix_for_concatenated_data_%s.npy'%(this_condition)),combined_dm)



	def weibull_prf_trialtimes(self, task_conditions=['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']):
		"""
		created to setup the correct trial_times.txt files for the weibullPRF experiment.
		will read the 2 stim .txt files, concatenate, order, and add direction of the bar
		based on the trial type in the file
		"""	
		for this_condition in task_conditions:

			self.logger.info('creating trial_times.txt for condition %s'%this_condition)

			for ri,r in enumerate([self.runList[i] for i in self.conditionDict[this_condition]]):
				#setup empty arrays to concatenate/vstack with
				trial_times_directions = np.empty((0,4))
				trial_times = np.empty((0,3))

				#get the list of files in the folder
				subfile_list = os.listdir(self.runFolder(stage = 'processed/behavior', run = r))

				if "trial_times.txt" in subfile_list:
					#remove the trial_times.txt from the list if it alreay exists.
					subfile_list.remove("trial_times.txt")
				if "trial_times.txt~" in subfile_list:
					subfile_list.remove("trial_times.txt~")
					
					#read in each subfile, attach to trial_times
				for subfile in subfile_list:
					subfile = os.path.join(self.runFolder(stage = 'processed/behavior', run = r), subfile)
					trial = np.genfromtxt(subfile, dtype = 'float', delimiter =',')
					trial_times = np.concatenate([trial_times, trial])

				# order the trial times based on their starttime (1st row, 0th index)
				# dont really understand what this does, but hey, it works
				trial_times = sorted(trial_times, key=lambda trial_times_entry: trial_times_entry[0])
				#lets add the directions as the 4th row (index 3)
				for row in trial_times:

					if row[2] == 1.0:
						direction_bar = 315.000
					elif row[2] == 2.0:
						direction_bar = 270.000
					elif row[2] == 3.0:
						direction_bar = 225.000
					elif row[2] == 4.0:
						direction_bar = 180.000
					elif row[2] == 5.0:
						direction_bar = 135.000
					elif row[2] == 6.0:
						direction_bar = 90.000
					elif row[2] == 7.0:
						direction_bar = 45.000				
					elif row[2] == 8.0:
						direction_bar = 0.000
					elif row[2] == 999.0:
						direction_bar = 999.000

					row = np.append(row, direction_bar)
					trial_times_directions = np.vstack([trial_times_directions, row])

				#now save in the correct folder as trial_times.txt with tabs
				out_name = os.path.join(self.runFolder(stage = 'processed/behavior', run = r), 'trial_times.txt')
				np.savetxt(out_name, X = trial_times_directions, delimiter = '\t', fmt="%s")


	def stimulus_timings(self, this_condition):

		"""
		this reads out the stimulus timing and directions from trial_times.txt for each condition.
		this version has been modified for the WeibullPRF experiment setup.

		"""

		for ri,r in enumerate([self.runList[i] for i in self.conditionDict[this_condition]]):

			# load trial_times.txt and get duration and orientation from it.
			this_trial_time_file = np.loadtxt(os.path.join(self.runFolder(stage = 'processed/behavior', run = r),'trial_times.txt'), delimiter="\t")
			r.trial_times = this_trial_time_file
			r.trial_duration = r.trial_times[:,1]
			r.orientations = r.trial_times[:,3]




	def setup_fit_PRF_on_concatenated_data(self, mask_file_name = 'bet_mask', postFix = [], n_jobs = 1, fit_on_all_data=True,n_slices = 30,
				n_vox_per_ROI=100,plotbool=False, model='OG',hrf_type='canonical',n_pixel_elements_raw = 101,n_pixel_elements_convolved=31,r_squared_threshold=0.1,
				use_shared_memory = False,slice_no=0,change_type='all_data',run_num=0,condition_index=0,all_conditions=np.array(['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10'])):
		"""
		docstring for fit_PRF
		"""

		# set tempdir 
		if socket.gethostname() == 'aeneas':
			#tempdir = '/home/vanes/temp'
			tempdir = '/home/voppel/temp'
		else:
			tempdir = '/dev/shm'
			#
		self.logger.info('starting PRF fit procedure')
		total_elapsed_time = 0

		# load anatomical mask - not large enough to store in shared data
		mask_file = NiftiImage(os.path.join(self.stageFolder('processed/mri/masks/anat'), mask_file_name +  '.nii.gz'))
		cortex_mask = np.array(mask_file.data, dtype = bool)

		# in case of not all fit, let's add stats mask
		if not fit_on_all_data:
			# and the results frames
			filename = os.path.join(self.stageFolder('processed/mri/PRF/'), 'frames.pickle')
			with open(filename) as f:
				picklefile = pickle.load(f)
			stats_frames = picklefile['stats_frames']
			results_frames = picklefile['results_frames']
			stats = NiftiImage(os.path.join(self.stageFolder('processed/mri/PRF/'), 'corrs_' + mask_file_name + '_' + '_'.join(postFix)  + '_' + 
				model + '_All_hrf_' + hrf_type+ '.nii.gz')).data[stats_frames['r_squared'],:]
			cortex_mask_size = cortex_mask.sum()
			cortex_mask *= (stats>r_squared_threshold)
			stat_mask_size = cortex_mask.sum()
			self.logger.info('%d/%d voxels selected(%.1f%s), R2 threshold = %.6f'%(stat_mask_size,cortex_mask_size,stat_mask_size/cortex_mask_size*100,'%',r_squared_threshold))
			del stats

		# load and concatenate the data
		all_data = []
		n_TRs = []
		import time as t
		t0=t.time()
		
		for condition in all_conditions[condition_index]:#['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']:
			# for ri,r in enumerate([self.runList[i] for i in self.conditionDict[condition]]): 
			r = [self.runList[i] for i in self.conditionDict[condition]][0]
			this_nii_file = self.runFile(stage = 'processed/mri', run = r, postFix=postFix)
			n_TRs.append(NiftiImage(this_nii_file).timepoints)
			self.logger.info('loading %s'%this_nii_file)
			all_data.append(NiftiImage(this_nii_file).data[:,cortex_mask])
		n_slices = np.min(NiftiImage(this_nii_file).extent)
		all_data = np.vstack(all_data)
		self.logger.info('Loading data lasted %.3f seconds'%(t.time()-t0))

		# the v1 mask has 847 voxels as true
		#TRs of each file = 300, 10x 300 = 3000
		#all_data.shape = 3000, 847

		# for ri,r in enumerate([self.runList[i] for i in self.conditionDict['PRF_01','PRF_02','PRF_03','PRF_04','PRF_05','PRF_06','PRF_07','PRF_08','PRF_09','PRF_10']]): 
		# 	this_nii_file = self.runFile(stage = 'processed/mri', run = r, postFix=postFix)
		# 	n_TRs.append(NiftiImage(this_nii_file).timepoints)
		# 	self.logger.info('loading %s'%this_nii_file)
		# 	all_data.append(NiftiImage(this_nii_file).data[:,cortex_mask])
		# n_slices = np.min(NiftiImage(this_nii_file).extent)
		# all_data = np.vstack(all_data)
		# self.logger.info('Loading data lasted %.3f seconds'%(t.time()-t0))

		# convert it into a shared variable
		# filename = os.path.join(tempdir, 'all_data.dat')
		# fp = np.memmap(filename, dtype='float32', mode='write', shape=all_data.shape)
		# fp[:] = all_data[:]
		# all_data_shared = np.memmap(filename, dtype='float32', mode='readonly', shape=all_data.shape)
		# del all_data, fp


		if hrf_type != 'canonical':
			hrf_nifti_filename = os.path.join(self.stageFolder('pRFocessed/mri/PRF/'), 'mean_hrf_parameters.nii.gz') 
			all_hrf_parameters = NiftiImage(hrf_nifti_filename).data[:,cortex_mask]
			# if use_shared_memory:		
			# convert it into a shared variable
			# filename = os.path.join(tempdir, 'all_hrf_parameters.dat')
			# fp = np.memmap(filename, dtype='float32', mode='write', shape=all_hrf_parameters.shape)
			# fp[:] = all_hrf_parameters[:]
			# all_hrf_parameters_shared = np.memmap(filename, dtype='float32', mode='readonly', shape=all_hrf_parameters.shape)
			# del all_hrf_parameters, fp
		else:
			all_hrf_parameters = np.tile([1,0,0],(cortex_mask.sum(),1)).T

		# number slices
		slices = (np.ones(cortex_mask.shape).T * np.arange(cortex_mask.shape[0])).T[cortex_mask]
		slices_in_full = (np.ones(cortex_mask.shape).T * np.arange(cortex_mask.shape[0])).T

		# figure out roi label per voxel and how many there are for when plotting individual voxels in Dumoulin_fit
		anatRoiFileNames = subprocess.Popen('ls ' + self.stageFolder( stage = 'processed/mri/masks/anat/' ) + '*' + standardMRIExtension, shell=True, stdout=PIPE).communicate()[0].split('\n')[0:-1]
		anatRoiFileNames = [anRF for anRF in anatRoiFileNames if np.all([np.any(['bh' in anRF,'lh' in anRF,'rh' in anRF]),'cortex' not in anRF])]
		roi_names = np.zeros_like(slices_in_full).astype('string')
		for this_roi in anatRoiFileNames:
			roi_nifti = NiftiImage(this_roi).data.astype('bool')
			roi_names[roi_nifti] = (this_roi.split('/')[-1]).split('.')[1]
		roi_names[roi_names=='0.0'] = 'unkown_roi'
		roi_names = roi_names[cortex_mask]
		roi_count = {}
		for roi in np.unique(roi_names):
			roi_count[roi] = np.size(roi_names[roi_names==roi]) 

		if fit_on_all_data:
			#task_conditions = ['All']
			# task_conditions = ['Aĺl']#
			if len(condition_index) ==1:
				task_conditions = all_conditions[condition_index]
			else:
				task_conditions = ['All_%d_%d'%(condition_index[0],condition_index[-1])]
			all_params = None
			results_frames = None
		else:
			# load in the results from the All condition
			all_params = NiftiImage(os.path.join(self.stageFolder('processed/mri/PRF/'), 'results_' + mask_file_name + '_' + '_'.join(postFix)  + '_' + model + 
				'_All_%d_%d_hrf_'%(condition_index[0],condition_index[-1])+ hrf_type+ '.nii.gz')).data[:,cortex_mask]
			# # put all_params in shared mem
			# filename = os.path.join(tempdir, 'all_params.dat')
			# fp = np.memmap(filename, dtype='float32', mode='write', shape=all_params.shape)
			# fp[:] = all_params[:]	 
			# all_params_shared = np.memmap(filename, dtype='float32', mode='readonly', shape=all_params.shape)
			# del all_params, fp

		# set up empty arrays for saving the data
		all_results = {}
		all_corrs = {}
		for condition in task_conditions:
			all_results[condition] = np.zeros([26] + list(cortex_mask.shape))
			all_corrs[condition] = np.zeros([5] + list(cortex_mask.shape))

		# create plot dir for individual voxel plots
		if plotbool: 
			plotbase = os.path.join(self.stageFolder('processed/mri/figs/PRF_fit_plots'))
			if not os.path.isdir(plotbase): os.mkdir(plotbase)
			if change_type == 'all_data':
				this_plot_dir = os.path.join(plotbase,'%s_%s_%s_%s_%s'%(mask_file_name,model,'_'.join(postFix),'_'.join(task_conditions),time_module.strftime("%d.%m.%Y")))
			else:
				this_plot_dir = os.path.join(plotbase,'%s_%s_%s_%s_%s_%d_%s'%(mask_file_name,model,'_'.join(postFix),'_'.join(task_conditions),change_type,run_num,time_module.strftime("%d.%m.%Y")))
			try:
				os.mkdir(this_plot_dir)
			except:
				if os.path.isdir(this_plot_dir):
					self.logger.info('not making plot dir as it exists already')
				else:
					self.logger.info('unkown error in plotdir generation')
		else:
			plotdir = []

		# load raw dms
		raw_dms = {}
		for this_condition in task_conditions:
			self.logger.info('loading raw design matrix for condition %s'%this_condition)
			if change_type == 'all_data':
				all_dms = []
				for sub_condition in all_conditions[condition_index]:
					all_dms.extend(np.load(os.path.join(self.stageFolder('processed/mri/PRF/design_matrices/'),'raw_design_matrix_for_concatenated_data_%s.npy'%(sub_condition))))
				raw_dms[this_condition] = np.array(all_dms)
			else:
				raw_dms[this_condition] = np.load(os.path.join(self.stageFolder('processed/mri/PRF/design_matrices/'),'raw_design_matrix_for_concatenated_data_%s_%s_%s.npy'%(this_condition,change_type,run_num)))



		if fit_on_all_data:
			convolved_dms_to_concatenate = []

			self.logger.info('loading convolved design matrix')
			for this_condition in task_conditions:
				for sub_condition in all_conditions[condition_index]:
					convolved_dms_to_concatenate.extend(np.load(os.path.join(self.stageFolder('processed/mri/PRF/design_matrices/'),'convolved_design_matrix_for_concatenated_data_%s.npy'%(sub_condition))))
					self.logger.info('loading convolved design matrix for condition %s'%sub_condition)
				convolved_dms_to_concatenate = np.array(convolved_dms_to_concatenate)

				# convolved_dm_this_condition= np.load(os.path.join(self.stageFolder('processed/mri/PRF/design_matrices/'),'convolved_design_matrix_for_concatenated_data_%s.npy'%(this_condition)))
				# convolved_dms_to_concatenate.append(convolved_dm_this_condition)
			#convolved_dm = np.vstack(convolved_dms_to_concatenate)

			convolved_dm = convolved_dms_to_concatenate
			del convolved_dms_to_concatenate
			convolved_dm = convolved_dm.reshape(-1,n_pixel_elements_convolved**2)
			valid_regressors = convolved_dm.sum(axis = 0) != 0
			#705 valid regressors van 961
			convolved_dm = convolved_dm[:,valid_regressors]
			# dus 705 valid points to work with.

			# self.logger.info('loading convolved design matrix')
			# convolved_dm = np.load(os.path.join(self.stageFolder('processed/mri/PRF/design_matrices/'),'convolved_design_matrix_for_concatenated_data_All.npy')).reshape(-1,n_pixel_elements_convolved**2)
			# valid_regressors = convolved_dm.sum(axis = 0) != 0
			# convolved_dm = convolved_dm[:,valid_regressors]

			# put in shared mem
			# filename = os.path.join(tempdir, 'convolved_dm.dat')
			# fp = np.memmap(filename, dtype='float32', mode='write', shape=convolved_dm.shape)
			# fp[:] = convolved_dm[:]
			# convolved_dm_shared = np.memmap(filename, dtype='float32', mode='readonly', shape=convolved_dm.shape)
			# del convolved_dm, fp
		else:
			valid_regressors = []
			convolved_dm = []

		fit_start_time = time_module.time()
		if slice_no != None:
			if (slices == slice_no).sum() > 0:
				slice_iterable = [slice_no]
			else: 
				slice_iterable = None
		else:
			slice_iterable = np.unique(slices)

		# estimate fit duration
		if fit_on_all_data:
			minute_per_voxel = 0.07
		else:
			minute_per_voxel = 0.35

		if slice_iterable != None:
			num_voxels = np.sum(cortex_mask[np.array(slice_iterable).astype(int)])
		else:
			num_voxels = cortex_mask.sum()
		estimated_fit_duration = num_voxels * minute_per_voxel
		self.logger.info('starting PRF model fits on %d voxels total'%(int(num_voxels)))
		self.logger.info('estimated total duration: %dm (%.1fh)' % (estimated_fit_duration,estimated_fit_duration/60.0))
		
		# now loop over slices and fit PRF model for voxels parallel
		if slice_iterable != None:
			for sl in slice_iterable:

				plotdir = os.path.join(this_plot_dir,'slice_%d'%sl)
				if os.path.isdir(plotdir): shutil.rmtree(plotdir)
				os.mkdir(plotdir)

				voxels_in_this_slice = (slices == sl)
				voxels_in_this_slice_in_full = (slices_in_full == sl)
				these_roi_names = roi_names[voxels_in_this_slice]

				if fit_on_all_data:
					this_slice_convolved_dm = convolved_dm[int(sl)::int(n_slices)]

				else:
					this_slice_convolved_dm = []
				randints_for_plot = [(np.random.randint(roi_count[these_roi_names[voxno]])<n_vox_per_ROI) for voxno in range(voxels_in_this_slice.sum())]
				self.logger.info('now fitting pRF models on slice %d, with %d voxels' % (sl, voxels_in_this_slice.sum()))
				res = Parallel(n_jobs = np.min([n_jobs,voxels_in_this_slice.sum()]), verbose = 9)(
							delayed(fit_PRF_on_concatenated_data)(
							data_shared 					= all_data,
							voxels_in_this_slice 			= voxels_in_this_slice,
							n_TRs 							= n_TRs,
							n_slices 						= n_slices,
							fit_on_all_data 				= fit_on_all_data,
							plotbool 						= plotbool,
							raw_design_matrices 			= raw_dms, 
							dm_for_BR 						= this_slice_convolved_dm, 
							valid_regressors 				= valid_regressors,
							n_pixel_elements_convolved		= self.n_pixel_elements_convolved,
							n_pixel_elements_raw 			= self.n_pixel_elements_raw,
							plotdir 						= plotdir,
							voxno 							= voxno,
							slice_no 						= sl,
							randint 						= randints_for_plot[voxno],
							roi 							= these_roi_names[voxno],
							TR 								= self.TR,
							model 							= model,
							hrf_params_shared				= all_hrf_parameters,
							all_results_shared				= all_params,
							conditions 						= task_conditions,
							results_frames 					= results_frames,
							postFix 						= postFix,
							)
						for voxno in range(voxels_in_this_slice.sum()))	

				for condition in task_conditions:
					# insert this slice's results in the whole brain variables
					all_corrs[condition][:, cortex_mask * voxels_in_this_slice_in_full] = np.array([rs[1].values() for rs in res]).T
					all_results[condition][:, cortex_mask * voxels_in_this_slice_in_full] = np.array([rs[0][condition].values() for rs in res]).T

			results_frames = {}
			for ki,key in enumerate(res[0][0][task_conditions[0]].keys()):
				results_frames[key] = ki
			stats_frames = {}
			for ki,key in enumerate(res[0][1].keys()):
				stats_frames[key] = ki				
			# save frames
			filename = os.path.join(self.stageFolder('processed/mri/PRF/'), 'frames.pickle')
			with open(filename, 'w') as f:
				pickle.dump({'results_frames':results_frames,'stats_frames':stats_frames}, f)

			fit_end_time = time_module.time()
			fit_time_this_condition = (fit_end_time-fit_start_time)/60.
			self.logger.info('fitting this condition lasted: %dm'%(fit_time_this_condition))
			total_elapsed_time += fit_time_this_condition
			self.logger.info('total fit time: %dm'%(total_elapsed_time))
			
			# and save for every condition
			# not that the corrs file will be the same for the fix/color/speed condition, but this is necessary for file correspondence with average fit
			for condition in task_conditions:

				self.logger.info('saving coefficients and correlations of PRF fits for condition %s'%condition)
				if change_type == 'all_data':
					filename = '%s_%s_%s_%s_hrf_%s'%(mask_file_name,'_'.join(postFix),model,condition,hrf_type)
				else:
					filename = '%s_%s_%s_%s_hrf_%s_%s_%d'%(mask_file_name,'_'.join(postFix),model,condition,hrf_type,change_type,run_num)

				if slice_no == None:
					stage_dir = self.stageFolder('processed/mri/PRF/')
					save_filename = filename
				else:
					stage_dir = os.path.join(self.stageFolder('processed/mri/PRF/'),filename)
					save_filename = filename + '_sl_%d'%slice_no

				if not os.path.isdir(stage_dir): os.mkdir(stage_dir)


				# replace infs in correlations with the maximal value of the rest of the array.
				all_corrs[condition][np.isinf(all_corrs[condition])] = all_corrs[condition][-np.isinf(all_corrs[condition])].max() + 1.0
				corr_nii_file = NiftiImage(all_corrs[condition])
				corr_nii_file.header = mask_file.header
				corr_nii_file.save(os.path.join(stage_dir,'corrs_'+save_filename+'.nii.gz'))

				results_nii_file = NiftiImage(all_results[condition])
				results_nii_file.header = mask_file.header
				results_nii_file.save(os.path.join(stage_dir,'results_'+save_filename+'.nii.gz'))
		else:
			self.logger.info('no voxels in slice %d'%slice_no)
