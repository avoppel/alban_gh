# !/usr/bin/env python
# encoding: utf-8

# import python packages:
import datetime, os, sys
from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as pl

from IPython import embed as shell

class PRFModelTrial(object):
	"""docstring for PRFModelTrial"""
	def __init__(self, orientation, n_elements, n_samples, sample_duration, bar_width = 0.2):
		super(PRFModelTrial, self).__init__()
		
		# the orientation needs to be flipped in order to match the experiment. 
		# First we need to invert it, as the angle moved clockwise in the experiment, while it moves counter-clockwise in the rotation matrix
		# then, we need to subtract pi/2, because 0 is shifted 90 degrees in relation to the experiment
		#self.orientation = orientation + np.pi/2.0 
		self.orientation = orientation - np.pi# - np.pi/2.0

		self.n_elements = n_elements
		self.n_samples = n_samples
		self.sample_duration = sample_duration
		self.bar_width = bar_width
		
		self.rotation_matrix = np.matrix([[cos(self.orientation), sin(self.orientation)],[-sin(self.orientation), cos(self.orientation)]])

		x, y = np.meshgrid(np.linspace(-1,1,self.n_elements), np.linspace(-1,1,self.n_elements))
		self.xy = np.matrix([x.ravel(), y.ravel()]).T  
		self.rotated_xy = np.array(self.xy * self.rotation_matrix)
		self.ecc_test = (np.array(self.xy) ** 2).sum(axis = 1) <= 1.0
	
	def in_bar(self, time = 0):
		"""in_bar, a method, not Ralph."""
		# a bar of self.bar_width width
		position = 2.0 * ((time * (1.0 + self.bar_width / 2.0)) - (0.5 + self.bar_width / 6.0))

		#original
		#position = 2.0 * ((time * (1.0 + self.bar_width / 2.0)) - (0.5 + self.bar_width / 4.0))


		# position = 2.0 * ((time * (1.0 + self.bar_width)) - (0.5 + self.bar_width / 2.0))
		extent = [-self.bar_width/2.0 + position, self.bar_width/2.0 + position] 
		# rotating the xy matrix itself allows us to test only the x component 
		return ((self.rotated_xy[:,0] >= extent[0]) * (self.rotated_xy[:,0] <= extent[1]) * self.ecc_test).reshape((self.n_elements, self.n_elements))
	
	def pass_through(self):
		"""pass_through models a single pass-through of the bar, 
		with padding as in the padding list for start and end."""
		
		self.pass_matrix = np.array([self.in_bar(i) for i in np.linspace(0.0, 1.0, self.n_samples, endpoint = True)])

class PRFModelRun(object):
	"""docstring for PRFModelRun"""
	def __init__(self, run,n_samples, n_pixel_elements, sample_duration = 0.6, bar_width = 0.2):
		super(PRFModelRun, self).__init__()
		self.run = run
		self.n_samples = n_samples
		self.sample_duration = sample_duration

		self.n_pixel_elements = n_pixel_elements
		self.bar_width = bar_width
		
		self.orientation_list = self.run.orientations

	def simulate_run(self, save_images_to_file = None):
		"""docstring for simulate_run"""
		self.sample_times = np.arange(0, self.n_samples * self.sample_duration, self.sample_duration)[:self.n_samples] 
		
		self.run_matrix = np.zeros((self.sample_times.shape[0], self.n_pixel_elements, self.n_pixel_elements))
		
		for i in range(len(self.run.trial_times)):
		#0-11 = 12 keer (12 trials)

			# pak de goede lines eruit.
			#eerst de begintijd; daarna de eindtijd
			samples_in_trial = (self.sample_times >= (self.run.trial_times[i][0])) * (self.sample_times < ((self.run.trial_times[i][0]) + self.run.trial_times[i][1]))
			
			#origineel:
			#samples_in_trial = (self.sample_times >= (self.run.trial_times[i][1])) * (self.sample_times < (self.run.trial_times[i][2]))

 			#if (self.run.trial_times[i][0] != 'Fix_no_stim') * (self.run.trial_times[i][0] != 3):
 			#ik denk dat dit nog eleganter kan. eerste zal altijd true zijn denk ik.
			if (self.run.trial_times[i][0] != 'Fix_no_stim') * (self.run.trial_times[i][2] != 999.0):
				pt = PRFModelTrial(orientation = self.orientation_list[i], n_elements = self.n_pixel_elements, n_samples = samples_in_trial.sum(), sample_duration = self.sample_duration, bar_width = self.bar_width)
				pt.pass_through()
				self.run_matrix[samples_in_trial] = pt.pass_matrix
		if save_images_to_file != None:
			for i in range(self.run_matrix.shape[0]):
				if i < 50:
					f = pl.figure()
					s = f.add_subplot(111)
					pl.imshow(self.run_matrix[i])
					pl.savefig(save_images_to_file + '_' + str(i) + '.pdf')

	def simulate_run_original(self, save_images_to_file = None):
		"""docstring for simulate_run"""
		self.sample_times = np.arange(0, self.n_samples * self.sample_duration, self.sample_duration)[:self.n_samples] 
		
		self.run_matrix = np.zeros((self.sample_times.shape[0], self.n_pixel_elements, self.n_pixel_elements))
		
		for i in range(len(self.run.trial_times)):
		#0-11 = 12 keer (12 trials)

			samples_in_trial = (self.sample_times >= (self.run.trial_times[i][1])) * (self.sample_times < (self.run.trial_times[i][2]))

			if (self.run.trial_times[i][0] != 'Fix_no_stim') * (self.run.trial_times[i][2] != 999.0):
				pt = PRFModelTrial(orientation = self.orientation_list[i], n_elements = self.n_pixel_elements, n_samples = samples_in_trial.sum(), sample_duration = self.sample_duration, bar_width = self.bar_width)
				pt.pass_through()
				self.run_matrix[samples_in_trial] = pt.pass_matrix
				
		if save_images_to_file != None:
			for i in range(self.run_matrix.shape[0]):
				if i < 50:
					f = pl.figure()
					s = f.add_subplot(111)
					pl.imshow(self.run_matrix[i])
					pl.savefig(save_images_to_file + '_' + str(i) + '.pdf')



					