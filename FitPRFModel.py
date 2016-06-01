# !/usr/bin/env python
# encoding: utf-8
"""
Session.py

Created by Tomas HJ Knapen on 2009-11-26.
Copyright (c) 2009 TK. All rights reserved.
"""

# import python packages:
from sklearn.linear_model import Ridge
from IPython import embed as shell
import numpy as np
import pylab as pl
from scipy.stats import pearsonr,spearmanr,kendalltau
from scipy import ndimage
from scipy import interpolate
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit, report_errors
from skimage.morphology import disk
from matplotlib import cm 
import seaborn as sn
from sympy.solvers import solve
from sympy import Symbol, exp
import hrf_estimation as he
import time as t
sn.set(style="ticks")

# import toolbox funcitonality
from Tools.Operators.ImageOperator import *
# from Tools.other_scripts.plotting_tools import *

class gpf(object):
	def __init__(self, design_matrix, max_eccentricity, n_pixel_elements, ssr, rtime,slice_no):#add_empty_trs=0, tr_per_trial=0, n_orientations=9):
		self.design_matrix = design_matrix
		self.max_eccentricity = max_eccentricity
		self.n_pixel_elements = n_pixel_elements
		self.ssr = ssr
		self.rtime = rtime	
		self.slice = slice_no

		X = np.linspace(-max_eccentricity, max_eccentricity, n_pixel_elements)
		Y = np.linspace(-max_eccentricity, max_eccentricity, n_pixel_elements)
		self.MG = np.meshgrid(X, Y)

	#define model function and pass independent variables x and y as a list
	def twoD_Gaussian(self,  xo, yo, sigma):
		(x,y) = self.MG
		theta=0
		a = (np.cos(theta)**2)/(2*sigma**2) + (np.sin(theta)**2)/(2*sigma**2)
		b = -(np.sin(2*theta))/(4*sigma**2) + (np.sin(2*theta))/(4*sigma**2)
		c = (np.sin(theta)**2)/(2*sigma**2) + (np.cos(theta)**2)/(2*sigma**2)
		gauss = np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
		gauss[disk((self.n_pixel_elements-1)/2)==0] = 0 
		return gauss

	def raw_model_prediction(self,  xo, yo, sigma):
		g = self.twoD_Gaussian(xo, yo, sigma).reshape(self.n_pixel_elements**2)
		return np.dot(self.design_matrix.reshape(-1,self.n_pixel_elements**2), g)

	# def hrf_model_prediction(self, xo, yo, sigma,hrf_a1,hrf_a2,hrf_b1,hrf_b2,hrf_c):
	def hrf_model_prediction(self, xo, yo, sigma,hrf_params):

		rmp = self.raw_model_prediction( xo, yo, sigma)
		original_size = len(rmp)
		rmp = np.repeat(rmp, self.ssr, axis=0)
		xx = np.arange(0,32,self.rtime/float(self.ssr))
		self.hrf_kernel = hrf_params[0] * he.hrf.spmt(xx) +hrf_params[1]* he.hrf.dspmt(xx) +hrf_params[2] * he.hrf.ddspmt(xx)
		# self.hrf_kernel = doubleGamma(np.arange(0,32,self.rtime/float(self.ssr)),hrf_a1,hrf_a2,hrf_b1,hrf_b2,hrf_c)
		if self.hrf_kernel.shape[0] % 2 == 1:
			self.hrf_kernel = np.r_[self.hrf_kernel, 0]
		self.hrf_kernel /= np.abs(self.hrf_kernel).sum()

		convolved_mp = fftconvolve( rmp, self.hrf_kernel, 'full' )[int(self.slice)::int(self.ssr)][:int(original_size)]
		# convolved_mp /= np.max(convolved_mp) # this ensures that the amplitude parameter reflects a multiplication factor from 1% signal change
		return convolved_mp, self.hrf_kernel

def fitRidge_for_Dumoulin(design_matrix, timeseries, n_iter = 50, compute_score = False, verbose = True,valid_regressors=[],n_pixel_elements=[], alpha = 1.0):
	"""fitRidge fits a design matrix to a given timeseries.
	It computes the coefficients and returns these coefficients
	plus the correlation between the model fit and timeseries.
	"""
	br = Ridge(alpha = alpha)
	br.fit(design_matrix, timeseries)
	predicted_signal = br.coef_ * design_matrix
	srp = list(spearmanr(timeseries, predicted_signal.sum(axis = 1)))
	srp = [srp[0], -np.log10(srp[1])]

	PRF = np.zeros(n_pixel_elements**2)
	PRF[valid_regressors] = br.coef_
	PRF = np.reshape(PRF,(n_pixel_elements,n_pixel_elements))
	maximum = ndimage.measurements.maximum_position(PRF)

	start_params = {}
	start_params['xo'], start_params['yo'] = maximum[1]/float(n_pixel_elements)*2-1, maximum[0]/float(n_pixel_elements)*2-1

	return start_params, PRF, predicted_signal.sum(axis = 1)


def fit_PRF_on_averaged_data(time_course,ci_time_course, design_matrix, n_pixel_elements_raw, n_pixel_elements_convolved, model='OG',plotbool=False, tr_times = [],TR = 1.5,  
						plotdir=[], voxno=[], dm_for_BR = [], valid_regressors = [], slice_no=[], randint=True, roi='unkown_roi',hrf_params=[],all_results=0,
						max_eccentricity = 1,ellipsoid=False,stim_duration_TR=24,n_slices=30,results_frames = [],max_ecc = 1.6,stim_radius=7.5,postFix=[]):
	""""""

	
	orientations = ['0','45','90','135','180','225','270','315','X']
	n_orientations = len(orientations)

	#########################################################################################################################################################################################################################
	#### Initiate Parameters
	#########################################################################################################################################################################################################################

	## initiate parameters:
	params = Parameters()

	# two location parameters
	params.add('xo', value= 0.0 )
	params.add('yo', value= 0.0)
	params.add('ecc',value=0.0,min=0,max=max_ecc,expr='sqrt(xo**2+yo**2)')

	# and a baseline
	params.add('baseline',value=0.0)

	# center parameters:
	# params.add('sigma_center',value=0.1,min=0.0000000001) # this means initialization at 0.1/2 * 15 = 0,75 degrees sd, which amounts to 6.6 degrees fwhm, 
	params.add('sigma_center',value=0.1,min=0.0000000001) # this means initialization at 0.1/2 * 15 = 0,75 degrees sd, which amounts to 6.6 degrees fwhm, 
	params.add('amp_center',value=0.05,min=0.0000000001) # this is initialized at 0.001

	# surround parameters
	params.add('sigma_surround',value=0.5,expr='sigma_center+delta_sigma') # surround size should roughly be 5 times that of the center
	params.add('amp_surround',value=-0.001,max=-0.0000000001,expr='-amp_center+delta_amplitude') # initialized at 10% of center amplitude
	params.add('delta_sigma',value=0.4,min=0.0000000001) # this difference parameter ensures that the surround is always larger than the center
	params.add('delta_amplitude',value=0.049,min=0.0000000001) # this difference parameter ensures that the surround is never deeper than the center is high


	# when fitting an OG model, set all surround and delta parameters to 0 and to not vary and set the expression to None, otherwise it will start to vary anyway
	if model == 'OG':	
		params['amp_surround'].value,params['amp_surround'].vary,params['amp_surround'].expr = 0, False, None
		params['delta_amplitude'].vary, params['delta_sigma'].vary,params['sigma_surround'].vary = False, False, False


	#########################################################################################################################################################################################################################
	#### Prepare data
	########################################################################################################################################################################################################## ###############

	# add empty periods between trials in dm in order to let the model prediction die down
	tr_per_trial = len(time_course)/n_orientations
	add_empty_trs = 20
	padded_dm = np.zeros((len(time_course)+add_empty_trs*n_orientations,n_pixel_elements_raw,n_pixel_elements_raw))
	padd_mask = np.zeros(len(padded_dm)).astype(bool)
	for i in range(n_orientations):
		padd_mask[i*tr_per_trial+add_empty_trs*i:(i+1)*tr_per_trial+add_empty_trs*i] = True		
		padded_dm[i*tr_per_trial+add_empty_trs*i:(i+1)*tr_per_trial+add_empty_trs*i,:,:] = design_matrix[i*tr_per_trial:(i+1)*tr_per_trial,:,:]


	#########################################################################################################################################################################################################################
	#### Prepare fit object and function
	#########################################################################################################################################################################################################################

	# initiate model prediction objec
	ssr = n_slices
	g = gpf(design_matrix = padded_dm, max_eccentricity = max_eccentricity, n_pixel_elements = n_pixel_elements_raw, rtime = TR, ssr = ssr,slice_no=slice_no)#, add_empty_trs=add_empty_trs,tr_per_trial=tr_per_trial,n_orientations=n_orientations)

	# initiate fit funcitonality
	def residual(params, time_course,padd_mask):

		center_model_prediction =  g.hrf_model_prediction(params['xo'].value, params['yo'].value, params['sigma_center'].value,hrf_params)[0] * params['amp_center'].value
		surround_model_prediction =  g.hrf_model_prediction(params['xo'].value, params['yo'].value, params['sigma_surround'].value,hrf_params)[0] * params['amp_surround'].value

		combined_model_prediction =  params['baseline'].value + center_model_prediction + surround_model_prediction 
		return time_course - combined_model_prediction[padd_mask]

	#########################################################################################################################################################################################################################
	#### initialize parameters
	#########################################################################################################################################################################################################################

	if np.size(all_results) == 1:

		## initiate search space with Ridge prefit
		Ridge_start_params, BR_PRF, BR_predicted = fitRidge_for_Dumoulin(dm_for_BR, time_course, valid_regressors=valid_regressors, n_pixel_elements=n_pixel_elements_convolved, alpha=1e14)
		params['xo'].value  = Ridge_start_params['xo']
		params['yo'].value = Ridge_start_params['yo']

	else:

		params['xo'].value = all_results[results_frames['xo']]
		params['yo'].value = all_results[results_frames['yo']]
		params['sigma_center'].value = all_results[results_frames['sigma_center']]
		params['sigma_surround'].value = all_results[results_frames['sigma_surround']]
		params['amp_center'].value = all_results[results_frames['amp_center']]
		params['amp_surround'].value = all_results[results_frames['amp_surround']]
		params['delta_sigma'].value = all_results[results_frames['delta_sigma']]
		params['delta_amplitude'].value = all_results[results_frames['delta_amplitude']]
		params['baseline'].value = all_results[results_frames['baseline']]

		surround_PRF =  g.twoD_Gaussian(params['xo'].value, params['yo'].value,params['sigma_surround'].value) * params['amp_surround'].value
		center_PRF =  g.twoD_Gaussian(params['xo'].value, params['yo'].value,params['sigma_center'].value) * params['amp_center'].value
		BR_PRF = center_PRF + surround_PRF

	#########################################################################################################################################################################################################################
	#### evalute fit
	#########################################################################################################################################################################################################################

	# find optimal parameters:
	minimize(residual, params, args=(), kws={'time_course':time_course,'padd_mask':padd_mask},method='powell')

	#########################################################################################################################################################################################################################
	#### Recreate resulting predictions and PRFs with optimized parameters
	#########################################################################################################################################################################################################################

	trimmed_center_mp = (g.hrf_model_prediction(params['xo'].value, params['yo'].value, params['sigma_center'].value,hrf_params)[0] * params['amp_center'].value)[padd_mask]
	trimmed_surround_mp = (g.hrf_model_prediction(params['xo'].value, params['yo'].value, params['sigma_surround'].value,hrf_params)[0] * params['amp_surround'].value)[padd_mask]
	trimmed_mp = params['baseline'].value + trimmed_center_mp + trimmed_surround_mp 


	raw_center_mp = (g.raw_model_prediction(params['xo'].value, params['yo'].value, params['sigma_center'].value)* params['amp_center'].value)[padd_mask]
	raw_surround_mp = (g.raw_model_prediction(params['xo'].value, params['yo'].value, params['sigma_surround'].value)* params['amp_surround'].value)[padd_mask]
	raw_mp = params['baseline'].value + raw_center_mp + raw_surround_mp

	surround_PRF = g.twoD_Gaussian(params['xo'].value, params['yo'].value,params['sigma_surround'].value) * params['amp_surround'].value
	center_PRF = g.twoD_Gaussian(params['xo'].value, params['yo'].value,params['sigma_center'].value) * params['amp_center'].value
	PRF = center_PRF + surround_PRF

	#########################################################################################################################################################################################################################
	#### Get fit diagnostics
	#########################################################################################################################################################################################################################

	## In a DoG model, the center region is determined by the subtraction of the positive and the negative gaussian. 
	## The size of the positive gaussian is therefore not directly linked to the size of the positive region. 
	## Therefore, the FWHM measure is more appropriate. To get it, we first create the PRF at center position,
	## and select from it the line that runs right through the middle:
	reconstruction_radius = 10
	this_ssr = 1000 
	t = np.linspace(-reconstruction_radius,reconstruction_radius,this_ssr*reconstruction_radius)
	PRF_2D =  params['amp_center'].value * np.exp(-t**2/(2*params['sigma_center'].value**2)) + params['amp_surround'].value * np.exp(-t**2/(2*(params['sigma_surround'].value)**2))
	## then, we fit a spline through this line, and get the roots (the fwhm points) of the spline:
	spline=interpolate.UnivariateSpline(range(len(PRF_2D)),PRF_2D-np.max(PRF_2D)/2,s=0)
	## and compute the distance between them
	try:
		fwhm = ((np.diff(spline.roots())/len(t)*reconstruction_radius) * stim_radius)[0]
	except:
		## when this procedure fails, set fwhm to 0:
		fwhm = 0

	## now find the surround size in the same way
	if (model == 'OG') + (params['amp_surround'].value == 0):
		surround_size = 0
	else:
		spline=interpolate.UnivariateSpline(range(len(PRF_2D)),PRF_2D+np.min(PRF_2D),s=0)
		surround_size = ((np.diff(spline.roots())/len(t)*reconstruction_radius) * stim_radius)[0]

	## EVALUATE FIT QUALITY
	# RSS = np.sum((time_course - trimmed_mp)**2)
	stats = {}
	stats['spearman'] = spearmanr(time_course,trimmed_mp)[0]
	stats['pearson'] = pearsonr(time_course,trimmed_mp)[0]
	stats['RSS'] = np.sum((time_course - trimmed_mp)**2)
	stats['r_squared'] = 1 - stats['RSS']/np.sum((time_course - np.mean(time_course)) ** 2) 
	stats['kendalls_tau'] = kendalltau(time_course,trimmed_mp)[0]

	## SETUP RESULTS DICT
	results={}
	for key in params.keys():
		results[key] = params[key].value
	
	results['ecc'] *= stim_radius# np.linalg.norm([params['xo'].value,params['yo'].value]) * stim_radius
	results['fwhm'] = fwhm
	results['surround_size'] = surround_size 
	results['polar'] = np.arctan2(params['yo'].value,params['xo'].value)
	# if the resulting PRF falls outside of the stimulus radius,
	# set the complex value here to 0 so that it falls off the retmaps	
	if results['ecc'] < (0.9*stim_radius):
		multiplier = stats['r_squared']
	else:
		multiplier = 0.001
	results['real_polar'] = np.cos(results['polar'])*np.arctanh(multiplier)
	results['imag_polar'] = np.sin(results['polar'])*np.arctanh(multiplier)
	results['real_eccen'] = np.cos(results['ecc'])*np.arctanh(multiplier)
	results['imag_eccen'] = np.sin(results['ecc'])*np.arctanh(multiplier)
	results['real_fwhm'] = np.cos(results['fwhm'])*np.arctanh(multiplier)
	results['imag_fwhm'] = np.sin(results['fwhm'])*np.arctanh(multiplier)
	results['SI'] = (params['amp_surround'].value * (params['sigma_surround'].value**2) ) / (params['amp_center'].value * (params['sigma_center'].value**2) )



	#########################################################################################################################################################################################################################
	#### Plot results
	#########################################################################################################################################################################################################################

	# print stats['r_squared']
	if plotbool:# + (stats['r_squared']> 0.1):# * randint:#:# * :# :#* (results['ecc'] < 3) :#:# * * randint ) #* :#* )

		plot_dir = os.path.join(plotdir, '%s'%roi)
		if not os.path.isdir(plot_dir): os.mkdir(plot_dir)

		f=pl.figure(figsize=(18,6)); ri = 2
		sn.set(style="ticks")
		minval = np.min(time_course - ci_time_course)
		maxval = np.max(time_course + ci_time_course)
		for di in range(9):

			this_timecourse = time_course[di*len(time_course)/len(orientations):(di+1)*len(time_course)/len(orientations)]
			this_ci = ci_time_course[di*len(time_course)/len(orientations):(di+1)*len(time_course)/len(orientations)]
			s=f.add_subplot(ri,len(orientations),di+1)
			pl.axhline(results['baseline'],linestyle='-',linewidth=1,color='k')
			pl.plot(tr_times,this_timecourse,'k',linewidth=2,label='data')
			pl.fill_between(tr_times,this_timecourse-this_ci,this_timecourse+this_ci,color='k',alpha=0.15)
			pl.plot(tr_times,trimmed_mp[di*len(trimmed_mp)/len(orientations):(di+1)*len(trimmed_mp)/len(orientations)],'m',linewidth=2,label = 'model')
			pl.plot(tr_times,results['baseline']+trimmed_surround_mp[di*len(trimmed_surround_mp)/len(orientations):(di+1)*len(trimmed_surround_mp)/len(orientations)],'b',linestyle='--',linewidth=1,label = 'surround mp')
			pl.plot(tr_times,results['baseline']+trimmed_center_mp[di*len(trimmed_center_mp)/len(orientations):(di+1)*len(trimmed_center_mp)/len(orientations)],'r',linestyle='--',linewidth=1,label = 'center mp')
			pl.ylim(minval,maxval)
			pl.xlim(0,np.max(tr_times))

			if di == 0:
				if 'psc' in postFix:
					pl.ylabel('% signal change')
				else:
					pl.ylabel('unkown unit')
			else:
				pl.yticks([])
				# pl.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 
			if di == (len(orientations)-1):
				leg = s.legend(fancybox = True, loc = 'best')
				leg.get_frame().set_alpha(0.5)
				if leg:
					for t in leg.get_texts():
					    t.set_fontsize(8)    # the legend text fontsize
			sn.despine(offset=10)

			pl.title(orientations[di])
			# if di == 4:
			# 	pl.xlabel('time (s)')
			# 	# pl.text(len(tr_times)/2.0,np.min(zip(*all_averaged_data))*0.8,'stimulus',verticalalignment='center',fontsize=10)
			# 	pl.xticks([0,int(stim_duration_TR*TR)],['0','%d'%int(stim_duration_TR*TR)])
			# else:
			pl.xticks([])

		s = f.add_subplot(ri,6,ri*6-5)
		# s = f.add_subplot(ri,2,3)
		pl.imshow(BR_PRF,origin='lowerleft',interpolation='nearest',cmap=cm.coolwarm)
		pl.axis('off')
		s.set_title('Ridge PRF')

		s = f.add_subplot(ri,6,ri*6-4)
		# s = f.add_subplot(ri,7,7)
		pl.imshow(PRF,origin='lowerleft',interpolation='nearest',cmap=cm.coolwarm)
		pl.axis('off')
		s.set_title('Direct model PRF')
		
		# pl.tight_layout()
		# pl.savefig(os.path.join(plot_dir, 'vox_%d_%d_%d.pdf'%(slice_no,voxno,n_pixel_elements_raw)))
		# pl.close()

		s = f.add_subplot(ri,6,ri*6-1)
		pl.imshow(np.ones((n_pixel_elements_raw,n_pixel_elements_raw)),cmap='gray')
		pl.clim(0,1)
		# s.text(n_pixel_elements_raw/2,n_pixel_elements_raw/2, '\nHRF parameters: \n\na1: %.2f\na2: %.2f\nb1: %.2f\nb2: %.2f\nc: %.2f'
		# 	 %(hrf_params['hrf_a1'],hrf_params['hrf_a2'],hrf_params['hrf_b1'],hrf_params['hrf_b2'],hrf_params['hrf_c']),horizontalalignment='center',verticalalignment='center',fontsize=12,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
		# s.text(n_pixel_elements_raw/2,n_pixel_elements_raw/2, '\nHRF parameters: \n\n1: %.2f\n2: %.2f\n3: %.2f\nb2: %.2f\nc: %.2f'
		# 	 %(hrf_params['hrf_a1'],hrf_params['hrf_a2'],hrf_params['hrf_b1'],hrf_params['hrf_b2'],hrf_params['hrf_c']),horizontalalignment='center',verticalalignment='center',fontsize=12,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
		pl.axis('off')

		s = f.add_subplot(ri,6,ri*6-2)
		pl.imshow(np.ones((n_pixel_elements_raw,n_pixel_elements_raw)),cmap='gray')
		pl.clim(0,1)
		s.text(n_pixel_elements_raw/2,n_pixel_elements_raw/2, '\nFIT PARAMETERS: \n\nsd center: %.2f\nsd surround: %.2f\namp center: %.6f\namp surround: %.6f\nbaseline: %.6f\n\nDERIVED QUANTIFICATIONS: \n\nr squared: %.2f\necc: %.2f\nFWHM: %.2f\nsurround size: %.2f\nsupression index: %.2f'
			 %(params['sigma_center'].value,params['sigma_surround'].value,params['amp_center'].value,params['amp_surround'].value,params['baseline'].value,stats['r_squared'],results['ecc'],results['fwhm'],results['surround_size'],results['SI']),horizontalalignment='center',verticalalignment='center',fontsize=12,bbox={'facecolor':'white', 'alpha':1, 'pad':10})

		pl.axis('off')

		with sn.axes_style("dark"):
			s = f.add_subplot(ri,6,ri*6-3)
			# pl.axhline(0,linestyle='--',linewidth=2,color='w')
			t = np.linspace(-1,1,n_pixel_elements_raw)
			PRF_2D =  params['amp_center'].value * np.exp(-t**2/(2*params['sigma_center'].value**2)) + params['amp_surround'].value * np.exp(-t**2/(2*(params['sigma_surround'].value)**2))
			PRF_2D_surround =  params['amp_surround'].value * np.exp(-t**2/(2*(params['sigma_surround'].value)**2))
			PRF_2D_center =  params['amp_center'].value * np.exp(-t**2/(2*params['sigma_center'].value**2))
			spline=interpolate.UnivariateSpline(range(len(PRF_2D)),PRF_2D-np.max(PRF_2D)/2,s=0)
			spline_surround=interpolate.UnivariateSpline(range(len(-PRF_2D_surround)),-PRF_2D_surround-np.max(-PRF_2D_surround)/2,s=0)
			spline_center=interpolate.UnivariateSpline(range(len(PRF_2D_center)),PRF_2D_center-np.max(PRF_2D_center)/2,s=0)
			pl.plot(PRF_2D,'k',linewidth=2)
			pl.plot(PRF_2D_surround,'--b',linewidth=1)
			pl.plot(PRF_2D_center,'--r',linewidth=1)
			pix_per_degree = n_pixel_elements_raw / stim_radius*2
			pl.xticks(np.array([-10,0,10])*pix_per_degree+n_pixel_elements_raw/2,[-10,0,10])
			# pl.yticks([-0.15,0,0.15])
			pl.fill_between(spline_center.roots(),np.min(PRF_2D),np.max(PRF_2D),color='r',alpha=0.1)
			pl.fill_between(spline_surround.roots(),np.min(PRF_2D),np.max(PRF_2D),color='b',alpha=0.1)
			pl.fill_between(spline.roots(),np.min(PRF_2D),np.max(PRF_2D),color='k',alpha=0.5)
			pl.ylabel('a.u.')
			pl.xlabel('visual degrees')
			# pl.text(n_pixel_elements_raw/2,-0.05,'FWHM',color='w',horizontalalignment='center',verticalalignment='center',fontsize=12,fontweight='bold')
			s.set_title('2D PRF profile')
			pl.ylim(np.min(PRF_2D),np.max(PRF_2D))

		with sn.axes_style("darkgrid"):
			xx = np.arange(0,32,TR/float(ssr))
			hrf_kernel = hrf_params[0] * he.hrf.spmt(xx) +hrf_params[1]* he.hrf.dspmt(xx) +hrf_params[2] * he.hrf.ddspmt(xx)
			# hrf_kernel = doubleGamma(np.arange(0,32,TR/float(ssr)),hrf_params['hrf_a1'],hrf_params['hrf_a2'],hrf_params['hrf_b1'],hrf_params['hrf_b2'],hrf_params['hrf_c'])
			hrf_kernel /= np.abs(hrf_kernel).sum()
			s = f.add_subplot(ri,6,ri*6)
			pl.plot(hrf_kernel)
			s.set_title('HRF-kernel')
			sn.despine(offset=5)
			pl.xticks(np.linspace(0,len(hrf_kernel),16),np.arange(0,32,2))
			pl.xlabel('time (s)')

		pl.savefig(os.path.join(plot_dir, 'vox_%d_%d_%d.pdf'%(slice_no,voxno,n_pixel_elements_raw)))
		pl.close()

	# self.results_frames = {'polar':0,'delta_amplitude':1,'ecc':2,'xo':3,'yo':4,'real_eccen':5,'amp_center':6,'surround_size':7,'imag_polar':8,'amp_surround':9,'sigma_surround':10,'real_fwhm':11,'imag_eccen':12,'imag_fwhm':13,'real_polar':14,'SI':15,'delta_sigma':16,:'sigma_center':17,'fwhm':18,'baseline':19}

	return results, stats

def fit_PRF_on_concatenated_data(data_shared,voxels_in_this_slice,n_TRs,n_slices,fit_on_all_data,plotbool,raw_design_matrices, dm_for_BR,
	valid_regressors, n_pixel_elements_convolved, n_pixel_elements_raw,plotdir,voxno,slice_no,randint,roi,TR,model,hrf_params_shared,all_results_shared,conditions,
	results_frames,	postFix=[],max_eccentricity=1,max_xy = 5,orientations=['0','45','90','135','180','225','270','315','X'],stim_radius = 7.5):
	"""
	stim_radius lijkt niet veel uit te maken.
	"""
	
	# grab data for this fit procedure from shared memory
	time_course = np.array(data_shared[:,voxels_in_this_slice][:,voxno])
	hrf_params = np.array(hrf_params_shared[:,voxels_in_this_slice][:,voxno])

	n_orientations = len(orientations)

	#to plot the time course:
	#%pylab
	#shell()
	# then input:
	#pl.plot(range(0,3000), time_course)



	# already initialize the final PRF dict
	PRFs = {}

	if fit_on_all_data:

		#########################################################################################################################################################################################################################
		#### Instantiate parameters 
		#########################################################################################################################################################################################################################

		## initiate search space with Ridge prefit
		Ridge_start_params, PRFs['Ridge'], BR_predicted = fitRidge_for_Dumoulin(dm_for_BR, time_course, valid_regressors=valid_regressors, n_pixel_elements=n_pixel_elements_convolved, alpha=1e14)
		# params['xo_%s'%conditions[0]].value = Ridge_start_params['xo']
		# params['yo_%s'%conditions[0]].value = Ridge_start_params['yo']

		## initiate parameters:
		params = Parameters()
		
		# one baseline parameter
		params.add('baseline',value=0.0)

		# two location parameters
		# xo_yo_search_width_in_degrees = 2

		# these lines work with PRF_01 etc.
		params.add('xo_%s'%conditions[0], value = Ridge_start_params['xo'])
		params.add('yo_%s'%conditions[0], value = Ridge_start_params['yo'])

		# # fit method with ecc boundary
		# params.add('xo_%s'%conditions[0], value = 0.0,min=-max_ecc,max=max_ecc)# if xo_%s>0 else -(sqrt(max_ecc**2-abs(yo_%s)**2) - abs(delta_xo_%s))'%(tuple(np.repeat(conditions[0],5))))
		# params.add('yo_%s'%conditions[0], value = 0.0,min=0,expr='(sqrt(max_ecc**2-abs(xo_%s)**2) - delta_yo_%s'%(tuple(np.repeat(conditions[0],2))))# if yo_%s>0 else -(sqrt(max_ecc**2-abs(xo_%s)**2) - abs(delta_yo_%s))'%(tuple(np.repeat(conditions[0],5))))

		# # these parameters ensure a maximum ecc
		# params.add('delta_yo_%s'%conditions[0], value=0.0,min=0)#,expr='sqrt(max_ecc**2-abs(xo_%s)**2)*2 - delta_delta_yo_%s'%(tuple(np.repeat(conditions[0],2))))
		# params.add('max_ecc',value=max_ecc,vary=False)
		# params.add('sign_yo_%s'%conditions[0],value=0.01)

		# V3 like eccen-sd relation
		# intercept = 0.7 /stim_radius
		# slope = 0.3 / stim_radius
		# start_size = intercept + np.linalg.norm(Ridge_start_params['xo'],Ridge_start_params['yo']) * slope
		params.add('sigma_center_%s'%conditions[0],value=0.1,min=0.0)#min=0.01 # this means initialization at 0.1 * 7.5 = 0.75 degrees, with minimum of 0.075 degrees
		params.add('amp_center_%s'%conditions[0],value=0.05,min=0.0)#min=0.01 # this is initialized at 0.001


		# surround parameters
		params.add('delta_sigma_%s'%conditions[0],value=0.4,min=0.0) # this difference parameter ensures that the surround is always larger than the center#,min=0.0000000001
		params.add('sigma_surround_%s'%conditions[0],value=0.3,expr='sigma_center_%s+delta_sigma_%s'%(conditions[0],conditions[0])) # surround size should roughly be 5 times that of the center
		params.add('delta_amplitude_%s'%conditions[0],value=0.045,min=0.0) # this difference parameter ensures that the surround is never deeper than the center is high,min=0.0000000001
		params.add('amp_surround_%s'%conditions[0],value=-0.005,max=0.0,expr='-amp_center_%s+delta_amplitude_%s'%(conditions[0],conditions[0])) # initialized at 10% of center amplitude #max=-0.0000000001,

		# when fitting an OG model, set all surround and delta parameters to 0 and to not vary and set the expression to None, otherwise it will start to vary anyway
		if model == 'OG':	
			params['amp_surround_%s'%conditions[0]].value,params['amp_surround_%s'%conditions[0]].vary,params['amp_surround_%s'%conditions[0]].expr = 0, False, None
			params['delta_amplitude_%s'%conditions[0]].vary, params['delta_sigma_%s'%conditions[0]].vary,params['sigma_surround_%s'%conditions[0]].vary = False, False, False

		# params['delta_yo_%s'%conditions[0]].value = sqrt(max_ecc**2-abs(Ridge_start_params['xo'])**2) - Ridge_start_params['yo']
		# params['delta_delta_yo_%s'%conditions[0]].value =  sqrt(max_ecc**2-abs(Ridge_start_params['xo'])**2)*2 + params['delta_yo_%s'%conditions[0]].value

	else:

		#########################################################################################################################################################################################################################
		#### INITIATING PARAMETERS with all results
		#########################################################################################################################################################################################################################

		# grab data for this fit procedure from shared memory
		all_results = np.array(all_results_shared[:,voxels_in_this_slice][:,voxno])

		## initiate parameters:
		params = Parameters()

		# shared baseline param:
		params.add('baseline', value = all_results[results_frames['baseline']])

		# location parameters
		for condition in conditions:
			params.add('xo_%s'%condition, value = all_results[results_frames['xo']])
			params.add('yo_%s'%condition, value = all_results[results_frames['yo']])

			# center parameters:
			params.add('sigma_center_%s'%condition,value=all_results[results_frames['sigma_center']]/stim_radius,min=0.0) # this means initialization at 0.05/2 * 15 = 1.5 degrees, ,min=0.0084
			params.add('amp_center_%s'%condition,value=all_results[results_frames['amp_center']],min=0.0) # this is initialized at 0.001 ,min=0.0000000001

			# surround parameters
			params.add('sigma_surround_%s'%condition,value=all_results[results_frames['sigma_surround']]/stim_radius,expr='sigma_center_%s+delta_sigma_%s'%(condition,condition)) # surround size should roughly be 5 times that of the center
			params.add('amp_surround_%s'%condition,value=all_results[results_frames['amp_surround']],max=0.0,expr='-amp_center_%s+delta_amplitude_%s'%(condition,condition)) # initialized at 10% of center amplitudemax=-0.0000000001
			params.add('delta_sigma_%s'%condition,value=all_results[results_frames['delta_sigma']],min=0.0) # this difference parameter ensures that the surround is always larger than the centermin=0.0000000001
			params.add('delta_amplitude_%s'%condition,value=all_results[results_frames['delta_amplitude']],min=0.0) # this difference parameter ensures that the surround is never deeper than the center is highmin=0.0000000001

			# when fitting an OG model, set all surround and delta parameters to 0 and to not vary and set the expression to None, otherwise it will start to vary anyway
			if model == 'OG':	
				params['amp_surround_%s'%condition].value,params['amp_surround_%s'%condition].vary,params['amp_surround_%s'%condition].expr = 0, False, None
				params['delta_amplitude_%s'%condition].vary, params['delta_sigma_%s'%condition].vary,params['sigma_surround_%s'%condition].vary = False, False, False


		g = gpf(design_matrix = raw_design_matrices[conditions[0]], max_eccentricity = max_eccentricity, n_pixel_elements = n_pixel_elements_raw, rtime = TR, ssr = 1,slice_no=slice_no)
		
		# recreate PRFs
		this_surround_PRF = g.twoD_Gaussian(all_results[results_frames['xo']],all_results[results_frames['yo']],
			all_results[results_frames['sigma_surround']]/stim_radius) * all_results[results_frames['amp_surround']]
		this_center_PRF = g.twoD_Gaussian(all_results[results_frames['xo']], all_results[results_frames['yo']],
			all_results[results_frames['sigma_center']]/stim_radius) * all_results[results_frames['amp_center']]
		PRFs['All_fit'] = this_center_PRF + this_surround_PRF

	#########################################################################################################################################################################################################################
	#### Prepare fit object and function
	#########################################################################################################################################################################################################################

	# initiate model prediction object
	ssr = np.round(1/(TR/float(n_slices)))
	
	gpfs = {}
	for condition in conditions:
		gpfs[condition] = gpf(design_matrix = raw_design_matrices[condition], max_eccentricity = max_eccentricity, n_pixel_elements = n_pixel_elements_raw, rtime = TR, ssr = ssr,slice_no=slice_no)

	def residual(params):
		
		# initiate model prediction at baseline value
		combined_model_prediction = np.ones_like(time_course) * params['baseline'].value

		# now loop over conditions, create prediction and add to total prediction
		for condition in conditions:
			combined_model_prediction +=  gpfs[condition].hrf_model_prediction(params['xo_%s'%condition].value, params['yo_%s'%condition].value, 
				params['sigma_center_%s'%condition].value,hrf_params)[0] * params['amp_center_%s'%condition].value
			combined_model_prediction +=  gpfs[condition].hrf_model_prediction(params['xo_%s'%condition].value, params['yo_%s'%condition].value, 
				params['sigma_surround_%s'%condition].value,hrf_params)[0] * params['amp_surround_%s'%condition].value
	
		return time_course - combined_model_prediction

	#########################################################################################################################################################################################################################
	#### evalute fit
	#########################################################################################################################################################################################################################
	
	# optimize parameters
	minimize(residual, params, args=(), kws={},method='powell')

	#########################################################################################################################################################################################################################
	#### Recreate resulting predictions and PRFs with optimized parameters
	#########################################################################################################################################################################################################################

	# initiate model prediction at baseline value
	combined_model_prediction = np.ones_like(time_course) * params['baseline'].value

	# now loop over conditions, create prediction and add to total prediction
	model_predictions = {}
	for condition in conditions:
		this_center_model_prediction = gpfs[condition].hrf_model_prediction(params['xo_%s'%condition].value, params['yo_%s'%condition].value, 
			params['sigma_center_%s'%condition].value,hrf_params)[0] * params['amp_center_%s'%condition].value
		this_surround_model_prediction = gpfs[condition].hrf_model_prediction(params['xo_%s'%condition].value, params['yo_%s'%condition].value, 
			params['sigma_surround_%s'%condition].value, hrf_params)[0] * params['amp_surround_%s'%condition].value
		model_predictions[condition] = this_center_model_prediction + this_surround_model_prediction
		combined_model_prediction += model_predictions[condition]

		# recreate PRFs
		this_center_PRF = gpfs[condition].twoD_Gaussian(params['xo_%s'%condition].value, params['yo_%s'%condition].value,
			params['sigma_center_%s'%condition].value) * params['amp_center_%s'%condition].value
		this_surround_PRF = gpfs[condition].twoD_Gaussian(params['xo_%s'%condition].value, params['yo_%s'%condition].value,
			params['sigma_surround_%s'%condition].value) * params['amp_surround_%s'%condition].value
		PRFs[condition] = this_center_PRF + this_surround_PRF

	#########################################################################################################################################################################################################################
	#### Get fit diagnostics
	#########################################################################################################################################################################################################################

	reconstruction_radius = 10
	this_ssr = 1000 
	t = np.linspace(-reconstruction_radius,reconstruction_radius,this_ssr*reconstruction_radius)
	
	fwhms = {}
	surround_sizes = {}
	for condition in conditions:
		PRF_2D =  params['amp_center_%s'%condition].value * np.exp(-t**2/(2*params['sigma_center_%s'%condition].value**2)) + params['amp_surround_%s'%condition].value * np.exp(-t**2/(2*(params['sigma_surround_%s'%condition].value)**2))
		## then, we fit a spline through this line, and get the roots (the fwhm points) of the spline:
		spline=interpolate.UnivariateSpline(range(len(PRF_2D)),PRF_2D-np.max(PRF_2D)/2,s=0)
		## and compute the distance between them
		try:
			fwhms[condition] = ((np.diff(spline.roots())/len(t)*reconstruction_radius) * stim_radius)[0]
		except:
			## when this procedure fails, set fwhm to 0:
			fwhms[condition] = 0
		
		## now find the surround size in the same way
		if (model == 'OG') + (params['amp_surround_%s'%condition].value == 0):
			surround_sizes[condition] = 0
		else:
			spline = interpolate.UnivariateSpline(range(len(PRF_2D)),PRF_2D+np.min(PRF_2D),s=0)
			surround_sizes[condition] = ((np.diff(spline.roots())/len(t)*reconstruction_radius) * stim_radius)[0]

	## EVALUATE OVERALL MODEL FIT QUALITY
	stats = {}
	stats['spearman'] = spearmanr(time_course,combined_model_prediction)[0]
	stats['pearson'] = pearsonr(time_course,combined_model_prediction)[0]
	stats['RSS'] = np.sum((time_course - combined_model_prediction)**2)
	stats['r_squared'] = 1 - stats['RSS']/np.sum((time_course - np.mean(time_course)) ** 2) 
	stats['kendalls_tau'] = kendalltau(time_course,combined_model_prediction)[0]

	## CREATE SEPERATE RESULTS DICT PER CONDITION
	results = {}
	for condition in conditions:
		results[condition] = {}
		results[condition]['baseline'] = params['baseline'].value
		# params from fit
		for key in params.keys():
			if condition in key:
				if condition in key:
					# leave out the condition in the keys (as the results frames are identical across conditions)
					new_key = key[:-len(condition)-1]
				else:
					new_key = key
				results[condition][new_key] = params[key].value

		results[condition]['ecc'] = np.linalg.norm([params['xo_%s'%condition].value,params['yo_%s'%condition].value]) * stim_radius
		results[condition]['sigma_center'] *= stim_radius
		results[condition]['sigma_surround'] *= stim_radius

		# derived params
		results[condition]['polar'] = np.arctan2(params['yo_%s'%condition].value,params['xo_%s'%condition].value)
		results[condition]['fwhm'] = fwhms[condition]
		results[condition]['surround_size'] = surround_sizes[condition]
		results[condition]['SI'] = ((params['amp_surround_%s'%condition].value * (params['sigma_surround_%s'%condition].value**2) ) 
			/ (params['amp_center_%s'%condition].value * (params['sigma_center_%s'%condition].value**2) ))
		
		# if the resulting PRF falls outside of the stimulus radius,
		# set the multiplier here to 0 so that it falls off the retmaps
		if results[condition]['ecc'] < (stim_radius):
			multiplier = stats['r_squared']
		else:
			multiplier = 0.001

		# here for only voxels within stim region:
		results[condition]['real_polar_stim_region'] = np.cos(results[condition]['polar'])*np.arctanh(multiplier)
		results[condition]['imag_polar_stim_region'] = np.sin(results[condition]['polar'])*np.arctanh(multiplier)
		results[condition]['real_eccen_stim_region'] = np.cos(results[condition]['ecc'])*np.arctanh(multiplier)
		results[condition]['imag_eccen_stim_region'] = np.sin(results[condition]['ecc'])*np.arctanh(multiplier)
		results[condition]['real_fwhm_stim_region'] = np.cos(results[condition]['fwhm'])*np.arctanh(multiplier)
		results[condition]['imag_fwhm_stim_region'] = np.sin(results[condition]['fwhm'])*np.arctanh(multiplier)
		
		# and for all voxels:
		results[condition]['real_polar'] = np.cos(results[condition]['polar'])*np.arctanh(stats['r_squared'])
		results[condition]['imag_polar'] = np.sin(results[condition]['polar'])*np.arctanh(stats['r_squared'])
		results[condition]['real_eccen'] = np.cos(results[condition]['ecc'])*np.arctanh(stats['r_squared'])
		results[condition]['imag_eccen'] = np.sin(results[condition]['ecc'])*np.arctanh(stats['r_squared'])
		results[condition]['real_fwhm'] = np.cos(results[condition]['fwhm'])*np.arctanh(stats['r_squared'])
		results[condition]['imag_fwhm'] = np.sin(results[condition]['fwhm'])*np.arctanh(stats['r_squared'])

	#########################################################################################################################################################################################################################
	#### Plot results
	#########################################################################################################################################################################################################################

	if plotbool * (stats['r_squared']>0):# (np.random.randint(10)<10):#* (stats['r_squared']>0.1):#(stats['r_squared']>0.1):# * :# :#* (results['ecc'] < 3) :#:# * * randint ) #* :#* )

		n_TRs = n_TRs[0]
		n_runs = int(len(time_course) / n_TRs)
		if fit_on_all_data:
			plot_conditions = ['Ridge',conditions[0]]
		else:
			plot_conditions = conditions + ['All_fit']
		plot_dir = os.path.join(plotdir, '%s'%roi)
		if not os.path.isdir(plot_dir): os.mkdir(plot_dir)

		f=pl.figure(figsize=(20,8)); rowi = (n_runs+4)

		import colorsys
		colors = np.array([colorsys.hsv_to_rgb(c,0.6,0.9) for c in np.linspace(0,1,3+1)])[:-1]

		for runi in range(n_runs):
			s = f.add_subplot(rowi,1,runi+1)
			pl.plot(time_course[n_TRs*runi:n_TRs*(runi+1)],'-ok',linewidth=0.75,markersize=2.5)#,label='data'
			if not fit_on_all_data:
				for ci, condition in enumerate(conditions):
					pl.plot(model_predictions[condition][n_TRs*runi:n_TRs*(runi+1)]+params['baseline'].value,color=colors[ci],label='%s model'%condition,linewidth=2)				
				pl.plot([0,n_TRs],[params['baseline'].value,params['baseline'].value],color=colors[0],linewidth=1)	
			else:
				pl.plot(combined_model_prediction[n_TRs*runi:n_TRs*(runi+1)],color=colors[0],label='model',linewidth=2)	
			sn.despine(offset=10)
			pl.xlim(0,n_TRs*1.1)
			if runi == (n_runs-1):
				pl.xlabel('TRs')
			else:
				pl.xticks([])
			if runi == (n_runs/2):
				pl.legend(loc='best',fontsize=8)
				if 'psc' in postFix:
					pl.ylabel('% signal change')
				else:
					pl.ylabel('unkown unit')	
			pl.yticks([int(np.min(time_course)),0,int(np.max(time_course))])	
			pl.ylim([int(np.min(time_course)),int(np.max(time_course))])


		rowi = (n_runs+2)/2
		k = 0
		for ci, condition in enumerate(plot_conditions):
			k+= 1
			s = f.add_subplot(rowi,len(plot_conditions)*2,(rowi-1)*len(plot_conditions)*2+k,aspect='equal')
			pl.imshow(PRFs[condition],origin='lowerleft',interpolation='nearest',cmap=cm.coolwarm)

			pl.axis('off')
			s.set_title('%s PRF'%condition)
			
			k+= 1
			if not (condition == 'Ridge') + (condition == 'All_fit'):
				s = f.add_subplot(rowi,len(plot_conditions)*2,(rowi-1)*len(plot_conditions)*2+k)
				pl.imshow(np.ones((n_pixel_elements_raw,n_pixel_elements_raw)),cmap='gray')
				pl.clim(0,1)
				if model == 'OG':
					s.text(n_pixel_elements_raw/2,n_pixel_elements_raw/2, "\n%s PARAMETERS: \n\nbaseline: %.2f\nsize: %.2f\namplitude: %.6f\n\n\nDERIVED QUANTIFICATIONS: \n\nr-squared: %.2f\necc: %.2f\nFWHM: %.2f"%
						(condition,results[condition]['baseline'],results[condition]['sigma_center'],results[condition]['amp_center'],
							stats['r_squared'],results[condition]['ecc'],results[condition]['fwhm']),
						horizontalalignment='center',verticalalignment='center',fontsize=10,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
				elif model == 'DoG':
					s.text(n_pixel_elements_raw/2,n_pixel_elements_raw/2, "\n%s PARAMETERS: \n\nbaseline: %.2f\nsd center: %.2f\nsd surround: %.2f\namp center: %.6f\namp surround: %.6f\n\nDERIVED QUANTIFICATIONS: \n\nr squared: %.2f\necc: %.2f\nFWHM: %.2f\nsurround size: %.2f\nsupression index: %.2f"
						%(condition,results[condition]['baseline'],results[condition]['sigma_center'],results[condition]['sigma_surround'],results[condition]['amp_center'],
						results[condition]['amp_surround'],stats['r_squared'],results[condition]['ecc'],results[condition]['fwhm'],results[condition]['surround_size'],
						results[condition]['SI']),horizontalalignment='center',verticalalignment='center',fontsize=10,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
				pl.axis('off')

		# pl.tight_layout()
		pl.savefig(os.path.join(plot_dir, 'vox_%d_%d_%d.pdf'%(slice_no,voxno,n_pixel_elements_raw)))
		pl.close()

	return results, stats


	
