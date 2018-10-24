import matplotlib.pyplot as plt
import numpy as np
import os,logging
import read_mesa as rm
import subprocess
from scipy.interpolate import interp1d


def corot_to_inertial(f_c,f_rot,azimuthal_order):
        return f_c + azimuthal_order*f_rot


def write_gyre_in( gyre_in_file, mesa_pulsation_file, gyre_out_file,
                   freq_min_inertial=-1.0, freq_max_inertial=-1.0,
                   npg_min=-120,npg_max=-5, omega_rot=-1.0, degrees=[1], azimuthal_orders=[1],
                   gyre_base_file = os.path.expandvars('$VSC_DATA/python/mesa/GYRE_MAMSIE_BASE')
                  ):

        freq_min_corot = freq_min_inertial - azimuthal_orders[0] * omega_rot
        freq_max_corot = freq_max_inertial - azimuthal_orders[0] * omega_rot

        with open(gyre_base_file, 'r') as f:
                lines = f.readlines()
        replacements = {
                        'FILENAME': '{}'.format("'"+mesa_pulsation_file+"'"),
                        'OUTPUT'  : '{}'.format("'"+gyre_out_file+"'"),
                        'N_PG_MIN': '{:1.0f}'.format(npg_min),
                        'N_PG_MAX': '{:1.0f}'.format(npg_max),
                        'FREQ_MIN': '{:8.6f}'.format(max(freq_min_corot,0.1)),
                        'FREQ_MAX': '{:8.6f}'.format(freq_max_corot),
                        'VROT'    : '{:12.10}'.format(omega_rot)
                        }

        new_lines = []
        for line in lines:
                new_line = line
                for key in replacements:
                        if (replacements[key] != ''):
                                new_line = new_line.replace(key, replacements[key])
                new_lines.append(new_line)

        with open(gyre_in_file, 'w') as f:
                f.writelines(new_lines)
        return


def run_gyre(gyre_in_file, gyre_out_file ):

        print '\tIN RUN GYRE:'
        gyre_com = [os.path.expandvars('$GYRE_DIR/bin/gyre') , gyre_in_file ]
        print '\t\tGYRE_COM: ',gyre_com
        gyre_proc = subprocess.Popen(gyre_com,bufsize=-1,stdout=subprocess.PIPE ,stderr= subprocess.PIPE)
        gyre_out = gyre_proc.stdout.readlines()

        if gyre_proc.returncode is not None:
                print '\t\t\t\tProcess Return code: ',gyre_proc.returncode
                try:
                        gyre_proc.kill()
                except:
                        print 'couldnt kill gyre process'
                return -np.inf,None
        else:
                try:
                        gyre_proc.kill()
                        print '\t\tProcess Terminated'
                except:
                        print 'Gyre Proc alrady terminated'


        ## Read frequencies, n_pgs and translate back to inertial frame
        gyre_data      = np.genfromtxt(gyre_out_file,skip_header=5,names=True)
        freqs_inertial = gyre_data['Refreq']
        n_g            = gyre_data['n_g']

        return freqs_inertial,n_g


def generate_obs_series(periods,errors):
    observed_spacings        = []
    observed_spacings_errors = []

    for kk,prd_k in enumerate(periods[:-1]):
        prd_k_p_1 = periods[kk+1]
        observed_spacings.append( abs( prd_k - prd_k_p_1 )*86400. )
        observed_spacings_errors.append(np.sqrt( errors[kk]**2 + errors[kk+1]**2  )*86400.)
    return observed_spacings,observed_spacings_errors


def generate_thry_series(periods):

    theoretical_spacings = []

    for kk,prd_k in enumerate(periods[:-1]):
        prd_k_p_1 = periods[kk+1]
        theoretical_spacings.append( abs(prd_k-prd_k_p_1)*86400. )
    return theoretical_spacings


def chisq_period_series_interp(tperiods,tspacings,operiods,ospacings,ospacing_errors):
    interp = interp1d(tperiods[1:][::-1],np.array(tspacings[::-1]),kind='linear')(np.array(operiods[1:][::-1]))[::-1]
    chi = (np.array(interp)-np.array(ospacings))**2 / np.array(ospacing_errors)**2

    return np.sum(chi), interp


def chisq_experimental(tperiods,tspacings,orders,operiods,operiod_errors,ospacings,ospacing_errors):

	chisqs = []
	for ii,per in enumerate(operiods):
		chisqs.append([])
		for jj,tper in enumerate(tperiods):
			chi2 = (per-tper)**2 / operiod_errors[ii]**2
			chisqs[ii].append(chi2)

	for ii,per in enumerate(operiods):
		ind = np.where(chisqs[ii]==min(chisqs[ii]))[0]
		#print ind
		print 'Obs/Thr/chi2/order: %f\t%f\t%f\t%i'%(per,tperiods[ind],chisqs[ii][ind],orders[ind])


def chisq_longest_sequence_SAVE(tperiods,orders,operiods,operiods_errors):

        # Generate two series
        dP,e_dP = generate_obs_series(operiods,operiods_errors)
        deltaP  = generate_thry_series(tperiods)

        # Find the best matches per observed period
        pairs_orders = []
        for ii,period in enumerate(operiods):

                chisqs = np.array([ ( (period-tperiod)/operiods_errors[ii] )**2 for tperiod in tperiods  ])

		## Locate the theoretical frequency (and accompanying order) with the best chi2
                min_ind = np.where( chisqs == min( chisqs ) )[0]
                best_match = tperiods[min_ind][0]
                best_order = orders[min_ind][0]

		## Toss everything together for bookkeeping
                pairs_orders.append([period,best_match,int(best_order),chisqs[min_ind]])

	pairs_orders = np.array(pairs_orders)

	#plt.figure(1,figsize=(6.6957,6.6957))
	#plt.subplot(211)
	#plt.plot(pairs_orders[:,0],pairs_orders[:,1],'o')
	#plt.ylabel('$\\mathrm{Period \\,[d]}$',fontsize=20)
	#plt.subplot(212)
	#plt.plot(pairs_orders[:,0],pairs_orders[:,2],'o')
	#plt.ylabel('$\\mathrm{Radial \\, Order}$',fontsize=20)
	#plt.xlabel('$\\mathrm{Period \\,[d]}$',fontsize=20)

	#plt.show()


	sequences = []
	## Look through all pairs of obs and theoretical frequencies and
	## check if the next obs freqency has a corresponding theoretical frequency
	## with the consecutive radial order
	current = []
        for ii,sett in enumerate(pairs_orders[:-1]):
                if abs(sett[2]) == abs(pairs_orders[ii+1][2])+1:
			#print sett,'-->',pairs_orders[ii+1]
			current.append(sett)
		else:
			current.append(sett)
			#print 'BREAK', np.array(current)[:,2]
			sequences.append(np.array(current).reshape(len(current),4))
			current = []

	len_list = np.array([len(x) for x in sequences])
	longest = np.where(len_list == max(len_list))[0] #[0]

	## Test if there really is one longest sequence
	if len(longest) == 1:
		lseq = sequences[longest[0]]

	## if not, pick, of all the sequences with the same length, the best based on chi2
	else:
		scores = [ np.sum(sequences[ii][:,-1])/len(sequences[ii]) for  ii in longest]
		min_score = np.where(scores == min(scores))[0][0]
		lseq = sequences[longest[min_score]]

	obs_ordering_ind = np.where(operiods == lseq[:,0][0])[0][0]
	thr_ordering_ind = np.where(tperiods == lseq[:,1][0])[0][0]

	thr_ordering_start = thr_ordering_ind - len(operiods[:obs_ordering_ind])
	thr_ordering_stop  = thr_ordering_ind + len(operiods[obs_ordering_ind:])

	ordered_theoretical_periods   = []
	corresponding_orders          = []

	for ii,oper in enumerate(operiods[:obs_ordering_ind]):
		tper  = tperiods[thr_ordering_start+ii]
		order = orders[thr_ordering_start+ii]
		#ordered_theoretical_periods_a.append(tper)
		ordered_theoretical_periods.append(tper)
		corresponding_orders.append(order)

	for ii,oper in enumerate(operiods[obs_ordering_ind:]):
		tper  = tperiods[thr_ordering_ind+ii]
		order = orders[thr_ordering_ind+ii]
		#ordered_theoretical_periods_b.append(tper)
		ordered_theoretical_periods.append(tper)
		corresponding_orders.append(order)

	#final_theoretical_periods = np.sort(np.hstack([ordered_theoretical_periods_a,ordered_theoretical_periods_b]))[::-1]
	final_theoretical_periods = np.array(ordered_theoretical_periods)

	obs_series,obs_series_errors = generate_obs_series(operiods,operiods_errors)
	thr_series = generate_thry_series(final_theoretical_periods)

	obs_series        = np.array(obs_series)
	obs_series_errors = np.array(obs_series_errors)
	thr_series        = np.array(thr_series)

	series_chi2 = np.sum( (obs_series-thr_series)**2 /obs_series_errors**2 ) / len(obs_series)
    # print 'orders: %i - %i'%(corresponding_orders[0],corresponding_orders[-1])

	fig = plt.figure(1,figsize=(6.6957,6.6957))
	fig.suptitle('$\mathrm{Longest \\ Sequence}$',fontsize=20)
	axT = fig.add_subplot(211)
	axT.errorbar(operiods[1:],obs_series,yerr=obs_series_errors,marker='x',color='black',label='Obs')
	axT.plot(final_theoretical_periods[1:],thr_series,'rx-',label='Theory')
	axT.set_ylabel('$\mathrm{Period \\ Spacing \\ (s)}$',fontsize=20)
	axT.legend(loc='best')
	axB = fig.add_subplot(212)
	axB.errorbar(operiods[1:],obs_series-thr_series,yerr=obs_series_errors,marker='',color='black')
	axB.set_ylabel('$\mathrm{Residuals \\ (s)}$',fontsize=20)
	axB.set_xlabel('$\mathrm{Period \\ (d^{-1})}$',fontsize=20)
	axB.text(0.75,0.85,'$\chi^2 = %.2f$'%series_chi2,fontsize=15,transform=axB.transAxes)

	plt.show()

	series_chi2 = np.sum( ( (obs_series-thr_series) /obs_series_errors )**2 ) / len(obs_series)
	return series_chi2,final_theoretical_periods,corresponding_orders


def chisq_longest_sequence(tperiods,orders,operiods,operiods_errors):
    if len(tperiods)<len(operiods):
        return 1e16, [-1. for i in range(len(operiods))], [-1 for i in range(len(operiods))]
    else:
        # Generate two series
        dP,e_dP = generate_obs_series(operiods,operiods_errors)
        deltaP  = generate_thry_series(tperiods)

        # Find the best matches per observed period
        pairs_orders = []
        for ii,period in enumerate(operiods):
            chisqs = np.array([ ( (period-tperiod)/operiods_errors[ii] )**2 for tperiod in tperiods  ])

            ## Locate the theoretical frequency (and accompanying order) with the best chi2
            min_ind = np.where( chisqs == min( chisqs ) )[0]
            best_match = tperiods[min_ind][0]
            best_order = orders[min_ind][0]

            ## Toss everything together for bookkeeping
            pairs_orders.append([period,best_match,int(best_order),chisqs[min_ind]])

        pairs_orders = np.array(pairs_orders)

        plt.figure(1,figsize=(6.6957,6.6957))
        plt.subplot(211)
        plt.plot(pairs_orders[:,0],pairs_orders[:,1],'o')
        plt.ylabel('$\\mathrm{Period \\,[d]}$',fontsize=20)
        plt.subplot(212)
        plt.plot(pairs_orders[:,0],pairs_orders[:,2],'o')
        plt.ylabel('$\\mathrm{Radial \\, Order}$',fontsize=20)
        plt.xlabel('$\\mathrm{Period \\,[d]}$',fontsize=20)

        # plt.show()


        sequences = []
        ## Look through all pairs of obs and theoretical frequencies and
        ## check if the next obs freqency has a corresponding theoretical frequency
        ## with the consecutive radial order
        current = []
        lp = len(pairs_orders[:-1])
        for ii,sett in enumerate(pairs_orders[:-1]):
                if abs(sett[2]) == abs(pairs_orders[ii+1][2])+1:
                        current.append(sett)
                else:
                       	current.append(sett)
                        sequences.append(np.array(current).reshape(len(current),4))
                        current = []
                if (ii==lp-1):
                     	current.append(sett)
                        sequences.append(np.array(current).reshape(len(current),4))
                        current = []
        len_list = np.array([len(x) for x in sequences])
        longest = np.where(len_list == max(len_list))[0] #[0]

        ## Test if there really is one longest sequence
        if len(longest) == 1:
            lseq = sequences[longest[0]]

        ## if not, pick, of all the sequences with the same length, the best based on chi2
        else:
            scores = [ np.sum(sequences[ii][:,-1])/len(sequences[ii]) for  ii in longest]
            min_score = np.where(scores == min(scores))[0][0]
            lseq = sequences[longest[min_score]]

        obs_ordering_ind = np.where(operiods == lseq[:,0][0])[0][0]
        thr_ordering_ind = np.where(tperiods == lseq[:,1][0])[0][0]

        ordered_theoretical_periods   = []
        corresponding_orders          = []

        thr_ind_start = thr_ordering_ind - obs_ordering_ind
        thr_ind_current = thr_ind_start

        for i,oper in enumerate(operiods):
            thr_ind_current = thr_ind_start + i
            if (thr_ind_current < 0):
                tper = -1
                ordr = -1
            elif (thr_ind_current >= len(tperiods)):
                tper = -1
                ordr = -1
            else:
                tper = tperiods[thr_ind_current]
                ordr = orders[thr_ind_current]
            ordered_theoretical_periods.append(tper)
            corresponding_orders.append(ordr)

        #final_theoretical_periods = np.sort(np.hstack([ordered_theoretical_periods_a,ordered_theoretical_periods_b]))[::-1]
        final_theoretical_periods = np.array(ordered_theoretical_periods)

        obs_series,obs_series_errors = generate_obs_series(operiods,operiods_errors)
        thr_series = generate_thry_series(final_theoretical_periods)

        obs_series        = np.array(obs_series)
        obs_series_errors = np.array(obs_series_errors)
        thr_series        = np.array(thr_series)

        series_chi2 = np.sum( (obs_series-thr_series)**2 /obs_series_errors**2 ) / len(obs_series)
        # print 'orders: %i - %i'%(corresponding_orders[0],corresponding_orders[-1])

        fig = plt.figure(2,figsize=(6.6957,6.6957))
        fig.suptitle('$\mathrm{Longest \\ Sequence}$',fontsize=20)
        axT = fig.add_subplot(211)
        # axT.errorbar(operiods[1:],obs_series,yerr=obs_series_errors,marker='x',color='black',label='Obs')
        # axT.plot(final_theoretical_periods[1:],thr_series,'rx-',label='Theory')
        axT.errorbar(range(len(obs_series)),obs_series,yerr=obs_series_errors,marker='x',color='black',label='Obs')
        axT.plot(range(len(thr_series)),thr_series,'rx-',label='Theory')
        axT.set_ylabel('$\mathrm{Period \\ Spacing \\ (s)}$',fontsize=20)
        axT.legend(loc='best')
        axB = fig.add_subplot(212)
        axB.errorbar(operiods[1:],obs_series-thr_series,yerr=obs_series_errors,marker='',color='black')
        axB.set_ylabel('$\mathrm{Residuals \\ (s)}$',fontsize=20)
        axB.set_xlabel('$\mathrm{Period \\ (d^{-1})}$',fontsize=20)
        axB.text(0.75,0.85,'$\chi^2 = %.2f$'%series_chi2,fontsize=15,transform=axB.transAxes)

        plt.show()

        for ii,oper in enumerate(operiods):
            print oper, final_theoretical_periods[ii], corresponding_orders[ii]

        series_chi2 = np.sum( ( (obs_series-thr_series) /obs_series_errors )**2 ) / len(obs_series)
        return series_chi2,final_theoretical_periods,corresponding_orders


def chisq_best_sequence(tperiods,orders,operiods,operiods_errors):

        # Generate two series
        dP,e_dP = generate_obs_series(operiods,operiods_errors)
        deltaP  = generate_thry_series(tperiods)

        # Find the best matches per observed period
        pairs_orders = []
        for ii,period in enumerate(operiods):

                chisqs = np.array([ ( (period-tperiod)/operiods_errors[ii] )**2 for tperiod in tperiods  ])

		## Locate the theoretical frequency (and accompanying order) with the best chi2
                min_ind = np.where( chisqs == min( chisqs ) )[0]#[0]
                best_match = tperiods[min_ind][0]
                best_order = orders[min_ind][0]

		## Toss everything together for bookkeeping
                pairs_orders.append([period,best_match,int(best_order),chisqs[min_ind]])

	pairs_orders = np.array(pairs_orders)

	sequences = []
	## Look through all pairs of obs and theoretical frequencies and
	## check if the next obs freqency has a corresponding theoretical frequency
	## with the consecutive radial order
	current = []
        for ii,sett in enumerate(pairs_orders[:-1]):
                if abs(sett[2]) == abs(pairs_orders[ii+1][2])+1:
			#print sett,'-->',pairs_orders[ii+1]
			current.append(sett)
		else:
			current.append(sett)
			#print 'BREAK', np.array(current)[:,2]
			sequences.append(np.array(current).reshape(len(current),4))
			current = []


	scores = [ np.sum(sequence[:,-1])/len(sequence) for sequence in sequences]
	min_score = np.where(scores == min(scores))[0][0]
	lseq = sequences[min_score]

	obs_ordering_ind = np.where(operiods == lseq[:,0][0])[0][0]
	thr_ordering_ind = np.where(tperiods == lseq[:,1][0])[0][0]

	thr_ordering_start = thr_ordering_ind - len(operiods[:obs_ordering_ind])
	thr_ordering_stop  = thr_ordering_ind + len(operiods[obs_ordering_ind:])

	ordered_theoretical_periods   = []
	corresponding_orders          = []

	for ii,oper in enumerate(operiods[:obs_ordering_ind]):
		tper  = tperiods[thr_ordering_start+ii]
		order = orders[thr_ordering_start+ii]
		ordered_theoretical_periods.append(tper)
		corresponding_orders.append(order)

	for ii,oper in enumerate(operiods[obs_ordering_ind:]):
		tper  = tperiods[thr_ordering_ind+ii]
		order = orders[thr_ordering_ind+ii]
		ordered_theoretical_periods.append(tper)
		corresponding_orders.append(order)

	final_theoretical_periods = np.array(ordered_theoretical_periods)


	obs_series,obs_series_errors = generate_obs_series(operiods,operiods_errors)
	thr_series = generate_thry_series(final_theoretical_periods)

	obs_series        = np.array(obs_series)
	obs_series_errors = np.array(obs_series_errors)
	thr_series        = np.array(thr_series)
	series_chi2 = np.sum( (obs_series-thr_series)**2 /obs_series_errors**2 ) / len(obs_series)



	fig = plt.figure(1,figsize=(6.6957,6.6957))
	fig.suptitle('$\mathrm{Best \\ \chi^2}$',fontsize=20)
	axT = fig.add_subplot(211)
	axT.errorbar(operiods[1:],obs_series,yerr=obs_series_errors,marker='x',color='black',label='Obs')
	axT.plot(final_theoretical_periods[1:],thr_series,'rx-',label='Theory')
	axT.set_ylabel('$\mathrm{Period \\ Spacing \\ (s)}$',fontsize=20)
	axT.legend(loc='best')
	axB = fig.add_subplot(212)
	axB.errorbar(operiods[1:],obs_series-thr_series,yerr=obs_series_errors,marker='',color='black')
	axB.set_ylabel('$\mathrm{Residuals \\ (s)}$',fontsize=20)
	axB.set_xlabel('$\mathrm{Period \\ (d^{-1})}$',fontsize=20)
	axB.text(0.75,0.85,'$\chi^2 = %.2f$'%series_chi2,fontsize=15,transform=axB.transAxes)

	plt.show()

	return series_chi2,final_theoretical_periods,corresponding_orders


def chisq_period_series_pairwise(tperiods,tspacings,orders,operiods,operiod_errors,ospacings,ospacing_errors):

	print tperiods
	tpairs     = []
	tpairs_ind = []
	for kk,period in enumerate(operiods):
		chi2 = [ ( (tper-period)/operiod_errors[kk])**2 for tper in tperiods ]
		print min(chi2),period,tperiods[np.where(chi2 == min(chi2) )[0]]
		tpairs.append( tperiods[np.where( chi2 == min(chi2) )][0] )
		tpairs_ind.append(      np.where( chi2 == min(chi2) )[0]  )
	tpairs_ind = np.hstack([pind for pind in tpairs_ind])

	count_lens  = []
	start_stops = []
	ncounter    = 0
	start,stop  = 0,0
	print orders[tpairs_ind]
	for ii,n in enumerate(orders[tpairs_ind][:-1]):
		print abs(int(n)),int(abs(orders[tpairs_ind][ii+1])+1)
		print '\t',start,stop
		if int(abs(n)) == int(abs(orders[tpairs_ind][ii+1])+1):
			stop     += 1
			ncounter += 1
		else:
			#print start,stop
			stop +=1
			start_stops.append((start,stop))
			count_lens.append(ncounter)
			start,stop = stop,ii+1
			ncounter   = 0

	start_stops.append((start,stop))
	count_lens.append(ncounter)

	count_lens = np.array(count_lens)
	print count_lens

	best_ind     = np.where(count_lens == max(count_lens))[0]
	print best_ind
	bstart,bstop = start_stops[best_ind]
	best_seq     = operiods[bstart:bstop]
	best_thry    = tperiods[::-1][:]
	print best_thry
	print 'BEST SEQUNCE: '
	for jj,per in enumerate(best_seq[:-1]):
		print best_seq[jj+1],per,abs(best_seq[jj+1]-per)*86400.,best_thry[jj+1],best_thry[jj]







	plt.subplot(211)
	plt.ylabel('Period [s]',fontsize=20)
	plt.plot(operiods,tpairs,'ko')
	plt.plot(operiods,operiods,'r--')
	plt.subplot(212)
	plt.xlabel('Period [s]',fontsize=20)
	plt.ylabel('Order',fontsize=20)
	plt.plot(operiods,orders[tpairs_ind],'ko')

	segments = []
	for ii,cnt in enumerate(count_lens):
		if cnt>2:
			start,stop  = start_stops[ii]
			mark_orders = orders[tpairs_ind][start:stop]
			operiods=np.array(operiods)
			plt.plot(operiods[start:stop],mark_orders,'s',ms=5)
			segments.append(operiods[start:stop])
	series = []
	new_pers = []
	for seg in segments:
		for jj,per in enumerate(seg[:-1]):
			series.append(abs(per-seg[jj+1])*86400.)
			new_pers.append(seg[jj+1])
		#for per in seg:
		#	new_pers.append(per)
	#plt.show()

	#chi = (np.array(interp)-np.array(ospacings))**2 / np.array(ospacing_errors)**2
	return new_pers,series
	#return np.sum(chi), interp

#if __name__=='__main__':

    #obs_freqs = [1.37295,1.33454,1.29974,1.26772,1.240146,1.216652,1.194916,
                 #1.17327,1.15451,1.13626,1.12089,1.10739, 1.09304, 1.08028,
                 #1.06837,1.05692,1.04602,1.03778,1.02992, 1.02072]
    #obs_periods = [0.72836,0.74932,0.76938,0.78882,0.806357,0.821928,0.836879,
                   #0.85232,0.86617,0.88008,0.89215,0.90302 ,0.91488, 0.92569,
                   #0.93601,0.94614,0.95600,0.96360,0.97095 ,0.97970]
    #sig_periods = [3e-5,3e-5,2e-5,3e-5,4e-6,4e-6,4e-6,1e-5,3e-5,4e-5,5e-5,
                   #3e-5,4e-5,3e-5,3e-5,4e-5,5e-5,4e-5,6e-5,6e-5]

    #gfile = '/STER/colej/astro/dOB/kic493/single_work/gyre_test.out'
    #dP,sigma_dP,periods,Delta_P,theoretical_periods = generate_series(gfile,obs_periods, sig_periods)
    #interp = interp1d(theoretical_periods[1:][::-1],Delta_P[::-1],kind='linear')(periods[1:][::-1])
    ##print interp(periods[1:-1])
    ##print len(interp),len(periods)
    #plt.errorbar(periods[1:],dP,yerr=sigma_dP,marker='',color='red')
    #plt.plot(periods[1:],dP,'k-x')
    #plt.plot(theoretical_periods[1:],Delta_P,'r-^')
    #plt.plot(periods[1:],interp[::-1],'b*-')
    #plt.show()
