import numpy as np
import itertools
import scipy as sp
import os


def get_adjpars_grid(file,delimeter=' '):
    '''
    This function reads in a generic file with a given delimeter of a format
    Parameter - lower limit - upper limit - step size - format multiplier - string format - position index
    '''
    adjpars = {}
    for line in open(file,'r'):
        if line.strip():
            vals = line.split(delimeter)
            par = vals[0]
            vals = vals[1:]
            adjpars[par.strip()]={'values': sp.arange(float(vals[0]),float(vals[1]),float(vals[2]),dtype=np.float64), 
                                  'multiplier':int(vals[3]), 'zfill':int(vals[4]), 
                                  'index':int(vals[5]), 'ftype':vals[6] }

    return adjpars



def make_filenames(adjpars, summary_file='./grid.summary', inlist_prefix=None, grid_prefix=None, grid_type=None, worker_prefix=None):
    parts              = []
    combinations       = []
    directory_segments = []
    keys               = adjpars.keys()
    indices            = [ adjpars[key]['index'] for key in keys ]
    sortkeys           = np.array(sorted(zip(keys,indices),key=lambda x: x[1]))

    history_array      = []
    inlist_array       = []
    combinations_array = []
    logdirs_array      = []
    point_name_array   = []

    summary_file_open = open(summary_file,'w')

    for ii,key in enumerate(sortkeys[:,0]):

        parts.append(  [ key+str( int(adjpars[key]['multiplier']*value) ).zfill(adjpars[key]['zfill']) for value in adjpars[key]['values'] ] )
        combinations.append(  adjpars[key]['values'] )

        if 'D' in adjpars[key]['ftype']:
            directory_segments.append( key )

    list_all = []
    list_tmp = []
    list_str = []


    for combination in itertools.product(*combinations):
        str_combination = [str(cmbntn) for cmbntn in combination]
        list_tmp.append(np.hstack(combination))
        list_str.append(np.hstack(str_combination))

    for ii,combination in enumerate(itertools.product(*parts)):
        stack         = np.hstack(combination)

        if grid_type is not None: history_name = grid_type+'_'
        else: history_name = 'gridpoint_'

        if grid_type is not None: inlist_name = 'inlist_'+grid_type+'_'
        else: inlist_name   = 'inlist_'

        point_name    = ''

        if grid_prefix is not None: dir_name = grid_prefix
        else: dir_name = './'

        if worker_prefix is not None: worker_name = worker_prefix + 'worker-'
        else: worker_name ='./worker-'


        #for cc,dir_seg in enumerate(directory_segments):
        #    dir_name += '%s%s/'%(dir_seg,stack[cc])
        #    worker_name += '%s%s_'%(dir_seg,stack[cc])
        for cc,dir_key in enumerate(sortkeys[:,0]):
            if 'D' in adjpars[dir_key]['ftype']:
                dir_name += '%s/'%(stack[cc])
                worker_name += '%s_'%(stack[cc])

        worker_name.strip('_')


        print 'DIR: ',dir_name
        logdirs_array.append(dir_name)

        ## Make histories dir

        if not os.path.isdir('%shistories'%dir_name):
            os.makedirs('%shistories'%dir_name)
        else:
            print '%shistories EXISTS'%dir_name

        if not os.path.isdir('%sprofiles'%dir_name):
            os.makedirs('%sprofiles'%dir_name)
        else:
            print '%sprofiles EXISTS'%dir_name

        if not os.path.isdir('%sfrequencies'%dir_name):
            os.makedirs('%sfrequencies'%dir_name)
        else:
            print '%sfrequencies EXISTS'%dir_name

        if not os.path.isdir('%smodels'%dir_name):
            os.makedirs('%smodels'%dir_name)
        else:
            print '%smodels EXISTS'%dir_name

        '''
        try: 
            os.system('mkdir -p %shistories'%dir_name)
        except:
            print 'DIR: %shistories already exists'%dir_name
        ## Make profiles dir
        try: 
            os.system('mkdir -p %sprofiles'%dir_name)
        except:
            print 'DIR: %sprofiles already exists'%dir_name
        ## Make frequencies dir
        try: 
            os.system('mkdir -p %sfrequencies'%dir_name)
        except:
            print 'DIR: %sfrequencies already exists'%dir_name
        ## Make tams_models dir
        try: 
            os.system('mkdir -p %stams_models'%dir_name)
        except:
            print 'DIR: %stams_models already exists'%dir_name
        '''

        for jj,string in enumerate(stack):
            history_name  += string+'_'
            inlist_name   += string+'_'
            point_name    += string+'_'

        if inlist_prefix  is not None: inlist_name  = inlist_prefix  + inlist_name

        point_name    = point_name.strip('_')
        inlist_name   = inlist_name.strip('_')
        history_name  = history_name.strip('_')
        history_name  = 'histories/' + history_name + '.history'

        line =  '%s\t%s\t%s %.5f %.5f %.5f %.5f %.5f\n'%(dir_name,inlist_name,point_name,list_tmp[ii][0],
                                                   list_tmp[ii][1],list_tmp[ii][2],list_tmp[ii][3],list_tmp[ii][4])

        summary_file_open.write(line)

	worker_line = '%s,%s,%s\n'%(inlist_name,dir_name,point_name)
        worker_file_open = open(worker_name,'a')
        worker_file_open.write(worker_line)
        worker_file_open.close()        

        inlist_array.append(inlist_name)
        history_array.append(history_name)
        combinations_array.append(list_tmp[ii])
        point_name_array.append(point_name)

    summary_file_open.close()
    return logdirs_array,inlist_array,history_array,combinations_array,point_name_array


def write_mesa_inlist_vsc_grid(  inlist, history_file, logsdir, 
                                 initial_Y, initial_Z, initial_mass,
                                 overshoot_f, mlt_alpha, log_minDmix,
                                 xctrl_1 = 0.705, xctrl_2 = 0.9975,
                                 tams_model_file = 'tams.model',
                                 base_inlist=os.path.expandvars('$VSC_DATA/python/mesa/inlist_MAMSIE_BASE_VSCGRID')):

        print 'WRITING INLIST: {}'.format(inlist)
        with open(base_inlist, 'r') as f:
                lines = f.readlines()
        replacements = {
		'TAMS_MODEL_FILENAME'           : '{}'.format("'"+tams_model_file+"'"),
                'STAR_HISTORY_NAME'             : '{}'.format("'"+history_file+"'"),
                'LOGS_DIR'                      : '{}'.format("'"+logsdir+"'"),
                'INITIAL_Y'                     : '{:6.5f}d0'.format(initial_Y),
                'INITIAL_Z'                     : '{:6.5f}d0'.format(initial_Z),
                'INITIAL_MASS'                  : '{:8.6f}d0'.format(initial_mass),
                'OVERSHOOT_F_ABOVE_BURN_H_CORE' : '{:5.4f}d0'.format(overshoot_f),
                'MIXING_LENGTH_ALPHA'           : '{:5.4f}d0'.format(mlt_alpha),
                'MIN_D_MIX'                     : '{:8.4f}d0'.format(10**log_minDmix),
                'X_CTRL_1'                      : '{:5.4f}d0'.format(xctrl_1),
                'X_CTRL_2'                      : '{:5.4f}d0'.format(xctrl_2),
                        }

        new_lines = []
        for line in lines:
                new_line = line
                for key in replacements:
                        if (replacements[key] != ''):
                                new_line = new_line.replace(key, replacements[key])
                new_lines.append(new_line)

        with open(inlist, 'w') as f:
                f.writelines(new_lines)
        return
