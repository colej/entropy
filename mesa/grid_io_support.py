import numpy as np
import itertools
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
            adjpars[par.strip()]={'values': np.arange(float(vals[0]),float(vals[1]),float(vals[2])), 
                                  'multiplier':float(vals[3]), 'zfill':int(vals[4]), 
                                  'index':int(vals[5]), 'ftype':vals[6] }
    return adjpars



def make_filenames(adjpars, summary_file='./grid.summary', inlist_prefix=None, grid_prefix=None):
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

    summary_file_open = open(summary_file,'w')

    for ii,key in enumerate(sortkeys[:,0]):

        parts.append(  [ str( int(adjpars[key]['multiplier']*value) ).zfill(adjpars[key]['zfill']) for value in adjpars[key]['values'] ] )
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
        stack = np.hstack(combination)
        history_name  = 'gridpoint_'
        inlist_name  = 'inlist_'

        if grid_prefix is not None: dir_name = grid_prefix
        else: dir_name = './'

        for cc,dir_seg in enumerate(directory_segments):
            #dir_name += '%s%s/'%(dir_seg,list_str[ii][cc])
            dir_name += '%s%s/'%(dir_seg,stack[cc])
        print 'DIR: ',dir_name
        logdirs_array.append(dir_name)
        try: 
            os.system('mkdir -p %shistories'%dir_name)
        except:
            print 'DIR: %shistories already exists'%dir_name

        try: 
            os.system('mkdir -p %sprofiles'%dir_name)
        except:
            print 'DIR: %sprofiles already exists'%dir_name

        try: 
            os.system('mkdir -p %sfrequencies'%dir_name)
        except:
            print 'DIR: %sfrequencies already exists'%dir_name


        for jj,string in enumerate(stack):
            history_name  += string+'_'
            inlist_name   += string+'_'

        #if history_prefix is not None: history_name = history_prefix + history_name
        #history_name = dir_name + history_name
        if inlist_prefix  is not None: inlist_name  = inlist_prefix  + inlist_name


        inlist_name   = inlist_name.strip('_')
        history_name  = history_name.strip('_')
        history_name += '.history'

        line =  '%s\t%s\t%s %.5f %.5f %.5f %.5f %.5f\n'%(dir_name,inlist_name,history_name,list_tmp[ii][0],
                                                   list_tmp[ii][1],list_tmp[ii][2],list_tmp[ii][3],list_tmp[ii][4])
        summary_file_open.write(line)
        
        inlist_array.append(inlist_name)
        history_array.append(history_name)
        combinations_array.append(list_tmp[ii])

    return logdirs_array,inlist_array,history_array,combinations_array


def write_mesa_inlist_vsc_grid(  inlist, history_file, logsdir,
                                 initial_Y, initial_Z, initial_mass,
                                 overshoot_f, mlt_alpha, log_minDmix,
                                 xctrl_1 = 0.705, xctrl_2 = 0.9975,
                                 base_inlist='/home/cole/python/mesa/inlist_MAMSIE_BASE_VSCGRID'):

        print 'WRITING INLIST: {}'.format(inlist)
        with open(base_inlist, 'r') as f:
                lines = f.readlines()
        replacements = {
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
    
