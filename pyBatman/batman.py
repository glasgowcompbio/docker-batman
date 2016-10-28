import datetime
import matplotlib
import numpy as np
import nmrglue as ng
import pandas as pd
import os
from IPython.display import display, HTML
from collections import OrderedDict

import multiprocessing
import json
import errno
import os
import tempfile
import shutil
from distutils.dir_util import copy_tree
from collections import namedtuple
import sys

from random import randint
import matplotlib.pyplot as plt
import pylab as plt
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
import plotly

import readline
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

READ_SPECTRA_NMRGLUE = 1
READ_SPECTRA_BATMAN = 2

class PyBatman(object):

    def __init__(self, input_dirs, pattern, working_dir, db, verbose=False,
        read_spectra_using=READ_SPECTRA_NMRGLUE):

        self.db = db
        self.input_dirs = input_dirs
        self.pattern = pattern
        self.working_dir = working_dir
        self.read_spectra_using = read_spectra_using
        self.verbose = verbose

        importr('batman')
        pandas2ri.activate() # to convert between R <-> pandas dataframes

        # see http://stackoverflow.com/questions/3925096/how-to-get-only-the-last-part-of-a-path-in-python
        self.labels = []
        for input_dir in self.input_dirs:
            last = os.path.basename(os.path.normpath(input_dir))
            self.labels.append(last)

        # extract spectra name
        self.matching = self._get_matching_paths(self.input_dirs,
                            self.pattern)
        if self.verbose:
            print "Found spectra matching pattern '%s':" % self.pattern
            for s in self.matching:
                print '-', s

        # create working directory
        self._mkdir_p(self.working_dir)
        prefix = datetime.datetime.now().strftime("%G%m%d_%H%M%S_")
        self.temp_dir = tempfile.mkdtemp(prefix=prefix, dir=self.working_dir)
        batman_dir = os.path.join(self.temp_dir, 'runBATMAN')
        self.batman_input = os.path.join(batman_dir, 'BatmanInput')
        self.batman_output = os.path.join(batman_dir, 'BatmanOutput')
        self._mkdir_p(self.batman_input)
        self._mkdir_p(self.batman_output)
        if self.verbose:
            print 'Working directory =', self.working_dir
            print '- batman_input =', self.batman_input
            print '- batman_output =', self.batman_output

        # extract spectra
        batman_data = os.path.join(batman_dir, 'data')
        temp_data = self._prepare_data(self.matching, batman_data)
        if self.read_spectra_using == READ_SPECTRA_NMRGLUE:
            self.spectra = self._load_data_with_nmrglue(temp_data)
        elif self.read_spectra_using == READ_SPECTRA_BATMAN:
            self.spectra = self._load_data_with_batman(batman_data)
        else:
            raise ValueError('Invalid value for options.read_spectra_using')

        # write spectra
        columns = ['ppm']
        columns.extend(self.labels)
        self.df = pd.DataFrame(self.spectra, columns=columns)
        self.spectra_file = os.path.join(self.batman_input, 'NMRdata_temp.txt')
        self.df.to_csv(self.spectra_file, index=False, sep='\t')
        if self.verbose:
            print 'Spectra written to', self.spectra_file

    def plot_spectra(self):
        plt.figure()
        x = self.spectra[:, 0]
        for n in range(1, self.spectra.shape[1]):
            y = self.spectra[:, n]
            plt.plot(x, y)
        plt.xlim([-0.1, 5])
        plt.gca().invert_xaxis()
        plt.xlabel('ppm')
        plt.ylabel('intensity')

    def plotly_spectra(self):
        data = []
        x = self.spectra[:, 0]
        for n in range(1, self.spectra.shape[1]):
            y = self.spectra[:, n]
            trace = go.Scatter(
                x = x,
                y = y
            )
            data.append(trace)
        plotly.offline.iplot(data)

    def get_default_params(self, names):

        # check if it's TSP
        is_tsp = False
        if 'TSP' in names:
            if len(names) > 1:
                raise ValueError('Not recommended to fit TSP with other metabolites')
            is_tsp = True

        # select only the metabolites we need
        selected = self._select_metabolites(names)
        options = PyBatmanOptions(selected)

        # and their ranges
        ranges = []
        for m in selected:
            if m.ppm_range > 0:
                ranges.extend(m.ppm_range) # a list of tuples
        ranges = ['(%s, %s)' % r for r in ranges] # convert to a list of string
        ranges_str = ' '.join(ranges)

        # set no. of spectra and no. of processors to use
        n_row, n_col = self.df.shape
        spec_no = n_col-1 # first column is the ppm
        para_proc = multiprocessing.cpu_count()
        random_seed = randint(0, 1e6)

        # set common parameters
        options = options.set('ppmRange', ranges_str)             \
                             .set('specNo', '1-%d' % spec_no)     \
                             .set('paraProc', para_proc)          \
                             .set('randSeed', random_seed)        \
                             .set('nItBurnin', 9000)              \
                             .set('nItPostBurnin', 1000)          \
                             .set('thinning', 5)                  \
                             .set('tauMean', -0.05)               \
                             .set('tauPrec', 2)                   \
                             .set('rdelta', 0.002)                \
                             .set('csFlag', 0)

        # special parameters for TSP since it's so different from the rest
        if is_tsp:
            default = options.set('muMean', 2)                    \
                                 .set('muVar', 0.1)               \
                                 .set('muVar_prop', 0.002)        \
                                 .set('nuMVar', 0)                \
                                 .set('nuMVarProp', 0.1)
        else:
            default = options.set('muMean', 0)                    \
                                 .set('muVar', 0.01)              \
                                 .set('muVar_prop', 0.0002)       \
                                 .set('nuMVar', 0.0025)           \
                                 .set('nuMVarProp', 0.01)

        return default

    def run(self, options, parallel=False, seed=None, plot=True):

        if not parallel:
            options = options.set('paraProc', 1)

        if seed is not None:
            options = options.set('randSeed', seed)

        # write out batman input files
        selected = options.selected
        metabolites_list_df = self._write_metabolites_list(selected, 'metabolitesList.csv')
        multi_data_user_df = self._write_multiplet_data(selected, 'multi_data_user.csv')
        chem_shift_per_spec_df = self._write_chem_shift(selected, 'chemShiftPerSpec.csv')
        options_lines = self._write_batman_options(options, 'batmanOptions.txt')

        # actually runs batman here using rpy2
        batman_r = robjects.r['batman']
        if self.verbose:

            display(metabolites_list_df)
            display(multi_data_user_df)
            display(chem_shift_per_spec_df)
            print
            print '----------------------------------------------------------------'
            print 'Parameters'
            print '----------------------------------------------------------------'
            for line in options_lines:
                print line
            print
            bm = batman_r(runBATMANDir=self.temp_dir, txtFile=self.spectra_file, figBatmanFit=False)

        else:

            sink_r = robjects.r['sink']
            sink_r('/dev/null')
            bm = batman_r(runBATMANDir=self.temp_dir, txtFile=self.spectra_file, figBatmanFit=False)
            sink_r()

        if plot:
            self._plot_batman_r(bm)
        output = BatmanOutput(bm)
        return output

    def cleanup(self):
        shutil.rmtree(self.temp_dir)
        if self.verbose:
            print 'Deleted', self.temp_dir

    def _select_metabolites(self, names):

        metabolite_names = [m.name for m in self.db.metabolites.values()]

        # activate only the metabolites we need
        self.db.deactivate_all()
        for name in names:
            assert name in metabolite_names
            self.db.activate(name)

        metabolites = self.db.get_active()
        return metabolites

    def _write_metabolites_list(self, metabolites, file_name):

        columns = ['Metabolite']
        data = [m.name for m in metabolites]
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(self.batman_input, file_name)
        if self.verbose:
            print 'metabolites list = %s' % out_file
        df.to_csv(out_file, header=False, index=False, sep=',')

        return df

    def _write_multiplet_data(self, metabolites, file_name):

        # create multi_data_user.csv
        columns = ['Metabolite', 'pos_in_ppm', 'couple_code', 'J_constant',
                   'relative_intensity', 'overwrite_pos', 'overwrite_truncation',
                   'Include_multiplet']
        data = []
        OVERWRITE_POS = 'n'
        OVERWRITE_TRUNCATION = 'n'
        INCLUDE_MULTIPLET = 1
        for m in metabolites:
            for u in m.multiplets:
                row = (u.parent.name, u.ppm, u.couple_code, u.j_constant, u.rel_intensity,
                        OVERWRITE_POS, OVERWRITE_TRUNCATION, INCLUDE_MULTIPLET)
                data.append(row)
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(self.batman_input, file_name)
        if self.verbose:
            print 'multiplet data = %s' % out_file
        df.to_csv(out_file, index=False, sep=',')

        return df

    def _write_chem_shift(self, metabolites, file_name):

        # create chemShiftPerSpec.csv
        columns = ['multiplets', 'pos_in_ppm']
        columns.extend(self.labels)
        data = []
        ns = ['n'] * len(self.labels)
        for m in metabolites:
            for u in m.multiplets:
                row = [u.parent.name, u.ppm]
                row.extend(ns)
                data.append(row)
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(self.batman_input, file_name)
        if self.verbose:
            print 'chem shift = %s' % out_file
        df.to_csv(out_file, index=False, sep=',')

        return df

    def _write_batman_options(self, options, file_name):

        # create batmanOptions.txt
        out_file = os.path.join(self.batman_input, file_name)
        if self.verbose:
            print 'batman options = %s' % out_file

        option_lines = []
        with open(out_file, 'w') as f:
            params = options.params
            params_desc = options.params_desc
            for key in params:
                value = str(params[key])
                description = params_desc[key]
                line = key + ' - ' + description + ': ' + value
                option_lines.append(line)
                f.write('%s\n' % line)

        return option_lines

    def _plot_batman_r(self, bm):

        plot_batman_fit_r = robjects.r['plotBatmanFit']
        plot_batman_fit_stack_r = robjects.r['plotBatmanFitStack']
        plot_rel_con_r = robjects.r['plotRelCon']
        plot_meta_fit_r = robjects.r['plotMetaFit']

        if self.verbose:
            plot_batman_fit_r(bm, showPlot=False)
            plot_batman_fit_stack_r(bm, offset=0.8, placeLegend='topleft', yto=5)
            plot_rel_con_r(bm, showPlot=False)
            plot_meta_fit_r(bm, showPlot=False)
        else:
            sink_r = robjects.r['sink']
            sink_r('/dev/null')
            plot_batman_fit_r(bm, showPlot=False)
            plot_batman_fit_stack_r(bm, offset=0.8, placeLegend='topleft', yto=5)
            plot_rel_con_r(bm, showPlot=False)
            plot_meta_fit_r(bm, showPlot=False)
            sink_r()

    def _get_matching_paths(self, input_dirs, pattern):

        matching = []
        for input_dir in input_dirs:

            # find all the child directories
            sub_dirs = [os.path.join(input_dir, x) for x in os.listdir(input_dir)]
            for path in sub_dirs:

                # check the pulse program if it exists
                if os.path.isdir(path):
                    pp = os.path.join(path, 'pulseprogram')
                    if os.path.isfile(pp): # if exists

                        # if it contains the pattern then store this path
                        with open(pp, 'r') as f:
                            head = [next(f) for x in xrange(2)]
                            if pattern in head[1]:
                                matching.append(path)

        return matching

    def _mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _prepare_data(self, matching, data_dir):
        out_dirs = []
        for f in range(len(matching)):
            path = matching[f]
            out_dir = os.path.join(data_dir, 'pdata', str(f))
            self._mkdir_p(out_dir)
            if self.verbose:
                print 'Copied spectra to', out_dir
            copy_tree(path, out_dir)
            out_dirs.append(out_dir)
        return out_dirs

    def _load_single_spectra(self, spectra_dir):

        p_data = os.path.join(spectra_dir, 'pdata/1')
        if self.verbose:
            print 'Processing', p_data

        dic, data = ng.bruker.read_pdata(p_data)
        udic = ng.bruker.guess_udic(dic, data)
        uc = ng.fileiobase.uc_from_udic(udic, 0)

        x = []
        y = []
        for ppm in uc.ppm_scale():
            x.append(ppm)
            y.append(data[uc(ppm, 'ppm')])
        x = np.array(x)
        y = np.array(y)

        return x, y

    def _load_data_with_nmrglue(self, input_dirs):
        spectra = []
        for input_dir in input_dirs:
            x, y = self._load_single_spectra(input_dir)
            spectra.append(y)
        combined = [x]
        combined.extend(spectra)
        combined = np.array(combined).transpose()
        if self.verbose:
            print 'Loaded', combined.shape
        return combined

    def _load_data_with_batman(self, input_dir):
        read_bruker = robjects.r['readBruker']
        spectra = read_bruker(input_dir)
        combined = np.array(spectra)
        if self.verbose:
            print 'Loaded', combined.shape
        return combined

class PyBatmanOptions(object):

    def __init__(self, selected, params=None):

        self.selected = selected
        if params is None: # default parameters

            # These parameters will be written to the batmanOptions.txt file, so
            # DO NOT CHANGE THE ORDER OF PARAMETERS HERE!!
            # since the parser in batman seems to break easily
            self.params = OrderedDict()
            self.params['ppmRange']      = None         # set inside PyBatman.get_default_params(), e.g. '(-0.05, 0.05)'
            self.params['specNo']        = None         # set inside PyBatman.get_default_params(), e.g. '1-4'
            self.params['paraProc']      = None         # set inside PyBatman.get_default_params(), e.g. 4
            self.params['negThresh']     = -0.5         # probably don't change
            self.params['scaleFac']      = 900000       # probably don't change
            self.params['downSamp']      = 3            # probably don't change
            self.params['hiresFlag']     = 1            # probably don't change
            self.params['randSeed']      = None         # set inside PyBatman.get_default_params(), e.g. 25
            self.params['nItBurnin']     = None         # set inside PyBatman.get_default_params(), e.g. 9000
            self.params['nItPostBurnin'] = None         # set inside PyBatman.get_default_params(), e.g.
            self.params['multFile']      = 2            # probably don't change
            self.params['thinning']      = None         # set inside PyBatman.get_default_params(), e.g. 5
            self.params['cfeFlag']       = 0            # probably don't change
            self.params['nItRerun']      = 5000         # unused
            self.params['startTemp']     = 1000         # probably don't change
            self.params['specFreq']      = 600          # probably don't change
            self.params['a']             = 0.00001      # probably don't change
            self.params['b']             = 0.000000001  # probably don't change
            self.params['muMean']        = None         # set inside PyBatman.get_default_params(), e.g. 0
            self.params['muVar']         = None         # set inside PyBatman.get_default_params(), e.g. 0.1
            self.params['muVar_prop']    = None         # set inside PyBatman.get_default_params(), e.g. 0.002
            self.params['nuMVar']        = None         # set inside PyBatman.get_default_params(), e.g. 0.0025
            self.params['nuMVarProp']    = None         # set inside PyBatman.get_default_params(), e.g. 0.1
            self.params['tauMean']       = None         # set inside PyBatman.get_default_params(), e.g. -0.05
            self.params['tauPrec']       = None         # set inside PyBatman.get_default_params(), e.g. 2
            self.params['rdelta']        = None         # set inside PyBatman.get_default_params(), e.g. 0.002
            self.params['csFlag']        = None         # set inside PyBatman.get_default_params(), e.g. 0

        else:
            self.params = params

        # load parameter descriptions
        with open('params_desc.json') as f:
            data = json.load(f)
        self.params_desc = {}
        for key, value in data:
            self.params_desc[key] = value

    # returns a new copy each time
    def set(self, key, val):
        copy = self.params.copy()
        copy[key] = val
        return PyBatmanOptions(self.selected, params=copy)

class Database(object):

    def __init__(self):
        self.metabolites = {}

    def add(self, m):
        self.metabolites[m.name] = m

    def activate(self, name):
        self.metabolites[name].active = True

    def deactivate_all(self):
        for name in self.metabolites:
            self.metabolites[name].active = False

    def get_active(self):
        metabolites = [m for m in self.metabolites.values() if m.active]
        return metabolites

    def __repr__(self):
        return 'Database of %d entries' % len(self.metabolites)

class Metabolite(object):

    def __init__(self, name, ppm_range, active=False):
        self.name = name
        self.ppm_range = ppm_range
        self.active = active
        self.multiplets = []

    def add(self, ppm, couple_code, j_constant, rel_intensity):
        # doesn't have to be integers, e.g. 1,1
        couple_code = str(couple_code)
        j_constant = str(j_constant)
        m = Multiplet(self, ppm, couple_code, j_constant, rel_intensity)
        self.multiplets.append(m)
        return self

    def __repr__(self):
        return 'name=%s, active=%s' % (self.name, self.active)

class Multiplet(object):

    def __init__(self, parent, ppm, couple_code, j_constant, rel_intensity):
        self.parent = parent
        self.ppm = ppm
        self.couple_code = couple_code
        self.j_constant = j_constant
        self.rel_intensity = rel_intensity

    def __repr__(self):
        output = '(ppm=%s, couple_code=%s, j_constant=%s, rel_intensity=%s)' % (self.ppm,
                    self.couple_code, self.j_constant, self.rel_intensity)
        return output

class BatmanOutput(object):

    def __init__(self, bm):

        self.bm = bm
        beta_r = bm[bm.names.index('beta')]
        self.beta_df = pandas2ri.ri2py(beta_r)