import datetime
import matplotlib
import numpy as np
from scipy import signal
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
import scipy.interpolate as interpolate

import readline
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

class PyBatman(object):

    def __init__(self, input_dirs, background_dirs, pattern, working_dir, db, verbose=False):

        self.db = db
        self.input_dirs = input_dirs
        self.background_dirs = background_dirs
        self.pattern = pattern
        self.working_dir = working_dir
        self.verbose = verbose

        importr('batman')
        pandas2ri.activate() # to convert between R <-> pandas dataframes

        # see http://stackoverflow.com/questions/3925096/how-to-get-only-the-last-part-of-a-path-in-python
        self.labels = []
        for input_dir in self.input_dirs:
            last = os.path.basename(os.path.normpath(input_dir))
            self.labels.append(last)

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

        # extract spectra and background
        self.spectra = self._load_data(self.input_dirs, self.pattern)
        self.background = self._load_data(self.background_dirs, self.pattern)

        # combine spectra to have the same ppm scales
        self.spectra_data, self.mean_bg = self._combine_spectra(self.spectra, self.background)
        self.spectra_file = os.path.join(self.batman_input, 'NMRdata_temp.txt')
        self._write_batman_input(self.spectra_file, self.spectra_data, self.labels)

    def plot_spectra(self, name):

        meta_ranges = self.db.metabolites[name].ppm_range

        for lower, upper in meta_ranges:

            ppm = self.spectra_data[:, 0]
            idx = (ppm > lower) & (ppm < upper)
            n_row, n_col = self.spectra_data.shape

            plt.figure()
            x = ppm[idx]
            for i in range(1, n_col):
                intensity = self.spectra_data[:, i]
                y = intensity[idx]
                plt.plot(x, y)
            plt.gca().invert_xaxis()
            plt.xlabel('ppm')
            plt.ylabel('intensity')
            plt.title('%s at (%.4f-%.4f)' % (name, lower, upper))
            plt.show()

    def plotly_spectra(self):
        data = []
        for spec in self.spectra:
            trace = go.Scatter(
                x = spec.ppm,
                y = spec.intensity
            )
            data.append(trace)
        plotly.offline.iplot(data)

    def plotly_data(self):
        data = []
        x = self.spectra_data[:, 0]
        n_row, n_col = self.spectra_data.shape
        for i in range(1, n_col):
            y = self.spectra_data[:, i]
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
        spec_no = len(self.spectra)
        para_proc = multiprocessing.cpu_count()

        # set common parameters
        options = options.set('ppmRange', ranges_str)             \
                             .set('specNo', '1-%d' % spec_no)     \
                             .set('paraProc', para_proc)          \
                             .set('nItBurnin', 19000)             \
                             .set('nItPostBurnin', 1000)          \
                             .set('thinning', 5)                  \
                             .set('tauMean', -0.01)               \
                             .set('tauPrec', 2)                   \
                             .set('rdelta', 0.005)                \
                             .set('csFlag', 0)

        # special parameters for TSP since it's so different from the rest
        if is_tsp:
            default = options.set('muMean', 2)                    \
                                 .set('muVar', 0.1)               \
                                 .set('muVar_prop', 0.002)        \
                                 .set('nuMVar', 0)                \
                                 .set('nuMVarProp', 0.1)
        else:
            default = options.set('muMean', 0)                     \
                                 .set('muVar', 0.01)               \
                                 .set('muVar_prop', 0.0002)        \
                                 .set('nuMVar', 0.0025)            \
                                 .set('nuMVarProp', 0.01)


        return default

    def run(self, options, parallel=False, seed=None, plot=True, verbose=False):
        self.verbose = verbose
        if not parallel:
            options = options.set('paraProc', 1)

        if seed is not None:
            random_seed = seed
        else:
            random_seed = randint(0, 1e6)
        options = options.set('randSeed', random_seed)

        # write out batman input files
        selected = options.selected
        metabolites_list_df = self._write_metabolites_list(selected, 'metabolitesList.csv')
        multi_data_user_df = self._write_multiplet_data(selected, 'multi_data_user.csv')
        chem_shift_per_spec_df = self._write_chem_shift(selected, 'chemShiftPerSpec.csv')
        options_lines = self._write_batman_options(options, 'batmanOptions.txt')

        # actually runs batman here using rpy2
        batman_r = robjects.r['batman']
        sink_r = robjects.r['sink']
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

            sink_r()
            bm = batman_r(runBATMANDir=self.temp_dir, txtFile=self.spectra_file, figBatmanFit=False)

        else:

            sink_r('/dev/null')
            bm = batman_r(runBATMANDir=self.temp_dir, txtFile=self.spectra_file, figBatmanFit=False)
            sink_r()

        if plot:
            self._plot_batman_r(bm)
        output = PyBatmanOutput(bm, options)
        return output

    def baseline_correct(self, options):

        print 'Doing baseline correction'
        ppm = self.spectra_data[:, 0]
        n_row, n_col = self.spectra_data.shape

        for i in range(1, n_col):

            intensity = self.spectra_data[:, i]
            selected = options.selected
            label = self.labels[i-1]

            for metabolite in selected:
                ranges = metabolite.ppm_range
                if len(ranges) > 0:
                    for lower, upper in ranges:
                        idx = (ppm > lower) & (ppm < upper)
                        selected_intensity = intensity[idx]
                        selected_ppm = ppm[idx]

                        # divide the array into halves and find the mininum indices from each portion
                        left, right = np.array_split(selected_intensity, 2)
                        to_find = [left.min(), right.min()]
                        to_find_idx = np.in1d(selected_intensity, to_find)
                        min_pos = np.where(to_find_idx)[0]

                        # node_list will be the indices of ...
                        first = 0
                        last = len(selected_intensity)-1
                        # node_list = [first, last]
                        node_list = [first]         # the first element
                        node_list.extend(min_pos)   # the min elements from left and right
                        node_list.append(last)      # the last element
                        corrected_intensity = ng.proc_bl.base(selected_intensity, node_list) # piece-wise baseline correction

                        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
                        title = '(%.4f-%.4f)' % (lower, upper)
                        ax1.plot(selected_ppm, selected_intensity, 'b', label='Spectra')
                        ax1.plot(selected_ppm, corrected_intensity, 'g', label='Baseline Corrected')
                        ax1.set_xlabel('ppm')
                        ax1.set_ylabel('intensity')
                        ax1.legend(loc='best')
                        ax1.invert_xaxis()
                        ax1.set_title(title)

                        intensity[idx] = corrected_intensity
                        idx = (ppm > lower-0.05) & (ppm < upper+0.05)
                        title = 'Corrected (%.4f-%.4f)' % (lower-0.05, upper+0.05)
                        ax2.plot(ppm[idx], intensity[idx], 'b')
                        ax2.set_xlabel('ppm')
                        ax2.invert_xaxis()
                        ax2.set_title(title)

                        plt.suptitle('%s %s' % (label, metabolite.name), fontsize=16, y=1.08)
                        plt.tight_layout()
                        plt.show()
                        assert np.array_equal(self.spectra_data[:, i], intensity)

        self._write_batman_input(self.spectra_file, self.spectra_data, self.labels)

    def background_correct(self, options):

        print 'Doing background correction'
        ppm = self.spectra_data[:, 0]
        background = self.mean_bg
        n_row, n_col = self.spectra_data.shape

        for i in range(1, n_col):

            intensity = self.spectra_data[:, i]
            selected = options.selected
            label = self.labels[i-1]

            for metabolite in selected:
                ranges = metabolite.ppm_range
                if len(ranges) > 0:
                    for lower, upper in ranges:
                        idx = (ppm > lower) & (ppm < upper)
                        selected_intensity = intensity[idx]
                        selected_ppm = ppm[idx]
                        selected_background = background[idx]
                        corrected_intensity = selected_intensity - selected_background

                        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
                        title = '(%.4f-%.4f)' % (lower, upper)
                        ax1.plot(selected_ppm, selected_intensity, 'b', label='Spectra')
                        ax1.plot(selected_ppm, selected_background, 'k--', label='Background')
                        ax1.plot(selected_ppm, corrected_intensity, 'g', label='Spectra - Background')
                        ax1.set_xlabel('ppm')
                        ax1.set_ylabel('intensity')
                        ax1.legend(loc='best')
                        ax1.invert_xaxis()
                        ax1.set_title(title)

                        intensity[idx] = corrected_intensity
                        idx = (ppm > lower-0.05) & (ppm < upper+0.05)
                        title = 'Corrected (%.4f-%.4f)' % (lower-0.05, upper+0.05)
                        ax2.plot(ppm[idx], intensity[idx], 'b')
                        ax2.set_xlabel('ppm')
                        ax2.invert_xaxis()
                        ax2.set_title(title)

                        plt.suptitle('%s %s' % (label, metabolite.name), fontsize=16, y=1.08)
                        plt.tight_layout()
                        plt.show()
                        assert np.array_equal(self.spectra_data[:, i], intensity)

        self._write_batman_input(self.spectra_file, self.spectra_data, self.labels)

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

    def _write_batman_input(self, spectra_file, spectra_data, labels):
        columns = ['ppm']
        columns.extend(labels)
        df = pd.DataFrame(spectra_data, columns=columns)
        df.to_csv(self.spectra_file, index=False, sep='\t')
        if self.verbose:
            print 'Spectra written to', spectra_file

    def _write_metabolites_list(self, metabolites, file_name):
        columns = ['Metabolite']
        data = [m.name for m in metabolites]
        df = pd.DataFrame(data, columns=columns)

        out_file = os.path.join(self.batman_input, file_name)
        if self.verbose:
            print 'metabolites list = %s' % out_file
        df.to_csv(out_file, header=False, index=False, sep=',')

        return df

    # create multi_data_user.csv
    def _write_multiplet_data(self, metabolites, file_name):
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

    # create chemShiftPerSpec.csv
    def _write_chem_shift(self, metabolites, file_name):
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

    # create batmanOptions.txt
    def _write_batman_options(self, options, file_name):
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

        sink_r = robjects.r['sink']
        if self.verbose:
            sink_r()
            plot_batman_fit_r(bm, showPlot=False)
            plot_batman_fit_stack_r(bm, offset=0.8, placeLegend='topleft', yto=5)
            plot_rel_con_r(bm, showPlot=False)
            plot_meta_fit_r(bm, showPlot=False)
        else:
            sink_r('/dev/null')
            plot_batman_fit_r(bm, showPlot=False)
            plot_batman_fit_stack_r(bm, offset=0.8, placeLegend='topleft', yto=5)
            plot_rel_con_r(bm, showPlot=False)
            plot_meta_fit_r(bm, showPlot=False)
            sink_r()

        # TODO: make traceplots etc

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

        return Spectra(x, y)

    def _load_data(self, input_dirs, pattern):
        matching = self._get_matching_paths(input_dirs, pattern)
        spectra_list = [self._load_single_spectra(input_dir) for input_dir in matching]
        return spectra_list

    def _resample(self, xs, ys, ppm):
        new_ys = []
        for x, y in zip(xs, ys):
            new_y = interpolate.interp1d(x, y)(ppm)
            new_ys.append(new_y)
        return new_ys

    def _get_ppms(self, spectra_list):
        xs = [spectra.ppm for spectra in spectra_list]
        minx = min([x[0] for x in xs])
        maxx = max([x[-1] for x in xs])
        ppm = np.linspace(minx, maxx, num=len(xs[0]))
        return ppm

    def _combine_spectra(self, spectra_list, background_list):

        ppm = self._get_ppms(spectra_list + background_list)

        # process the sample spectra first
        xs = [spectra.ppm for spectra in spectra_list]
        ys = [spectra.intensity for spectra in spectra_list]
        new_ys = self._resample(xs, ys, ppm)

        combined = [ppm]
        combined.extend(new_ys)
        combined = np.array(combined).transpose()
        if self.verbose:
            print 'Loaded', combined.shape

        # process the background too if available
        if len(background_list) > 0:

            xs = [spectra.ppm for spectra in background_list]
            ys = [spectra.intensity for spectra in background_list]
            new_ys = self._resample(xs, ys, ppm)

            # get the average of all the background spectra
            all_bg = np.array(new_ys)
            mean_bg = all_bg.mean(axis=0)

        else:
            mean_bg = None

        return combined, mean_bg

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
            self.params['downSamp']      = 1            # probably don't change
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

class Spectra(object):

    def __init__(self, ppm, intensity):
        self.ppm = ppm
        self.intensity = intensity

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

class PyBatmanOutput(object):

    def __init__(self, bm, options):

        self.output_dir = bm[bm.names.index('outputDir')][0]
        self.bm = bm
        self.options = options

        beta = bm[bm.names.index('beta')]
        # beta_sam = bm[bm.names.index('betaSam')]
        specfit = bm[bm.names.index('sFit')]

        self.beta_df = pandas2ri.ri2py(beta)
        # self.beta_sam_df = pandas2ri.ri2py(beta_sam).transpose()

        self.specfit_df = pandas2ri.ri2py(specfit)
        self.ppm = self.specfit_df['ppm'].values
        self.original_spectrum = self.specfit_df['Original.spectrum'].values
        self.metabolites_fit = self.specfit_df['Metabolites.fit'].values
        self.wavelet_fit = self.specfit_df['Wavelet.fit'].values
        self.overall_fit = self.specfit_df['Overall.fit'].values

    def plot_fit(self):

        for metabolite in self.options.selected:

            for lower, upper in metabolite.ppm_range:

                idx = (self.ppm > lower) & (self.ppm < upper)
                ppm = self.ppm[idx]
                original_spectrum = self.original_spectrum[idx]
                metabolites_fit = self.metabolites_fit[idx]
                wavelet_fit = self.wavelet_fit[idx]

                plt.figure()
                plt.plot(ppm, original_spectrum, color='blue', label='Original Spectra')
                plt.plot(ppm, metabolites_fit, color='green', label='Metabolites Fit')
                plt.plot(ppm, wavelet_fit, color='red', label='Wavelet Fit')
                # plt.plot(self.ppm, self.overall_fit, '--', color='black', label='Combined Fit')
                plt.legend(loc='best')
                plt.gca().invert_xaxis()
                title = 'Fit Results -- %s (%.4f-%.4f)' % (metabolite.name, lower, upper)
                plt.title('%s' % title)
                plt.show()

    # def plot_beta_sam(self):
    #     self.beta_sam_df.boxplot()

    def rmse(self):
        error = np.sqrt((self.error_vect() ** 2).mean())
        return error

    def error_vect(self):
        return self.original_spectrum - self.metabolites_fit