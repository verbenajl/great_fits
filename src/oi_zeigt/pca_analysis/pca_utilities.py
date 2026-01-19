from __future__ import print_function
import math
import copy
import os
import pickle
import yaml
import subprocess

import pyclass
import pgutils

import pandas as pd
import numpy as np

import cv2 as cv
import numpy.ma as ma
from jinja2 import Environment, PackageLoader
from scipy import stats

from .pca_errors import NoScienceSourceFound, MultipleScienceSourcesFound


def jinja_raise(msg):
    raise Exception(msg)


def debug(text):
    return ''


ENV = Environment(
    loader=PackageLoader('kosma_py_lib', 'templates'),
    trim_blocks=True,
    lstrip_blocks=True)
ENV.globals['jinja_raise'] = jinja_raise
ENV.filters['debug'] = debug

NON_SCIENCE_SOURCES = [
    "SKY-DIFF", "SKYCHOPDIFF",
    "TAU_SIG", "TAU_IMG", "TAU_AVG",
    "S-H_OBS", "S-H_FIT", "TSYS",
    "CAL_SIG", "CAL_IMG", "TREC (SSB)",
    "TCOLD_SIG", "TCOLD_IMG", "THOT_IMG",
    "THOT_SIG", "HOT-COLD"
]

PICKLE_SUBFOLDER = ".pca_pickled_objects"


def create_output_folders():
    if not os.path.isdir(PICKLE_SUBFOLDER):
        subprocess.call(["mkdir", PICKLE_SUBFOLDER])


def create_pickle_file_name(config):
    tag = config["tag"]
    if tag is None:
        with open(create_last_tag_file_name(), "r") as filehandler:
            tag = filehandler.read().strip()
    try:
        tag = tag.decode()
    except AttributeError:
        pass
    return "{}/{}___{}___{}.pkl".format(
        PICKLE_SUBFOLDER, tag, config["telescope"],
        config["mission_id"])


def create_last_tag_file_name():
    return "{}/last_decomposition_tag".format(
        PICKLE_SUBFOLDER)


def derive_channel_from_velocity(reference, central_velocity,
                                 channel_velocity, velo_res):
    channel = reference - math.floor(
        ((central_velocity - float(channel_velocity)) / velo_res)
    )
    return int(channel)


def detect_if_line_is_present(number, config):
    ''' Checking if a line is present in the input spectra

    '''
    # noise_cutoff turns this feature on
    if config["noise_cutoff"] is False:
        return False

    flat_line_windows = np.asarray(config["line_windows"]).flatten()
    pyclass.comm("sic message class s-i")
    pyclass.comm("get {}".format(number.decode()))
    pyclass.comm("sic message class s+i")
    if config["look_ahead_resolution"]:
        pyclass.comm(
            "resample * * * {} v".format(config["look_ahead_resolution"])
        )
    if config["extract_range"]:
        pyclass.comm("extract {} {} v".format(
            config["extract_range"][0], config["extract_range"][1]))

    pyclass.comm("set unit v")
    if config["line_windows"]:
        windows = [str(entry) for entry in flat_line_windows]
        try:
            pyclass.comm("set window {}".format(" ".join(windows)))
            print(windows)
            pyclass.comm("bas 3")
        except pgutils.PygildasError:
            raise

    reference = copy.deepcopy(
        pyclass.gdict.reference.__sicdata__)
    velocity = copy.deepcopy(
        pyclass.gdict.velocity.__sicdata__)
    velo_res = copy.deepcopy(
        pyclass.gdict.velo_step.__sicdata__)
    # channels = copy.deepcopy(
    #     pyclass.gdict.channels.__sicdata__)

    # Derive the channels that correspond to the line
    # velocities
    window_channels = [
        derive_channel_from_velocity(
            reference=reference,
            central_velocity=velocity,
            velo_res=velo_res,
            channel_velocity=window_entry
        ) for window_entry in flat_line_windows
    ]
    # TBD

    window_channels = sorted(window_channels)
    line_free_regions = [0] + window_channels + [1]
    noise_array = []
    for i in range(0, len(line_free_regions), 2):
        noise_array = np.append(
            noise_array,
            copy.deepcopy(
                pyclass.gdict.ry.__sicdata__[
                    line_free_regions[i]:line_free_regions[
                        i + 1
                    ]
                ]
            )
        )
    this_noise_measure = np.std(noise_array)
    line_array = []
    for i in range(0, len(window_channels), 2):
        line_array = np.append(
            line_array,
            copy.deepcopy(
                pyclass.gdict.ry.__sicdata__[
                    window_channels[i]:window_channels[i + 1]
                ]
            )
        )
    if np.count_nonzero(
        line_array > float(config["noise_cutoff"]) * this_noise_measure
    ) > config["channel_line_width"]:
        return True
    else:
        return False


def baseline_spectrum(config, baseline_order=0, associated_array=False):
    # TODO: we want that all spectra have the LINE associate array to be able
    # to do a general average of them. So some have a dummy array, let's deal
    # with that here, so that they are properly windowed.
    if associated_array:
        pyclass.comm("set window /associated")
        pyclass.comm("sic message class s-i")
        pyclass.comm("bas {}".format(baseline_order))
        pyclass.comm("sic message class s+i")
        return True

    elif config["line_windows"]:
        flat_line_windows = np.asarray(config["line_windows"]).flatten()
        windows = [str(entry) for entry in flat_line_windows]
        pyclass.comm("set window {}".format(" ".join(windows)))
        pyclass.comm("sic message class s-i")
        pyclass.comm("bas {}".format(baseline_order))
        pyclass.comm("sic message class s+i")
        return True

    else:
        return False


def prepare_spectrum(number, config, replace_ry=None,
                     window_mask=None, baseline_order=0,
                     z_score_blank=None):
    pyclass.comm("sic message class s-i")
    try:
        number = number.decode()
    except AttributeError:
        pass
    pyclass.comm("get {}".format(number))
    pyclass.comm("sic message class s+i")

    if config["final_resolution"]:
        pyclass.comm(
            "resample {} {} * {} v".format(config["nx_user"],
                                           config["xref_user"],
                                           config["final_resolution"]))
    if config["extract_range"]:
        pyclass.comm("extract {} {} v".format(
            config["extract_range"][0], config["extract_range"][1]))
    if replace_ry is not None:
        pyclass.gdict.ry = copy.deepcopy(replace_ry)
    associated_array = False
    # TODO: all spectra have the LINE array now, let's properly
    # baseline them.
    if window_mask is not None:
        try:
            pyclass.comm("del /var window_mask")
        except pgutils.PygildasError:
            pass
        pyclass.comm("def real window_mask /like ry")
        try:
            pyclass.gdict.window_mask = copy.deepcopy(window_mask)
        except ValueError:
            print(window_mask, window_mask[0])
            pyclass.gdict.window_mask = copy.deepcopy(window_mask[0])
        try:
            pyclass.comm("associate line /delete")
        except pgutils.PygildasError:
            pass
        pyclass.comm("associate line window_mask")
        associated_array = True
        baseline_spectrum(config, associated_array=associated_array,
                      baseline_order=baseline_order)
    # deal with the case where no window mask was detected.
    elif window_mask is None and config["line_window"]:
        window_mask= np.zeros_like(pyclass.gdict.ry.__sicdata__)
        try:
            pyclass.comm("del /var window_mask")
        except pgutils.PygildasError:
            pass
        pyclass.comm("def real window_mask /like ry")
        try:
            pyclass.gdict.window_mask = copy.deepcopy(window_mask)
        except ValueError:
            print(window_mask, window_mask[0])
            pyclass.gdict.window_mask = copy.deepcopy(window_mask[0])
        try:
            pyclass.comm("associate line /delete")
        except pgutils.PygildasError:
            pass
        pyclass.comm("associate line window_mask")
        baseline_spectrum(config, associated_array=False,
                          baseline_order=baseline_order)
    #pyclass.comm("associate line window_mask")
    #associated_array = True
    #baseline_spectrum(config, associated_array=associated_array,
    #                  baseline_order=baseline_order)
    spectrum = copy.deepcopy(pyclass.gdict.ry.__sicdata__)
    outlier_mask = None
    if config["z_score_cutoff"]:
        bad_value = float(
            copy.deepcopy(
                pyclass.gdict.r.head.spe.bad.__sicdata__
            )
        )
        z_score_spectrum = copy.deepcopy(spectrum)
        z_score_spectrum[
            np.where(
                np.isclose(z_score_spectrum, bad_value, atol=1e6))] = np.nan
        z = np.abs(stats.zscore(z_score_spectrum, nan_policy="omit"))
        outlier_channels = np.where(z >= int(config["z_score_cutoff"]))
        outlier_mask = np.where(z >= int(config["z_score_cutoff"]), 1, 0)

        spectrum[outlier_channels] = bad_value
        pyclass.gdict.ry = copy.deepcopy(spectrum)

    return spectrum, outlier_mask


def get_and_check_science_sources(group, scan, subscan, config=None):
    # Check for science sources
    non_science_sources = NON_SCIENCE_SOURCES

    if config:
        if (config["pca_source"] != "SKYDIFF") or (config["pca_source"] != "SKYCHOPDIFF"):
            # + [config["pca_source"]]
            non_science_sources = NON_SCIENCE_SOURCES

    science_sources = [source for source in
                       group.source.unique() if source not in
                       non_science_sources]

    if len(science_sources) > 1:
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            "Multiple science sources found for scan {}: {}".format(
                scan,
                subscan
            )
        )
        raise MultipleScienceSourcesFound

    elif len(science_sources) == 0:
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            "No science sources found for scan {}, subscan {}: {}".format(
                scan,
                subscan,
                science_sources
            )
        )
        raise NoScienceSourceFound
    else:
        science_source = science_sources[0]
    return science_source


def refill_pca_exclude_range(spectrum, original_spectrum, config):
    if config["pca_exclude_range"]:
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            ("Filling original data back into spectrum for "
             "the excluded range {}").format(config["pca_exclude_range"])
        )

        reference = copy.deepcopy(
            pyclass.gdict.reference.__sicdata__)
        velocity = copy.deepcopy(
            pyclass.gdict.velocity.__sicdata__)
        velo_res = copy.deepcopy(
            pyclass.gdict.velo_step.__sicdata__)
        channels = [derive_channel_from_velocity(
            reference=reference,
            central_velocity=velocity,
            velo_res=velo_res,
            channel_velocity=exclude_velo
        ) for exclude_velo in config["pca_exclude_range"]
        ]
        channels = sorted(channels)

        spectrum = np.concatenate(
            (spectrum[0:channels[0] - 1],
             original_spectrum[
                 channels[0] - 1:channels[1] - 1
            ],
                spectrum[channels[0] - 1:],
            )
        )
    return spectrum


def get_exclude_channels(config):
    channels = []
    if config["pca_exclude_range"]:
        reference = copy.deepcopy(
            pyclass.gdict.reference.__sicdata__)
        velocity = copy.deepcopy(
            pyclass.gdict.velocity.__sicdata__)
        velo_res = copy.deepcopy(
            pyclass.gdict.velo_step.__sicdata__)

        # Derive the channels that correspond to the line
        # velocities
        channels = [derive_channel_from_velocity(
            reference=reference,
            central_velocity=velocity,
            velo_res=velo_res,
            channel_velocity=exclude_velo
        ) for exclude_velo in config["pca_exclude_range"]
        ]
        channels = sorted(channels)
    return channels


def apply_pca_exclude_range(spectrum, config):
    if config["pca_exclude_range"]:
        channels = get_exclude_channels(config)
        spectrum = ma.concatenate(
            (spectrum[0:channels[0] - 1], spectrum[channels[1] - 1:]))
    return spectrum


def blank_line_fit(spectrum, config):
    fit_spec = copy.deepcopy(spectrum)
    fit_spec_tmp = copy.deepcopy(spectrum)

    if ((config["blank_line_fit"]
         and config["line_window"] and config["blank_line_cutoff"])):

        sigma = copy.deepcopy((pyclass.gdict.sigma.__sicdata__))
        # noise = np.random.normal(0, sigma, len(fit_spec))
        test_spike = np.where(
            fit_spec > float(config["blank_line_cutoff"]) * sigma)
        if np.any(test_spike):
            # TODO only do this in the velocity window....
            fit_spec_tmp[:] = False
            fit_spec_tmp[test_spike] = True
            reference = copy.deepcopy(
                pyclass.gdict.reference.__sicdata__)
            velocity = copy.deepcopy(
                pyclass.gdict.velocity.__sicdata__)
            velo_res = copy.deepcopy(
                pyclass.gdict.velo_step.__sicdata__)
            # channels = copy.deepcopy(
            #     pyclass.gdict.channels.__sicdata__)

            # Derive the channels that correspond to the line
            # velocities
            channels = [derive_channel_from_velocity(
                reference=reference,
                central_velocity=velocity,
                velo_res=velo_res,
                channel_velocity=window_entry
            ) for window_entry in config["line_window"]
            ]
            channels = sorted(channels)
            channels = [0] + channels + [-1]
            for i in range(0, len(channels), 2):
                print(channels[i], channels[i + 1])
                fit_spec_tmp[channels[i]:channels[i + 1]] = False
            fit_spec = ma.array(fit_spec, mask=fit_spec_tmp)
    return fit_spec


def get_line_windows(config):
    line_windows = []
    if config["line_window"]:
        if (len(config["line_window"]) == 0 or len(
                config["line_window"]) % 2 != 0):
            pyclass.message(pyclass.seve.e, "PCA", "Please add an even "
                            "list of window parameters")
            pyclass.sicerror()
            return
        for i in range(0, len(config["line_window"]), 2):
            line_windows += [[config["line_window"][i],
                              config["line_window"][i + 1]]]
    else:
        line_windows = None
    return line_windows


def set_line_windows(config):
    if config["line_windows"] is None:
        pyclass.message(pyclass.seve.e, "PCA", "No line windows defined")
    else:
        for window in config["line_windows"]:
            pyclass.comm('set window {} {}'.format(window[0], window[1]))
    return None


def create_index(config):
    redo_index = True
    if config["input_file"]:
        pyclass.comm('file in "{}"'.format(config["input_file"]))
        pyclass.comm('find')
        modification_times_yaml_file = (
            ".{}_modification_times_pca.yaml".format(
                config["input_file"].replace("/", "_"))
        )
        modification_time = os.path.getmtime(config["input_file"])

        if os.path.exists(modification_times_yaml_file):
            previous_modification_times = yaml.safe_load(
                open(modification_times_yaml_file))
        else:
            previous_modification_times = {}
        if previous_modification_times is None:
            previous_modification_times = {}
        if config["input_file"] in previous_modification_times.keys():
            previous_modification_time = previous_modification_times[
                config["input_file"]]
            if (round(previous_modification_time, 0) >= round(
                    modification_time, 0)):
                redo_index = False

        previous_modification_times[config["input_file"]] = modification_time
        with open(modification_times_yaml_file, 'w') as outfile:
            yaml.dump(previous_modification_times, outfile)

    # Get the index of the data
    pca_index_file = ".pca_index_{0}.pkl".format(
        os.path.basename(config["input_file"])
    )
    if not os.path.isfile(pca_index_file):
        redo_index = True
    if redo_index:
        idx_list = []
        idx_list.append([str(item) for item in pyclass.gdict.idx.num])
        idx_list.append([str(item.decode())
                        for item in pyclass.gdict.idx.teles])
        idx_list.append([str(item) for item in pyclass.gdict.idx.scan])
        idx_list.append([str(item) for item in pyclass.gdict.idx.subscan])
        idx_list.append([str(item.decode())
                        for item in pyclass.gdict.idx.sourc])
        idx_list.append([str(item.decode())
                        for item in pyclass.gdict.idx.line])
        idx_list.append([str(item) for item in pyclass.gdict.idx.boff])
        idx_list.append([str(item) for item in pyclass.gdict.idx.loff])
        idx_list.append([str(item) for item in pyclass.gdict.idx.ver])
        idx_list = np.asarray(idx_list)
        print(idx_list)
        df = pd.DataFrame(
            idx_list.T,
            columns=["number", "telescope", "scan", "subscan",
                     "source", "line", "boff", "loff", "version"]
        )

        pyclass.comm("set var user")
        pyclass.comm("import sofia")
        scan_group = df.groupby(["scan"])
        df["mission_id"] = None

        if config["no_sofia_user_section"]:
            df["mission_id"] = "MISSION_ID_DUMMY"
            df["aot_id"] = "AOT_ID_DUMMY"
        else:
            for name, group in scan_group:
                number = group.number.iloc[0]
                try:
                    pyclass.comm("get {}".format(number))
                except pgutils.PygildasError:
                    # TODO: This is a bad workaround. Sometimes some fields of
                    # older SOFIA user sections are emtpy and newer class versions
                    # give an error this circumvents this. Bit since PygildasError
                    # is not very specific this masks all other problems....
                    print("This has no USER Section! Has it? I don't know!")
                    pass

                mission_id = copy.deepcopy(
                    pyclass.gdict.r.user.sofia.mission_id.__sicdata__)

                df.loc[df["scan"] == group.scan.iloc[0], "mission_id"
                       ] = str(mission_id).replace(
                    "b'", ""
                ).replace(
                    "'", ""
                ).strip()
                df["mission_id"] = df["mission_id"].str.strip()
                aot_id = copy.deepcopy(
                    pyclass.gdict.r.user.sofia.aot_id.__sicdata__)
                df.loc[df["scan"] == group.scan.iloc[0], "aot_id"
                       ] = str(aot_id).replace(
                    "b'", ""
                ).replace(
                    "'", ""
                ).strip()

        df.source = df.source.astype(str).str.strip()
        with open(pca_index_file, "wb") as filehandler:
            pickle.dump(df, filehandler)
    else:
        with open(pca_index_file, "rb") as filehandler:
            pyclass.message(
                pyclass.seve.i, "PCA",
                "Loading pickled index"
            )
            df = pickle.load(filehandler)
    return df


def load_config_file(config_file):
    if os.path.exists(config_file):
        import toml
        config = toml.load(config_file)
        if not config:
            pyclass.message(
                pyclass.seve.e, "PCA",
                ("Configuration file `{}` is empty, edit it or "
                 "do not use it as an argument").format(config_file))
        return config
    else:
        pyclass.message(
            pyclass.seve.e, "PCA",
            ("Configuration File does not exist a "
             "template has been created in: {}").format(config_file))
        create_default_config(config_file)
        return False


def create_default_config(config_file):
    template = ENV.get_template("pca_config.toml")
    with open(config_file, "w") as config_file_out:
        config_file_out.write(template.render())


def gather_self_pca(threshold,):
    pass


def find_lines(input_data, kernel_size=51, title=None,
               plot=False, plot_image=None, config=None, cutoff_std=2):
    from matplotlib.figure import Figure

    # Ensure kernel_size is a positive odd integer
    if kernel_size <= 0 or kernel_size % 2 == 0:
        raise ValueError(f"Invalid kernel_size: {kernel_size}. Must be a positive odd integer.")

    if config:
        line_window = get_line_windows(config)
    input_data = copy.deepcopy(input_data)
    try:
        input_data = input_data.filled(0)
    except AttributeError:
        pass
    data = np.asarray(input_data)
    data_8uc = 255 * (data-np.min(data))/np.max((data-np.min(data)))

    if plot_image is not None:
        plot_image = np.asarray(plot_image)

    img = data_8uc.astype(np.uint8)
    img = cv.normalize(img, None, 0, 100, cv.NORM_MINMAX)
    gray = cv.GaussianBlur(img, (kernel_size, kernel_size), 0)

    mean = gray[np.where(gray != 0)].mean()
    std = gray[np.where(gray != 0)].std()

    _, threshold = cv.threshold(
        gray, mean + cutoff_std * std, 1, cv.THRESH_BINARY
    )

    new_mean = 0
    new_std = 0
    new_gray = copy.deepcopy(gray)
    while ((not np.isclose(new_mean, mean).all()) and
           (not np.isclose(new_std, std).all())):
        new_mean = mean
        new_std = std
        new_gray[np.where(threshold == 1)] = 0
        mean = new_gray[np.where(new_gray != 0)].mean()
        std = new_gray[np.where(new_gray != 0)].std()
        _, threshold = cv.threshold(
            gray, mean + cutoff_std * std, 1, cv.THRESH_BINARY)

    contours, hierarchy = cv.findContours(
        threshold,
        cv.RETR_TREE,
        cv.CHAIN_APPROX_NONE)
    if plot_image is not None:
        cv.drawContours(plot_image, contours, -1,
                        (int(data.max()), 0, int(data.max())), thickness=3)
        data = plot_image
    else:
        cv.drawContours(input_data, contours, -1,
                        (int(data.max()), 0, int(data.max())), thickness=3)

    return threshold, data


def baseline(y, deg=None, max_it=None, tol=None):
    """
    Computes the baseline of a given data.
    Iteratively performs a polynomial fitting in the data to detect its
    baseline. At every iteration, the fitting weights on the regions with
    peaks are reduced to identify the baseline only.
    Parameters
    ----------
    y : ndarray
        Data to detect the baseline.
    deg : int (default: 3)
        Degree of the polynomial that will estimate the data baseline. A low
        degree may fail to detect all the baseline present, while a high
        degree may make the data too oscillatory, especially at the edges.
    max_it : int (default: 100)
        Maximum number of iterations to perform.
    tol : float (default: 1e-3)
        Tolerance to use when comparing the difference between the current
        fit coefficients and the ones from the last iteration. The iteration
        procedure will stop when the difference between them is lower than
        *tol*.
    Returns
    -------
    ndarray
        Array with the baseline amplitude for every original point in *y*
    """
    # for not repeating ourselves in `envelope`
    import scipy.linalg as LA
    if deg is None:
        deg = 3
    if max_it is None:
        max_it = 100
    if tol is None:
        tol = 1e-3

    order = deg + 1
    coeffs = np.ones(order)

    # try to avoid numerical issues
    cond = math.pow(abs(y).max(), 1. / order)
    x = np.linspace(0., cond, y.size)
    base = y.copy()

    vander = np.vander(x, order)
    vander_pinv = LA.pinv2(vander)

    for _ in range(max_it):
        coeffs_new = np.dot(vander_pinv, y)

        if LA.norm(coeffs_new - coeffs) / LA.norm(coeffs) < tol:
            break

        coeffs = coeffs_new
        base = np.dot(vander, coeffs)
        y = np.minimum(y, base)

    return base


def add_decompose_options(parser):
    parser.add_option("-n", "--n_components",
                      dest="number_components", nargs=1, default=5)

    # Controls if skydiffs are added to the pca components sample if the
    # pca_source is not the SKY_DIFFS
    parser.add_option("-a", "--add_sky_diff", dest="add_sky_diff",
                      action="store_true", default=False)
    parser.add_option("-j", "--noise_cutoff",
                      dest="noise_cutoff", nargs=1, default=False)
    parser.add_option("--scramble",
                      dest="scramble",
                      action="store_true", default=False)
    return parser


def add_correct_options(parser):
    parser.add_option("-c", "--cutoff", dest="cutoff", nargs=1, default=False)
    parser.add_option("-o", "--output", dest="output_folder",
                      nargs=1, default=False)
    # This option is active if the pca_source is not the SKY_DIFF spectra.
    parser.add_option("-d", "--line_detection",
                      dest="line_detection", action="store_true",
                      default=False)
    parser.add_option("-b", "--blank_line_cutoff",
                      dest="blank_line_cutoff", nargs=1, default=3)

    # If this option is set the data will be resampled
    # to the resolution in km/s

    parser.add_option("-z", "--blank_line_fit",
                      dest="blank_line_fit", action="store_true",
                      default=False)
    parser.add_option("--dump_exp_variance",
                      dest="dump_exp_variance", nargs=1, default=False)
    parser.add_option("--export_components",
                      dest="export_components", action="store_true",
                      default=False)
    parser.add_option("--cut_coefficients", nargs=1,
                      dest="cut_coefficients", default=0)
    parser.add_option("--noise_ratio_cutoff", nargs=1,
                      dest="noise_ratio_cutoff", default=0)
    parser.add_option("--global_noise_ratio_cutoff", nargs=1,
                      dest="global_noise_ratio_cutoff", default=0)
    parser.add_option("--dump_original_spec_array", nargs=1,
                      dest="dump_original_spec_array", default=False)
    parser.add_option("--custom_baseline_fit", nargs=1,
                      dest="custom_baseline_fit", default=False)
    parser.add_option("--line_kernel_size", nargs=1,
                      dest="line_kernel_size", default=False)
    parser.add_option("--component_map", dest="component_map",
                      nargs=1, default=None)
    parser.add_option("--no_subscan_grouping", action="store_true",
                      dest="no_subscan_grouping", default=False)
    return parser


def add_common_options(parser):
    parser.add_option("--limit_to_scan", dest="limit_to_scan",
                      nargs=1, default=None)
    parser.add_option("--limit_to_telescope", dest="limit_to_telescope",
                      nargs=1, default=None)
    parser.add_option("--config_file", dest="config_file",
                      nargs=1, default=None)
    parser.add_option("-i", "--input_file", dest="input_file",
                      nargs=1, default=None)
    parser.add_option("-p", "--pca_source", dest="pca_source",
                      nargs=1, default="SKY-DIFF")
    parser.add_option("-w", "--window", dest="line_window",
                      nargs="+", default=False)
    parser.add_option("-a", "--look_ahead",
                      dest="look_ahead_resolution", nargs=1, default=False)
    # If this option is set the data will be resampled
    # to the resolution in km/s
    parser.add_option("-f", "--final_res",
                      dest="final_resolution", nargs=1, default=False)
    parser.add_option("-x", "--xref",
                      dest="xref_user", nargs=1, default="*")
    parser.add_option("-r", "--nx",
                      dest="nx_user", nargs=1, default="*")

    # If this option is set the data will be cut before further processing
    parser.add_option("-e", "--e_range", dest="extract_range",
                      nargs=2, default=False)
    parser.add_option("--pca_exclude_range",
                      dest="pca_exclude_range",
                      nargs="+", default=False)
    parser.add_option("--tag",
                      dest="tag", default=None)
    parser.add_option("--channel_line_width",
                      dest="channel_line_width", nargs=1, default=1)
    parser.add_option("--create_self_pca",
                      dest="create_self_pca", nargs=1, default=False)
    parser.add_option("--smoothing_kernel_size", nargs=1,
                      dest="smoothing_kernel_size", default=False)
    parser.add_option("--rolling_noise_window", nargs=1,
                      dest="rolling_noise_window", default=False)
    parser.add_option("--decomposition_type", nargs=1,
                      dest="decomposition_type", default="PCA")
    parser.add_option("--do_scale", nargs=1,
                      dest="do_scale", default=True)
    parser.add_option("--z_score_cutoff", nargs=1,
                      dest="z_score_cutoff", default=False)
    parser.add_option("--debug", action="store_true", dest="debug",
                      default=False)
    parser.add_option("--no_sofia_user_section", action="store_true",
                      dest="no_sofia_user_section",
                      default=False)

    return parser
