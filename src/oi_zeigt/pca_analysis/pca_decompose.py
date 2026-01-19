from __future__ import print_function
import pickle
import random
import uuid
import copy

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.decomposition import SparsePCA
from sklearn.decomposition import FastICA
from sklearn import preprocessing
import numpy.ma as ma
from scipy import stats

import pyclass
import pgutils
from sicparse import OptionParser

from .pca_utilities import (
    detect_if_line_is_present,
    create_pickle_file_name,
    prepare_spectrum,
    apply_pca_exclude_range,
    create_index,
    get_line_windows,
    create_output_folders,
    create_last_tag_file_name,
    add_decompose_options,
    add_common_options,
    load_config_file,
)


from .pca_errors import InvalidOption, MissingMandatoryOption, NoDataFound


class InputSpectra(object):

    def __init__(self, spectra=None, scramble=None):
        self.scramble = scramble
        if spectra is None:
            self.spectra = None
        else:
            self.spectra = np.asarray([spectra])
            if self.scramble is not None:
                self.scramble_spectra()

    def __len__(self):
        if self.spectra is None:
            return 0
        else:
            return len(self.spectra)

    def scramble_spectra(self):
        spectra = []
        sign = 1
        for spec in self.spectra:
            spectra += [sign * random.random() * spec]
            sign = -1 * sign
        self.spectra = np.asarray(spectra)

    def add_spectrum(self, spectrum):
        # z = np.abs(stats.zscore(spectrum))
        # spectrum[np.where(z > 3)] = 0
        if self.spectra is None:
            self.spectra = np.asarray(spectrum)
        else:
            self.spectra = np.vstack([self.spectra, np.asarray(spectrum)])

    @property
    def normalized_array(self):
        if self.scramble is not None:
            self.scramble_spectra()
        # scaler = preprocessing.StandardScaler().fit(self.spectra)
        # sky_diff_array_normalized = scaler.transform(self.spectra)
        sky_diff_array_normalized = self.spectra
        return sky_diff_array_normalized

    @property
    def mean(self):
        return self.spectra.mean(0)


class PCADecomposition(object):
    def __init__(self, config=None, input_spectra=None, type="PCA"):
        self.input_spectra = input_spectra
        self.config = dict(config)
        self.global_structure = []
        self.number_components_used = self.config["number_components"]
        if input_spectra is not None:
            if type == "PCA":
                self.compute_PCA()
            if type == "SparsePCA":
                self.computeSparsePCA()
            if type == "ICA":
                self.compute_ICA()
            self.smooth_pca_components()
            self.evaluate_global_structure()
            # self.scale_to_running_noise()

    @property
    def telescope(self):
        return self.config["telescope"]

    @property
    def mission_id(self):
        return self.config["mission_id"]

    @property
    def group_name(self):
        return "{}_{}".format(self.__mission_id, self.__telescope)

    def check_number_of_components(self):
        """Check that enough spectra exist to derive the desired

        number of components
        """
        if self.number_components_used > len(self.input_spectra):
            self.number_components_used = len(self.input_spectra)
            pyclass.message(
                pyclass.seve.i,
                "PCA",
                (
                    "There are not enough spectra for this group: "
                    "{} to derive {} components. Maximum set to {} "
                    "for this component"
                ),
            ).format(
                self.group_name,
                self.config["number_components"],
                self.number_components_used,
            )
            decomposition_output_file = create_pickle_file_name(self.config)
            empty_decomposition = pd.DataFrame()
            with open(decomposition_output_file, "wb") as filehandler:
                pickle.dump(empty_decomposition, filehandler)

    def evaluate_global_structure(self):

        for component in self.pca.components_:
            if self.config["rolling_noise_window"]:
                running_rms_noise = np.std(
                    rolling_window(component, self.config["rolling_noise_window"]), 1
                )
                mean_running_rms_noise = running_rms_noise.mean()
                global_noise = component.std()
                self.global_structure += [global_noise / mean_running_rms_noise]
            else:
                self.global_structure += [np.nan]

    def scale_to_running_noise(self):

        for i, component in enumerate(self.pca.components_):
            if self.config["rolling_noise_window"]:
                running_rms_noise = np.std(
                    rolling_window(component, self.config["rolling_noise_window"]), 1
                )
            else:
                running_rms_noise = np.std(
                    rolling_window(component, self.config["rolling_noise_window"]), 1
                )
            global_noise = component.std()
            mean_running_rms_noise = running_rms_noise.mean()
            self.pca.components_[i] = component / global_noise
        self.pca_comp = np.vstack(
            [self.input_spectra.normalized_array.mean(0), self.pca.components_]
        )

    def smooth_pca_components(self):

        if self.config["smoothing_kernel_size"]:
            kernel_size = self.config["smoothing_kernel_size"]
            kernel = np.ones(kernel_size) / kernel_size
            smoothed_components = []
            for component in self.pca.components_:
                component_ = copy.deepcopy(component)
                smoothed_components += [np.convolve(component_, kernel, mode="same")]
            self.pca.components_ = smoothed_components
            self.pca_comp = np.vstack(
                [self.input_spectra.normalized_array.mean(0), self.pca.components_]
            )

    def compute_PCA(self):
        """Create the PCA decomposition of the spectra"""
        self.check_number_of_components()
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            (
                "Performing PCA decomposition of "
                "spectra from {} Data for group {}_{}"
            ).format(self.config["pca_source"], self.telescope, self.mission_id),
        )

        pca = PCA(self.number_components_used - 1)
        pca.fit(self.input_spectra.normalized_array)
        self.pca = pca
        self.pca_comp = np.vstack(
            [self.input_spectra.normalized_array.mean(0), pca.components_]
        )

    def computeSparsePCA(self):
        self.check_number_of_components()
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            (
                "Performing PCA decomposition of "
                "spectra from {} Data for group {}_{}"
            ).format(self.config["pca_source"], self.telescope, self.mission_id),
        )

        pca = SparsePCA(self.number_components_used - 1)
        pca.fit(self.input_spectra.normalized_array)
        self.pca = pca
        self.pca.explained_variance_ratio_ = np.arange(len(pca.components_))
        self.pca_comp = np.vstack(
            [self.input_spectra.normalized_array.mean(0), pca.components_]
        )

    def compute_ICA(self):
        """Create the PCA decomposition of the spectra"""
        self.check_number_of_components()
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            (
                "Performing PCA decomposition of "
                "spectra from {} Data for group {}_{}"
            ).format(self.config["pca_source"], self.telescope, self.mission_id),
        )
        pca = FastICA(self.number_components_used - 1)
        pca.fit(self.input_spectra.normalized_array)
        self.pca = pca
        self.pca_comp = np.vstack(
            [self.input_spectra.normalized_array.mean(), pca.components_]
        )


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def create_components(name, group, config):
    config = dict(config)
    telescope = str(name[0]).strip()
    mission_id = str(name[1]).strip()

    config["telescope"] = telescope
    config["mission_id"] = mission_id

    input_spectra = InputSpectra()

    pyclass.message(
        pyclass.seve.i,
        "PCA",
        "Using spectra from source {} to derive pca components".format(
            config["pca_source"]
        ),
    )
    dropped = []
    kept = []
    if not config["add_sky_diff"]:
        numbers = group.loc[group.source == config["pca_source"], "number"]
    else:
        numbers = group.loc[
            (group.source == config["pca_source"]) | (group.source == "SKY-DIFF"),
            "number",
        ]
    total_number_of_spectra = numbers.count()
    if numbers.count() == 0:
        pyclass.message(
            pyclass.seve.e,
            "PCA",
            ("No Data found pixel: {} " "mission_id: {}-> Skipping").format(
                telescope, mission_id
            ),
        )
        raise NoDataFound

    for number in numbers:

        source = group.loc[group.number == number, "source"]
        source = source.to_string(index=False)

        # The line detection tries to check the input components_
        # for the presence of astronomical lines. If one is found then
        # the component is dropped. STILL EXPERIMENTAL
        if source not in ["SKY-DIFF", "SKYCHOPDIFF"]:
            try:
                line_present = detect_if_line_is_present(
                    number=number,
                    config=config,
                )
            except pgutils.PygildasError:
                continue
            if line_present:
                dropped += [number]
                continue
            else:
                kept += [number]
        else:
            kept += [number]
        try:
            spectrum, outlier_mask = prepare_spectrum(number=number, config=config)
            spectrum = apply_pca_exclude_range(spectrum, config)
        except pgutils.PygildasError:
            continue
        spectrum_rms = copy.deepcopy(pyclass.gdict.sigma.__sicdata__)
        # spectrum = spectrum * 1000  # spectrum_rms
        input_spectra.add_spectrum(spectrum)
        # sky_diff_array += [spectrum]
    print(
        (
            "group {}:{} -- number_of_spectra: {}; "
            "dropped {} spectra; kept {} spectra"
        ).format(
            telescope, mission_id, total_number_of_spectra, len(dropped), len(kept)
        )
    )
    if len(input_spectra) == 0:
        print("no data")
        decomposition_output_file = create_pickle_file_name(config)
        empty_decomposition = pd.DataFrame()
        with open(decomposition_output_file, "wb") as filehandler:
            pickle.dump(empty_decomposition, filehandler)
        raise NoDataFound
    pyclass.comm("sic message class s+i")

    decomposition = PCADecomposition(
        config=config,
        input_spectra=input_spectra,
    )

    # Save the derived components to a pickle file
    decomposition_output_file = create_pickle_file_name(config)

    with open(decomposition_output_file, "wb") as filehandler:
        pickle.dump(decomposition, filehandler)
    print("write pickle file: {}".format(decomposition_output_file))

    with open(create_last_tag_file_name(), "w") as _file:
        _file.write(config["tag"])


def command_line_arguments():
    parser = OptionParser()

    parser = add_decompose_options(parser)
    parser = add_common_options(parser)

    try:
        (options, args) = parser.parse_args()
    except KeyError:
        pyclass.message(pyclass.seve.e, "PCA", "Invalid option")
        pyclass.sicerror()
        raise InvalidOption

    # Start with command line options
    config = options.__dict__

    if options.config_file:
        pyclass.message(
            pyclass.seve.i,
            "PCA",
            "Read configuration from {}".format(options.config_file),
        )
        config_from_file = load_config_file(options.config_file)
        if not config_from_file:
            raise InvalidOption
        if options.debug:
            print(config_from_file)

        ## Update with config file settings, but only if not set in command line
        if "common" in config_from_file:
            for key, value in config_from_file["common"].items():
        #        if key not in config:
                 config[key] = value
        #        else:
        #            pyclass.message(
        #                pyclass.seve.i,
        #                "PCA",
        #                "Command line option takes precedence over config file for: {}".format(
        #                    key
        #                ),
        #            )
#
        if "decompose" in config_from_file:
            for key, value in config_from_file["decompose"].items():
        #        if key not in config:
                 config[key] = value
        #        else:
        #            pyclass.message(
        #                pyclass.seve.i,
        #                "PCA",
        #                "Command line option takes precedence over config file for: {}".format(
        #                    key
        #                ),
        #            )


    if options.input_file is None:
        pyclass.message(
            pyclass.seve.e, "PCA", "The input file name is needed: /input_file FILENAME"
        )
        pyclass.sicerror()
        raise MissingMandatoryOption

    config["number_components"] = int(config["number_components"]) + 1
    line_windows = get_line_windows(config)
    config["line_windows"] = line_windows
    return config


def pca_decompose(config=None):

    if config is None:
        try:
            config = command_line_arguments()
        except InvalidOption:
            return
        except MissingMandatoryOption:
            return

    if not pyclass.gotgdict():
        pyclass.get(verbose=False)

    # Setup parameters

    create_output_folders()

    df = create_index(config)

    # The pca reduction has to tbe done per frontend, because each of the
    # frontends has a different response

    telescope_groups = df.groupby(["telescope", "mission_id"])
    if config["tag"] is None:
        config["tag"] = str(uuid.uuid1())

    for name, group in telescope_groups:
        telescope = str(name[0]).strip()
        mission_id = str(name[1]).strip()
        config["telescope"] = telescope
        config["mission_id"] = mission_id
        try:
            pyclass.message(
                pyclass.seve.i, "PCA", "Creating Components for group {}".format(name)
            )
            create_components(name=name, group=group, config=config)
        except NoDataFound:
            pyclass.message(pyclass.seve.i, "PCA", "Empty group {}".format(name))
            decomposition_output_file = create_pickle_file_name(config)
            empty_decomposition = pd.DataFrame()
            with open(decomposition_output_file, "wb") as filehandler:
                pickle.dump(empty_decomposition, filehandler)
            continue
